#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Michel Spectrum Parameter Scan: Reflectivity x H2oAttenuationLengthCoefficient
- 600-POINT FULL GRID SCAN (30×20)
- WITH PMT QUALITY CUT (>= 2 PE per PMT)

import argparse
import glob
import json
import os
import subprocess
import time
import random
import fcntl

import numpy as np
import awkward as ak
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import uproot

# ============================================================================
# Configuration
# ============================================================================

BASE_DIR = "/home/manoja450/G4WithoutLeadSheilding/MODULE2/CUSTOMOPTICALMODULE2/NEXTmodify/G4H2O"
DATA_DIR = os.path.join(BASE_DIR, "data")
MAC_DIR = os.path.join(BASE_DIR, "mac")
BEAMON_FILE = os.path.join(BASE_DIR, "beamOn.dat")
REAL_DATA_PATH = os.path.join(MAC_DIR, "all_histogramsgood.root")  # Updated to use the good file!
GEANT4_EXECUTABLE = os.path.join(BASE_DIR, "build", "G4d2o")

MIN_PE_THRESHOLD = 60.0

# ================================================================
# PMT QUALITY CUT
# ================================================================
MIN_HITS_PER_PMT = 2
N_PMTS = 12

# 600-POINT GRID: 30 × 20 = 600 points
REFLECTIVITY_GRID = np.linspace(0.90, 0.96, 30)
H2OREFL_GRID = np.linspace(0.15, 0.35, 20)

MU_WEIGHT = 2.0

RESULTS_DIR = "results_600"
PLOT_OUT = "parameter_scan_loss_surface_600.png"
RUN_NUMBER_OFFSET = 500
GEANT4_TIMEOUT_SECONDS = 5400  # 1.5 hours

# Add random delay to avoid race conditions
MIN_DELAY = 0.5
MAX_DELAY = 3.0


# ============================================================================
# Grid indexing
# ============================================================================

def grid_point_for_index(index):
    n_refl = len(REFLECTIVITY_GRID)
    n_h2o = len(H2OREFL_GRID)
    total = n_refl * n_h2o
    if not (0 <= index < total):
        raise ValueError(f"index {index} out of range")
    i, j = divmod(index, n_refl)
    return float(REFLECTIVITY_GRID[j]), float(H2OREFL_GRID[i]), i, j


def total_grid_points():
    return len(REFLECTIVITY_GRID) * len(H2OREFL_GRID)


# ============================================================================
# PMT QUALITY MASK 
# ============================================================================

def build_pmt_quality_mask(pmt_num, min_hits_per_pmt=MIN_HITS_PER_PMT, n_pmts=N_PMTS):
    """
    Per-event mask: True only if every PMT (0..n_pmts-1) registered at
    least min_hits_per_pmt hits in that event.
    pmt_num is a jagged awkward array (one variable-length list of PMT
    indices per event, from pmtHits/pmtHits.pmtNum).
    """
    mask = np.empty(len(pmt_num), dtype=bool)
    for i, evt in enumerate(pmt_num):
        evt_np = np.asarray(ak.to_numpy(evt), dtype=np.int64)
        counts_per_pmt = np.bincount(evt_np, minlength=n_pmts)[:n_pmts]
        mask[i] = np.all(counts_per_pmt >= min_hits_per_pmt)
    return mask


# ============================================================================
# WRITE beamOn.dat WITH FILE LOCKING
# ============================================================================

def write_parameters(R, attenuation, run_number):
    """
    Write a clean beamOn.dat file with FILE LOCKING.
    This prevents race conditions when multiple tasks write simultaneously.
    """
    
    content = f"""\
{run_number}          //Run-number
20000          //primaries-to-generate
0               //Visualization:0=no,1=yes
0               //Debug:0=no-debug,1=debugToScreen,2=DebugToFile
0               //RandomSeed:0=use-default-seed,1=use-clock
1               //PhysicsStatus:0=physicsOff,1=physicsOn
1               //NeutronTreatment:0=noNeutronHP,1=HPneutron
3               //PGAType:0=photonPrimaries,1=electronPrimaries,2=cosmicPrimaries,3=michelElectronPrimaries,4=mu-,5=nueOgen,6=N16,7=neutron,8=vD=e,9=vO=e,10=specialelectrons
0               //SideLining:0=teflon,1=photoCathode
0               //BottomPMTs:0=teflon,1=photoCathode
1               //PMTQE:0=StandardQE,1=HighQE
0               //BottomVeto:0=noVeto,1=useVeto
0               //BottomShielding:0=noShielding,1=useShielding
1               //UserSeed
-8              //PMTDiameter(inches),negNumUsesEllipsoidalPMT
10              //TailCatcherThickness(cm)
2               //ShieldThickness(in)
{attenuation:.6f}          //H2oAttenuationLengthCoefficient
{R:.6f}          //ReflectivityOfTyvek
10.0            //SpecialEnergy
./data          //outputDir
"""
    
    # Use file locking to prevent race conditions
    with open(BEAMON_FILE, 'w') as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        f.write(content)
        f.flush()
        fcntl.flock(f, fcntl.LOCK_UN)
    
    print(f"    [Wrote beamOn.dat: run={run_number}, R={R:.4f}, H2o={attenuation:.4f}]")


def get_current_primaries():
    try:
        with open(BEAMON_FILE, "r") as f:
            for line in f:
                if "//primaries-to-generate" in line:
                    return line.split("//")[0].strip()
    except Exception:
        pass
    return "unknown"


# ============================================================================
# GET REAL DATA μ AND σ (DIRECT WEIGHTED MEAN)
# ============================================================================

def get_real_data_mu_sigma():
    """
    Read real data histogram and compute μ and σ directly.
    Uses all_histogramsgood.root (with quality cuts applied).
    """
    with uproot.open(REAL_DATA_PATH) as f:
        h = f["michel_energy"]
        edges = h.axis().edges()
        counts = h.values()
        centers = (edges[:-1] + edges[1:]) / 2

    # Apply threshold (60 PE cut)
    first_bin = np.where(centers >= MIN_PE_THRESHOLD)[0][0]
    filtered_counts = counts[first_bin:]
    filtered_centers = centers[first_bin:]

    # Compute weighted mean and std (NO FITTING!)
    total = filtered_counts.sum()
    mu_data = np.sum(filtered_centers * filtered_counts) / total
    var_data = np.sum(filtered_counts * (filtered_centers - mu_data)**2) / total
    sigma_data = np.sqrt(var_data)
    
    print(f"Real data: mu={mu_data:.2f} sigma={sigma_data:.2f} (from {total:.0f} events)")
    return mu_data, sigma_data

def run_geant4(sim_output_path, run_number):
    """
    Run Geant4 simulation - PASSES RUN NUMBER AS COMMAND-LINE ARGUMENT.
    This avoids race conditions with beamOn.dat.
    """
    if os.path.exists(sim_output_path):
        os.remove(sim_output_path)

    try:
        result = subprocess.run(
            [GEANT4_EXECUTABLE, str(run_number)],
            cwd=BASE_DIR,
            timeout=GEANT4_TIMEOUT_SECONDS,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            print(f"    [Geant4] non-zero exit code {result.returncode}")
            if result.stderr:
                print(f"    [Geant4 stderr] {result.stderr[-2000:]}")
            return False
        if not os.path.exists(sim_output_path):
            print(f"    [Geant4] finished but expected output not found: {sim_output_path}")
            return False
        return True
    except subprocess.TimeoutExpired:
        print("    [Geant4] TIMED OUT")
        return False
    except Exception as e:
        print(f"    [Geant4] failed to launch: {e}")
        return False


# ============================================================================
# LOAD SIMULATION STATISTICS WITH PMT QUALITY CUT
# ============================================================================

def load_simulation_statistics(sim_output_path):
    """
    Load simulation and compute mean/std.
    APPLIES PMT QUALITY CUT (>= 2 PE per PMT) to match the real data.
    """
    try:
        sim_file = uproot.open(sim_output_path)
        tree = sim_file["Sim_Tree"]
        
        # Get numHits (total PE per event)
        num_hits = tree["eventData/numHits"].array().to_numpy()
        
        # Get PMT hit information for quality cut
        pmt_num = tree["pmtHits/pmtHits.pmtNum"].array()
        
        n_events_all = len(num_hits)
        
        # Apply PMT quality cut: every PMT must have >= MIN_HITS_PER_PMT hits
        pmt_quality_mask = build_pmt_quality_mask(pmt_num)
        
        # Apply threshold (60 PE cut) AND PMT quality cut
        pe_mask = num_hits >= MIN_PE_THRESHOLD
        final_mask = pe_mask & pmt_quality_mask
        num_hits_filtered = num_hits[final_mask]
        
        print(f"    Simulation events: {n_events_all} total -> "
              f"{pmt_quality_mask.sum()} pass {MIN_HITS_PER_PMT} PE/PMT -> "
              f"{len(num_hits_filtered)} pass both cuts")
        
        if len(num_hits_filtered) == 0:
            return None, None, 0

        mu_sim = np.mean(num_hits_filtered)
        sigma_sim = np.std(num_hits_filtered)
        return mu_sim, sigma_sim, len(num_hits_filtered)
        
    except Exception as e:
        print(f"    [ROOT read] failed: {e}")
        return None, None, 0


def loss(mu_sim, sigma_sim, mu_data, sigma_data, mu_weight=MU_WEIGHT):
    return float(np.sqrt((mu_weight * (mu_sim - mu_data)) ** 2 + (sigma_sim - sigma_data) ** 2))


# ============================================================================
# Mode: run
# ============================================================================

def do_run(index):
    reflectivity, h2orefl, i, j = grid_point_for_index(index)
    run_number = RUN_NUMBER_OFFSET + index
    sim_output_path = os.path.join(DATA_DIR, f"Sim_D2ODetector{run_number:03d}.root")

    # Add random delay to avoid race conditions
    delay = random.uniform(MIN_DELAY, MAX_DELAY)
    time.sleep(delay)

    print(f"Task {index}: grid[{i}][{j}] -> Reflectivity={reflectivity:.4f}, "
          f"H2oRefl={h2orefl:.4f}, run_number={run_number}")
    print(f"    primaries-to-generate: {get_current_primaries()}")
    print(f"    PMT quality cut: >= {MIN_HITS_PER_PMT} PE on each of {N_PMTS} PMTs")
    print(f"    Random delay: {delay:.2f}s")

    os.makedirs(RESULTS_DIR, exist_ok=True)
    out_path = os.path.join(RESULTS_DIR, f"point_{index}.json")

    mu_data, sigma_data = get_real_data_mu_sigma()

    record = {
        "index": index, "i": i, "j": j, "run_number": run_number,
        "reflectivity": reflectivity, "h2orefl": h2orefl,
        "mu_data": mu_data, "sigma_data": sigma_data,
        "pmt_quality_cut": f">= {MIN_HITS_PER_PMT} PE on each of {N_PMTS} PMTs",
    }

    try:
        # Write clean beamOn.dat (WITH FILE LOCKING)
        write_parameters(reflectivity, h2orefl, run_number)

        # Run the simulation (WITH COMMAND-LINE RUN NUMBER)
        t0 = time.time()
        ok = run_geant4(sim_output_path, run_number)
        elapsed = time.time() - t0
        print(f"    [Geant4] finished in {elapsed:.1f}s, success={ok}")

        if not ok:
            raise RuntimeError("Geant4 run failed")

        # Load simulation statistics (WITH PMT QUALITY CUT)
        mu_sim, sigma_sim, n_events_passed = load_simulation_statistics(sim_output_path)

        if mu_sim is None:
            raise RuntimeError("Simulation produced zero events above threshold")

        # Compute loss
        L = loss(mu_sim, sigma_sim, mu_data, sigma_data)
        record.update({
            "mu_sim": mu_sim,
            "sigma_sim": sigma_sim,
            "n_events_passed": n_events_passed,
            "loss": L,
            "status": "ok"
        })
        print(f"mu_sim={mu_sim:.2f} sigma_sim={sigma_sim:.2f} loss={L:.3f} (n_events={n_events_passed})")

    except Exception as e:
        record.update({"error": str(e), "status": "failed"})
        print(f"FAILED: {e}")

    # Delete ROOT file to save disk space
    if os.path.exists(sim_output_path):
        try:
            os.remove(sim_output_path)
            print(f"    [Deleted trial file: {os.path.basename(sim_output_path)}]")
        except OSError as e:
            print(f"    [Could not delete trial file: {e}]")

    with open(out_path, "w") as f:
        json.dump(record, f, indent=2)
    print(f"Wrote {out_path}")


# ============================================================================
# Mode: collect
# ============================================================================

def do_collect():
    loss_grid = np.full((len(H2OREFL_GRID), len(REFLECTIVITY_GRID)), np.nan)
    mu_data = sigma_data = None
    n_found = n_failed = 0

    for path in glob.glob(os.path.join(RESULTS_DIR, "point_*.json")):
        with open(path) as f:
            rec = json.load(f)
        i, j = rec["i"], rec["j"]
        if rec.get("status") == "ok":
            loss_grid[i, j] = rec["loss"]
            n_found += 1
        else:
            n_failed += 1
            print(f"  point {rec.get('index')} FAILED: {rec.get('error')}")
        if mu_data is None and "mu_data" in rec:
            mu_data, sigma_data = rec["mu_data"], rec["sigma_data"]

    expected = total_grid_points()
    print(f"Loaded {n_found} ok + {n_failed} failed = {n_found + n_failed} / {expected}")

    if np.all(np.isnan(loss_grid)):
        raise RuntimeError("No successful grid points to plot")

    best_flat_idx = np.nanargmin(loss_grid)
    bi, bj = np.unravel_index(best_flat_idx, loss_grid.shape)
    best_refl = REFLECTIVITY_GRID[bj]
    best_h2orefl = H2OREFL_GRID[bi]
    best_loss = loss_grid[bi, bj]

    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.pcolormesh(REFLECTIVITY_GRID, H2OREFL_GRID, loss_grid,
                        shading="nearest", cmap="cividis_r")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r"$L = \sqrt{%.0f(\mu_{sim}-\mu_{data})^2 + (\sigma_{sim}-\sigma_{data})^2}$" % (MU_WEIGHT ** 2))

    ax.plot(best_refl, best_h2orefl, marker="*", markersize=18,
            color="white", markeredgecolor="black", markeredgewidth=1.2, zorder=5)
    ax.annotate(
        f"Best-fit\nReflectivity={best_refl:.3f}\nH2oRefl={best_h2orefl:.3f}\nL={best_loss:.2f}",
        xy=(best_refl, best_h2orefl), xytext=(10, 10), textcoords="offset points",
        fontsize=9, bbox=dict(boxstyle="round", facecolor="white", alpha=0.85, edgecolor="black"),
    )

    ax.set_xlabel("Reflectivity")
    ax.set_ylabel("Water Attenuation Length Coefficient")
    title = f"Michel Parameter Scan: Loss Surface (600 points, {n_found} successful)"
    if mu_data is not None:
        title += f"\n(PMT cut: >= {MIN_HITS_PER_PMT} PE/PMT, μ_data={mu_data:.1f}, σ_data={sigma_data:.1f} PE)"
    ax.set_title(title, fontsize=12)

    plt.tight_layout()
    plt.savefig(PLOT_OUT, dpi=200, bbox_inches="tight")
    print(f"\nSaved plot to: {PLOT_OUT}")
    print(f"Best fit: Reflectivity={best_refl:.4f}, H2oRefl={best_h2orefl:.4f}, Loss={best_loss:.4f}")

    print("\n" + "="*60)
    print("SUMMARY - 600-POINT GRID SCAN")
    print("="*60)
    print(f"PMT quality cut:        >= {MIN_HITS_PER_PMT} PE on each of {N_PMTS} PMTs")
    print(f"Total grid points:        {expected}")
    print(f"Successful:               {n_found}")
    print(f"Failed:                   {n_failed}")
    print(f"Best Reflectivity:        {best_refl:.4f}")
    print(f"Best H2oRefl:             {best_h2orefl:.4f}")
    print(f"Best Loss:                {best_loss:.4f}")
    if mu_data is not None:
        print(f"Data μ:                   {mu_data:.2f} PE")
        print(f"Data σ:                   {sigma_data:.2f} PE")
    print(f"primaries-to-generate:    {get_current_primaries()}")
    print("="*60)


# ============================================================================
# Entry point
# ============================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["run", "collect"], required=True)
    parser.add_argument("--index", type=int, help="Grid index for --mode run")
    args = parser.parse_args()

    if args.mode == "run":
        if args.index is None:
            parser.error("--mode run requires --index")
        do_run(args.index)
    else:
        do_collect()
