#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TString.h>
#include <fstream>
#include <iostream>

// This script produces a color plot ("2D plot") from a .dat file containing the
// efficiency per cutoff (example, 10 Photoelectrons, 97% efficiency).
// The .dat file should be a CSV file, with lines like "10,.97" (no quotes) as
// per the example above.

void heatmap_from_grid(const char* fname="efficiency_by_energy1400.dat", int forceCSV=-1)
{
    // Auto-detect delimiter unless user forces it
    bool useCSV = false;
    if (forceCSV == -1) {
        std::ifstream f(fname);
        if (!f) { Error("heatmap_from_grid","Cannot open %s", fname); return; }
        std::string line; std::getline(f, line);
        useCSV = (line.find(',') != std::string::npos);
    } else useCSV = (forceCSV != 0);

    // Load into TTree: fields x,y are integers in [1..53], [1..200], z is float [0..1]
    auto T = new TTree("T","x y z grid");
    Long64_t nRead = 0;
    if (useCSV) {
        nRead = T->ReadFile(fname, "x/I:y/I:z/F", ',');
    } else {
        nRead = T->ReadFile(fname, "x/I:y/I:z/F");
    }
    if (nRead <= 0) { Error("heatmap_from_grid","Failed to parse %s", fname); return; }

    // Book a 53 x 200 histogram with bin centers exactly at integer x,y
    // (bins: 1..53, 1..200). Edges at 0.5/53.5 and 0.5/200.5 ensure integer centers.
    auto h2 = new TH2F("h2","f(x,y)=z; x index; y index; z", 
                       53, 0.5, 53.5,
                       200, 0.5, 200.5);

    // Fill: one entry per (x,y). We set the bin content directly (no summing).
    int x, y; float z;
    T->SetBranchAddress("x",&x);
    T->SetBranchAddress("y",&y);
    T->SetBranchAddress("z",&z);

    // Track duplicates or out-of-range rows
    int badRange = 0, dups = 0, total = (int)nRead;

    for (Long64_t i=0; i<nRead; ++i) {
        T->GetEntry(i);
        if (x < 1 || x > 53 || y < 1 || y > 200) { ++badRange; continue; }

        // Check if this bin was already filled (non-zero by default). If your z can be 0.0,
        // use a separate mask TH2C or initialize to NaN (see above) and check std::isnan.
        double prev = h2->GetBinContent(x, y);
        if (prev != 0.0) ++dups;

        h2->SetBinContent(x, y, z);
    }

    // Draw
    gStyle->SetOptStat(0);
    // gStyle->SetPalette(kBird); // uncomment your favorite palette
    auto c = new TCanvas("c_heat","Heatmap f(x,y)=z", 900, 750);
    c->SetRightMargin(0.16);
    h2->SetMinimum(0.0); // since z in [0,1], lock z range if you like
    h2->SetMaximum(1.0);
    h2->Draw("COLZ");

    if (badRange || dups) {
        Info("heatmap_from_grid","Loaded %d rows. Out-of-range: %d, duplicates: %d.",
             total, badRange, dups);
    }
}
