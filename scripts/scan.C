#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>

// G4d2o tree data
#include "../simEvent/simEvent.h"
#include "../simEvent/G4d2oGeom.h"

/*
 * Loop over simEvent data tree.
 *
 * Usage:
 * root load_dictionary.C
 * root[1] .x read.C
 */

// Bare bones version of read.C

void scan()
{
    // Load weight histogram into memory
    auto *weight_file = TFile::Open("../xscnData/fluxWeight.root", "read");
    auto *hweights = weight_file->Get<TH2D>("fluxW");
    hweights->Scale(1. / hweights->Integral());

    //    auto cw = new TCanvas("cw", "cw", 700, 500);
    //    hweights->Draw("ncolz");

    // Open simulation file and fetch tree
    auto *sim_file = TFile::Open("../data/Sim_D2ODetector113-.2-.941.root", "read");
    auto *tree = sim_file->Get<TTree>("Sim_Tree");

    // Load branch
    simEvent *event{nullptr};
    tree->SetBranchAddress("eventData", &event);

    auto h_nhits = new TH1D("", "flat", 100, 0, 700);
    auto h_nhits_weighted = new TH1D("weighted", "weighted", 100, 0, 700);

    h_nhits_weighted->SetLineColor(kRed);

    std::cout << "Total entries: " << tree->GetEntries() << std::endl;
    auto nabove = 0;
    // Loop over tree
    auto maxenergy = 0;
    auto thisangle = 0.0;
    for (int i = 0; i < 60; i++)
    {
        for (float j = 0; j < 101; j++)
        {
            auto thisangle = -1.0 + 2.0 * j / 100.0;
            auto thisone = hweights->FindBin(i, thisangle);
            std::cout << "this energy = " << i << " , thisangle= " << thisangle << ", this bin =" << hweights->GetBinContent(thisone) << std::endl;
            if (hweights->GetBinContent(thisone) > 0 && i > maxenergy)
            {
                maxenergy = i;
            }
        }

    }
    std::cout << "maxenergy= " << maxenergy << std::endl;

    for (int i = 0; i < 0; i++)
    {
        if (i % 10000 == 0)
        {
            std::cout << "\rProcessing event: " << i;
            std::cout.flush();
        }

        treeD->GetEntry(i);
        TVector3 vB(0.1182589088, -0.37975131476, -0.91749864818);
        auto cos_theta = ROOT::Math::VectorUtil::CosTheta(event->direction0, vB);
        auto new_cos_theta = cos_theta;

        auto energy = eventD->sourceParticleEnergy;

        auto weight_binD = hweightsD->FindBin(energy, new_cos_theta);
        auto weightD = hweightsD->GetBinContent(weight_binD);
        auto nhitsD = eventD->numHits;

        if (nhits >= 50)
        {
            nabove += 1;
        }
        h_nhits->Fill(nhits);

        auto hnw_bin = h_nhits->FindBin(nhits);
        auto hnw_bin_val = h_nhits_weighted->GetBinContent(hnw_bin);
        h_nhits_weighted->SetBinContent(hnw_bin, hnw_bin_val + weight);
    }
    std::cout << std::endl;

    auto c = new TCanvas("c", "c", 700, 500);
    double nabover = nabove;
    auto ratioofall = nabover / (tree->GetEntries());
    std::cout << "tree->GetEntries()=" << tree->GetEntries() << std::endl;
    std::cout << "nabove=" << nabove << std::endl;
    std::cout << "ratioofall=" << ratioofall << std::endl;

    //    h_nhits->Draw();
    //    h_nhits_weighted->Draw("sames");
    h_nhits_weighted->Draw();
    h_nhits_weighted->SetTitle("Flux weighted signal MC, .2, .941, D2O vD events; PEs; counts");
}
