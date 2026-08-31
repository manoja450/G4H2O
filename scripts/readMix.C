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
 * root[1] .x readHax.C
 */

   // Adapted from read.C
   
void readMix()
{
    // Load weight histogram into memory
    auto *weight_file = TFile::Open("../xscnData/haxtonhists.root", "read");
    auto *hweights = weight_file->Get<TH2D>("haxtondd2");
    hweights->Scale(1. / hweights->Integral());

    auto *weight_fileD = TFile::Open("../xscnData/fluxWeight.root", "read");
    auto *hweightsD = weight_fileD->Get<TH2D>("fluxW");
    hweightsD->Scale(1. / hweightsD->Integral());

//    auto cw = new TCanvas("cw", "cw", 700, 500);
//    hweights->Draw("ncolz");

    // Open simulation file and fetch tree
    auto *sim_file = TFile::Open("../data/H2O-from-ornl/Sim_D2ODetector1101.root", "read");
    auto *tree = sim_file->Get<TTree>("Sim_Tree");

    auto *sim_fileD = TFile::Open("../data/Sim_D2ODetector113-.2-.941.root", "read");
    auto *treeD = sim_fileD->Get<TTree>("Sim_Tree");

    // Load branch
    simEvent *event{nullptr};
    tree->SetBranchAddress("eventData", &event);

    simEvent *eventD{nullptr};
    treeD->SetBranchAddress("eventData", &eventD);

    auto h_nhits = new TH1D("", "flat", 100, 0, 500);
    auto h_nhits_weighted = new TH1D("h_nhits_weighted", "h_nhits_weighted", 100, 0, 500);

    auto h_nhitsD = new TH1D("", "flat", 100, 0, 500);
    auto h_nhits_weightedD = new TH1D("h_nhits_weightedD", "h_nhits_weightedD", 100, 0, 500);

    h_nhits_weighted->SetLineColor(kRed);
    h_nhits_weightedD->SetLineColor(kBlue);
    

    std::cout << "Total entries: " << tree->GetEntries() << std::endl;
    std::cout << "Total entries: " << treeD->GetEntries() << std::endl;


    // Loop over tree
    for (int i = 0; i <  tree->GetEntries(); i++)
    {
        if (i % 10000 == 0)
        {
            std::cout << "\rProcessing event: " << i;
            std::cout.flush();
        }

        tree->GetEntry(i);
        TVector3 vB(0.1182589088, -0.37975131476, -0.91749864818);
        auto cos_theta = ROOT::Math::VectorUtil::CosTheta(event->direction0, vB);
        auto new_cos_theta = cos_theta;
        auto energy = event->sourceParticleEnergy;

        auto weight_bin = hweights->FindBin(energy, new_cos_theta);
        auto weight = hweights->GetBinContent(weight_bin);
        auto nhits = event->numHits;
        h_nhits->Fill(nhits);

        auto hnw_bin = h_nhits->FindBin(nhits);
        auto hnw_bin_val = h_nhits_weighted->GetBinContent(hnw_bin);
        h_nhits_weighted->SetBinContent(hnw_bin, hnw_bin_val + weight);
    }
    std::cout << std::endl;

    for (int i = 0; i <  treeD->GetEntries(); i++)
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
        h_nhitsD->Fill(nhitsD);

        auto hnw_binD = h_nhitsD->FindBin(nhitsD);
        auto hnw_bin_valD = h_nhits_weightedD->GetBinContent(hnw_binD);
        h_nhits_weightedD->SetBinContent(hnw_binD, hnw_bin_valD + weightD);
    }
    std::cout << std::endl;

    auto c = new TCanvas("c", "c", 700, 500);
//    h_nhits->Draw();
//    h_nhits_weighted->Draw("sames");
    h_nhits_weighted->Scale(627.*.153/h_nhits_weighted->Integral());
    h_nhits_weightedD->Scale(627./h_nhits_weightedD->Integral());


    for (int i = 1; i <= h_nhits_weighted->GetNbinsX(); ++i) {
    h_nhits_weighted->SetBinError(i, 0.0);
}

    for (int i = 1; i <= h_nhits_weightedD->GetNbinsX(); ++i) {
    h_nhits_weightedD->SetBinError(i, 0.0);
}

    h_nhits_weightedD->Draw();
    h_nhits_weighted->Draw("SAME");
    
    h_nhits_weightedD->SetTitle("Signal for vD and vO; Number of PEs; counts per 10 PEs");
}
