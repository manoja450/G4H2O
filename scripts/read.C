#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TString.h>
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

 // This scripts allows you to make plots based on information that
 // is stored inside a .root file produced by G4d2o. Do not forget
 // to run the dictionary first, and then, with root open, run this
 // macro

void read()
{
    // Load weight histogram into memory
    auto *weight_file = TFile::Open("../xscnData/fluxWeight.root", "read");
    auto *hweights = weight_file->Get<TH2D>("fluxW");
    hweights->Scale(1. / hweights->Integral());

    // Open simulation file and fetch tree
    std::ofstream clearFile("./efficiency_by_energy.dat", std::ios::out | std::ios::trunc);
    for (int nn = 1; nn < 54; nn++)
    {
        int file_number = 1200 + nn;
        TString filename = Form("../data/Sim_D2ODetector%d.root", file_number);
        std::cout << "Opening: " << filename << std::endl;
        auto *sim_file = TFile::Open(filename, "READ");

        if (!sim_file || sim_file->IsZombie())
        {
            std::cerr << "Could not open file: " << filename << std::endl;
            return;
        }

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
        for (int i = 0; i < tree->GetEntries(); i++)
        {
            if (i % 1000 == 0)
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

        int var_a = nn;
        double var_b = ratioofall;

        std::ofstream outfile;

        outfile.open("./efficiency_by_energy.dat", std::ios::app); // "app" = append mode

        if (!outfile.is_open())
        {
            std::cerr << "Error: could not open file!" << std::endl;
            return;
        }

        // Write variables, comma separated
        outfile << var_a << "," << var_b << std::endl;

        outfile.close();

        //    h_nhits->Draw();
        //    h_nhits_weighted->Draw("sames");
        h_nhits_weighted->Draw();
        h_nhits_weighted->SetTitle("Flux weighted signal MC, .2, .941, D2O vD events; PEs; counts");

        sim_file->Close();
    }
}
