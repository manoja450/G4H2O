#include "G4d2oEventAction.hh"
#include "G4d2oDetectorHit.hh"
#include "G4d2oNeutrinoAlley.hh"
#include "G4d2oSensitiveDetector.hh"
#include "G4d2oDetector.hh"
#include "G4d2oRunAction.hh"
#include "G4d2oMaterialsDefinition.hh"
#include "G4d2oSteppingAction.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"
#include "G4RunManager.hh"
#include "Randomize.hh"

#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include "TClassTable.h"
#include "inputVariables.hh"

G4d2oEventAction::G4d2oEventAction(void)
{
    G4cerr << "\tConstructing G4d2oEventAction...";

    // Initialize pointers
    HCoE = NULL;
    for (G4int i = 0; i < MAX_SEN_DET; i++)
    {
        detHC[i] = NULL;
        sourceHC[i] = NULL;
        targHC[i] = NULL;
    }

    numEvents = 0;

    input = inputVariables::GetIVPointer();

    G4d2oNeutrinoAlley *detCon = (G4d2oNeutrinoAlley *)G4RunManager::GetRunManager()->GetUserDetectorConstruction();
    theDet = (G4d2oDetector *)detCon->GetDetectorPtr();

    G4d2oRunAction *runAct = (G4d2oRunAction *)G4RunManager::GetRunManager()->GetUserRunAction();
    fout = TFile::Open(runAct->GetOutFileName(), "recreate");

    if (!TClassTable::GetDict("simEvent"))
    {
        G4cout << "Error in G4d2oEventAction::G4d2oEventAction() - class simEvent not found\n"
               << G4endl;
        exit(0);
    }
    if (!TClassTable::GetDict("simHit"))
    {
        G4cout << "Error in G4d2oEventAction::G4d2oEventAction() - class simHit not found\n"
               << G4endl;
        exit(0);
    }
    if (!TClassTable::GetDict("simAreaHit"))
    {
        G4cout << "Error in G4d2oEventAction::G4d2oEventAction() - class simAreaHit not found\n"
               << G4endl;
        exit(0);
    }
    theEventData = new simEvent();

    outTree = new TTree("Sim_Tree", "tree");
    outTree->Branch("eventData", "simEvent", &theEventData);

    numPMTs = theDet->GetTotalPMTS();
    numRows = theDet->GetPMTRows();
    numInEachRow = theDet->GetPMTsInRow();

    if (numPMTs != numRows * numInEachRow)
    {
        printf("\nError in G4d2oEventAction::G4d2oEventAction() - %d x %d != %d for PMT array.\n\n", numRows, numInEachRow, numPMTs);
        exit(0);
    }

    // some histograms
    if (numPMTs > 0 && numRows > 0 && numInEachRow > 0)
    {
        hPMTArray = new TH2D("hPMTArray", "PMT hit pattern", numInEachRow, -0.5, numInEachRow - 0.5, numRows, -0.5, numRows - 0.5);
        hPMTNumVsTime = new TH2D("hPMTNumVsTime", "PMT number vs. time", 100, 0, 100, numPMTs, -0.5, numPMTs - 0.5);
    }
    else
    {
        hPMTArray = 0;
        hPMTNumVsTime = 0;
    }
    hPhotonEnergy = new TH1D("hPhotonEnergy", "Photon Energy", 500, 1.0, 10.0);
    hTotalPhotons = new TH1D("hTotalPhotons", "Total photons", 1000, 0, 10000);

    for (Int_t i = 0; i < 2; i++)
        hSidePMT[i] = 0;

    if (input->GetSideLining() == 1 && theDet->SidePMTX() > 0 && theDet->SidePMTZ() > 0)
    {
        for (Int_t i = 0; i < 2; i++)
        {
            hSidePMT[i] = new TH2D(Form("hSidePMT[%d]", i), Form("Side PMT %d", i),
                                   500, -theDet->SidePMTX() / 2.0, theDet->SidePMTX() / 2.0,
                                   500, -theDet->SidePMTZ() / 2.0, theDet->SidePMTZ() / 2.0);
        }
    }

    thePGA = (G4d2oPrimaryGeneratorAction *)G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();

    // ============================================================
    // Add instrumentation branches to Sim_Tree
    // ============================================================
    outTree->Branch("nReflections", &theEventData->nReflections, "nReflections/I");
    outTree->Branch("totalPathLength", &theEventData->totalPathLength, "totalPathLength/D");
    outTree->Branch("nPhotons", &theEventData->nPhotons, "nPhotons/I");
    outTree->Branch("meanReflections", &theEventData->meanReflections, "meanReflections/D");
    outTree->Branch("meanPathLength", &theEventData->meanPathLength, "meanPathLength/D");

    G4cerr << "done." << G4endl;
}

G4d2oEventAction::~G4d2oEventAction()
{
    G4cout << "Deleting G4d2oEventAction";

    if (outTree)
    {
        outTree->AutoSave();

        fout->cd();
        fout->mkdir("histos");
        fout->cd("histos");

        if (hPMTArray)
            hPMTArray->Write("pmtHits");
        if (hPMTNumVsTime)
            hPMTNumVsTime->Write("pmtNumVsTime");
        if (hPhotonEnergy)
            hPhotonEnergy->Write("photonEnergy");
        if (hTotalPhotons)
            hTotalPhotons->Write("photonsPerEvent");

        if (input->GetSideLining() == 1)
        {
            for (Int_t i = 0; i < 2; i++)
                if (hSidePMT[i])
                    hSidePMT[i]->Write(Form("sidePMTHits%d", i));
        }

        fout->mkdir("pmtPositions");
        fout->cd("pmtPositions");
        if (theDet->GetPMTPositionHistogram(0))
            theDet->GetPMTPositionHistogram(0)->Write("pmtPosX");
        if (theDet->GetPMTPositionHistogram(1))
            theDet->GetPMTPositionHistogram(1)->Write("pmtPosY");
        if (theDet->GetPMTPositionHistogram(2))
            theDet->GetPMTPositionHistogram(2)->Write("pmtPosZ");
    }

    G4String theOutFileName = fout->GetName();
    system("mkdir -p scripts");
    ofstream lastFile("scripts/lastFileName.txt");
    lastFile << fout->GetName() << endl;
    lastFile.close();

    delete fout;

    delete gBialkali_ev;

    G4cout << "...done" << G4endl << G4endl;

    G4cout << "\t**********************************************************************************" << G4endl;
    G4cout << "\t***   Output File: " << theOutFileName << " " << G4endl;
    G4cout << "\t**********************************************************************************" << G4endl;
    G4cout << G4endl << G4endl;
}

void G4d2oEventAction::BeginOfEventAction(const G4Event *thisEvent)
{
    fVetoEdep = 0.0;
    ResetNarrowVetoEdep();
    ResetWideVetoEdep();
}

void G4d2oEventAction::EndOfEventAction(const G4Event *thisEvent)
{
    HCoE = thisEvent->GetHCofThisEvent();
    if (HCoE)
    {
        GetHitsCollection();
        ProcessEvent();
    }
    for (int i = 0; i < kNarrowPanels; ++i)
        theEventData->narrow_veto_edep[i] = fNarrowVetoEdep[i];
    for (int i = 0; i < kWidePanels; ++i)
        theEventData->wide_veto_edep[i] = fWideVetoEdep[i];
}

void G4d2oEventAction::GetHitsCollection(void)
{
    for (Int_t i = 0; i < input->numDetHC; i++)
        detHC[i] = (G4d2oDetectorHitsCollection *)HCoE->GetHC(input->detCollID[i]);
    for (Int_t i = 0; i < input->numSrcHC; i++)
        sourceHC[i] = (G4d2oDetectorHitsCollection *)HCoE->GetHC(input->srcCollID[i]);
    for (Int_t i = 0; i < input->numTargHC; i++)
        targHC[i] = (G4d2oDetectorHitsCollection *)HCoE->GetHC(input->targCollID[i]);
}

void G4d2oEventAction::ProcessEvent(void)
{
    theEventData->ClearData();
    ZeroEventVariables();
    numEvents++;

    eventNumber = numEvents;

   
    G4d2oMaterialsDefinition::SetCurrentEventNumber(eventNumber);

    G4d2oSteppingAction* stepAct = (G4d2oSteppingAction*)G4RunManager::GetRunManager()->GetUserSteppingAction();
    if (stepAct) stepAct->SetEventID(eventNumber);

    G4d2oDetectorHit *thisHit;

    G4ThreeVector initDir = thePGA->GetInitialDirection();
    theEventData->direction0.SetXYZ(initDir.x(), initDir.y(), initDir.z());
    G4ThreeVector initPos = thePGA->GetOriginalPosition();
    theEventData->position0.SetXYZ(initPos.x(), initPos.y(), initPos.z());
    theEventData->sourceParticleEnergy = thePGA->GetSourceEnergy();
    theEventData->eventNumber = eventNumber;
    theEventData->vol0 = -1;

    G4Navigator *Navigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();

    G4VPhysicalVolume *volume0 = Navigator->LocateGlobalPointAndSetup(initPos);
    if (volume0->GetName() == "h2oPhysV")
    {
        theEventData->vol0 = 2;
    }
    else if (volume0->GetName() == "d2oPhysV")
    {
        theEventData->vol0 = 1;
    }
    else if (volume0->GetName() == "acrylicPhysV")
    {
        theEventData->vol0 = 3;
    }
    else if (volume0->GetName() == "pmtPhysV_0")
    {
        theEventData->vol0 = 4;
    }
    else if (volume0->GetName() == "pmtPhysV_1")
    {
        theEventData->vol0 = 5;
    }
    else if (volume0->GetName() == "pmtPhysV_2")
    {
        theEventData->vol0 = 6;
    }
    else if (volume0->GetName() == "pmtPhysV_3")
    {
        theEventData->vol0 = 7;
    }
    else if (volume0->GetName() == "pmtPhysV_4")
    {
        theEventData->vol0 = 8;
    }
    else if (volume0->GetName() == "pmtPhysV_5")
    {
        theEventData->vol0 = 9;
    }
    else if (volume0->GetName() == "pmtPhysV_6")
    {
        theEventData->vol0 = 10;
    }
    else if (volume0->GetName() == "pmtPhysV_7")
    {
        theEventData->vol0 = 11;
    }
    else if (volume0->GetName() == "pmtPhysV_8")
    {
        theEventData->vol0 = 12;
    }
    else if (volume0->GetName() == "pmtPhysV_9")
    {
        theEventData->vol0 = 13;
    }
    else if (volume0->GetName() == "pmtPhysV_10")
    {
        theEventData->vol0 = 14;
    }
    else if (volume0->GetName() == "pmtPhysV_11")
    {
        theEventData->vol0 = 15;
    }
    else if (volume0->GetName() == "pmtPhysV_12")
    {
        theEventData->vol0 = 16;
    }
    else if (volume0->GetName() == "pmtPhysV_13")
    {
        theEventData->vol0 = 17;
    }
    else if (volume0->GetName() == "pmtPhysV_14")
    {
        theEventData->vol0 = 18;
    }
    else if (volume0->GetName() == "pmtPhysV_15")
    {
        theEventData->vol0 = 19;
    }
    else if (volume0->GetName() == "pmtPhysV_16")
    {
        theEventData->vol0 = 20;
    }
    else if (volume0->GetName() == "pmtPhysV_17")
    {
        theEventData->vol0 = 21;
    }
    else if (volume0->GetName() == "pmtPhysV_18")
    {
        theEventData->vol0 = 22;
    }
    else if (volume0->GetName() == "pmtPhysV_19")
    {
        theEventData->vol0 = 23;
    }
    else if (volume0->GetName() == "pmtPhysV_20")
    {
        theEventData->vol0 = 24;
    }
    else if (volume0->GetName() == "pmtPhysV_21")
    {
        theEventData->vol0 = 25;
    }
    else if (volume0->GetName() == "pmtPhysV_22")
    {
        theEventData->vol0 = 26;
    }
    else if (volume0->GetName() == "pmtPhysV_23")
    {
        theEventData->vol0 = 27;
    }
    else if (volume0->GetName() == "pmtPhysV_24")
    {
        theEventData->vol0 = 28;
    }
    else if (volume0->GetName() == "pmtPhysV_25")
    {
        theEventData->vol0 = 29;
    }
    else if (volume0->GetName() == "pmtPhysV_26")
    {
        theEventData->vol0 = 30;
    }
    else if (volume0->GetName() == "pmtPhysV_27")
    {
        theEventData->vol0 = 31;
    }
    else if (volume0->GetName() == "pmtPhysV_28")
    {
        theEventData->vol0 = 32;
    }
    else if (volume0->GetName() == "pmtPhysV_29")
    {
        theEventData->vol0 = 33;
    }
    else if (volume0->GetName() == "pmtPhysV_30")
    {
        theEventData->vol0 = 34;
    }
    else if (volume0->GetName() == "pmtPhysV_31")
    {
        theEventData->vol0 = 35;
    }
    else if (volume0->GetName() == "pmtPhysV_32")
    {
        theEventData->vol0 = 36;
    }
    else if (volume0->GetName() == "pmtPhysV_33")
    {
        theEventData->vol0 = 37;
    }
    else if (volume0->GetName() == "pmtPhysV_34")
    {
        theEventData->vol0 = 38;
    }
    else if (volume0->GetName() == "pmtPhysV_35")
    {
        theEventData->vol0 = 39;
    }
    else if (volume0->GetName() == "pmtPhysV_36")
    {
        theEventData->vol0 = 40;
    }
    else if (volume0->GetName() == "pmtPhysV_37")
    {
        theEventData->vol0 = 41;
    }
    else if (volume0->GetName() == "pmtPhysV_38")
    {
        theEventData->vol0 = 42;
    }
    else if (volume0->GetName() == "pmtPhysV_39")
    {
        theEventData->vol0 = 43;
    }
    else if (volume0->GetName() == "pmtPhysV_40")
    {
        theEventData->vol0 = 44;
    }
    else if (volume0->GetName() == "pmtPhysV_41")
    {
        theEventData->vol0 = 45;
    }
    else if (volume0->GetName() == "pmtPhysV_42")
    {
        theEventData->vol0 = 46;
    }
    else if (volume0->GetName() == "pmtPhysV_43")
    {
        theEventData->vol0 = 47;
    }
    else if (volume0->GetName() == "pmtPhysV_44")
    {
        theEventData->vol0 = 48;
    }
    else if (volume0->GetName() == "pmtPhysV_45")
    {
        theEventData->vol0 = 49;
    }
    else if (volume0->GetName() == "pmtPhysV_46")
    {
        theEventData->vol0 = 50;
    }
    else if (volume0->GetName() == "pmtPhysV_47")
    {
        theEventData->vol0 = 51;
    }

    std::map<std::string, int> vmuid = {{"muVetoBI", 0}, {"muVetoBO", 1}, {"muVetoTI", 2}, {"muVetoTO", 3}, {"muVetoLI", 4}, {"muVetoLO", 5}, {"muVetoRI", 6}, {"muVetoRO", 7}, {"muVetoNI", 8}, {"muVetoNO", 9},
                                        {"muVetoFI", 10},
                                        {"muVetoFO", 11}};

    // Loop through the detector hits collection
    for (G4int iDetHC = 0; iDetHC < input->numDetHC; iDetHC++)
    {
        if (detHC[iDetHC])
        {
            TString SDname = detHC[iDetHC]->GetSDname().data();
            if (SDname.BeginsWith("muVetoInner"))
            {
                int muid = 0;
                for (G4int i = 0; i < detHC[iDetHC]->entries(); i++)
                {
                    thisHit = (*detHC[iDetHC])[i];
                    if (thisHit->GetPDGCharge() != 0)
                    {
                        theEventData->muVetoEnergy[muid] += thisHit->GetTotalEnergyDeposit();
                    }
                }
            }
            else if (SDname.BeginsWith("muVetoOuter"))
            {
                int muid = 1;
                for (G4int i = 0; i < detHC[iDetHC]->entries(); i++)
                {
                    thisHit = (*detHC[iDetHC])[i];
                    if (thisHit->GetPDGCharge() != 0)
                    {
                        theEventData->muVetoEnergy[muid] += thisHit->GetTotalEnergyDeposit();
                    }
                }
            }
            else if (SDname.BeginsWith("muVeto"))
            {
                int muid = vmuid[SDname.Data()];
                for (G4int i = 0; i < detHC[iDetHC]->entries(); i++)
                {
                    thisHit = (*detHC[iDetHC])[i];
                    if (thisHit->GetPDGCharge() != 0)
                    {
                        theEventData->muVetoEnergy[muid] += thisHit->GetTotalEnergyDeposit();
                    }
                }
            }
            else
            {
                for (G4int i = 0; i < detHC[iDetHC]->entries(); i++)
                {
                    thisHit = (*detHC[iDetHC])[i];

                    G4String hitParticleName = thisHit->GetParticleName();
                    if (hitParticleName == "opticalphoton")
                    {
                        G4int hitPMTNumber = thisHit->GetDetectorNumber();
                        G4double hitGlobalTime = thisHit->GetGlobalTime();
                        G4double hitKineticEnergy = thisHit->GetKineticEnergy();

                        // 02/14/2023: PMT QE scaling
#if 1
                        if (hitPMTNumber == 0) { if (G4UniformRand() > 0.8329) continue; }
                        if (hitPMTNumber == 1) { if (G4UniformRand() > 0.8197) continue; }
                        if (hitPMTNumber == 2) { if (G4UniformRand() > 0.9205) continue; }
                        if (hitPMTNumber == 3) { if (G4UniformRand() > 0.6747) continue; }
                        if (hitPMTNumber == 4) { if (G4UniformRand() > 0.8538) continue; }
                        if (hitPMTNumber == 5) { if (G4UniformRand() > 1.0000) continue; }
                        if (hitPMTNumber == 6) { if (G4UniformRand() > 0.9951) continue; }
                        if (hitPMTNumber == 7) { if (G4UniformRand() > 0.8934) continue; }
                        if (hitPMTNumber == 8) { if (G4UniformRand() > 0.7615) continue; }
                        if (hitPMTNumber == 9) { if (G4UniformRand() > 1.0000) continue; }
                        if (hitPMTNumber ==10) { if (G4UniformRand() > 0.6164) continue; }
                        if (hitPMTNumber ==11) { if (G4UniformRand() > 1.0000) continue; }
#endif

                        if (detHC[iDetHC]->GetSDname() == "pmt")
                        {
                            theEventData->AddPMTHit(hitPMTNumber, hitGlobalTime, hitKineticEnergy);
                            if (hPMTArray)
                                hPMTArray->Fill(hitPMTNumber / numRows, hitPMTNumber % numRows);
                            if (hPMTNumVsTime)
                                hPMTNumVsTime->Fill(hitGlobalTime, hitPMTNumber);
                        }
                        hPhotonEnergy->Fill(hitKineticEnergy / eV);
                    }
                }
            }

            if (detHC[iDetHC]->GetSDname() == "topVetoSD")
            {
                if (detHC[iDetHC]->GetSize() > 0)
                {
                    theEventData->veto_tag = true;
                    theEventData->veto_edep = fVetoEdep;
                }
                else
                {
                    theEventData->veto_tag = false;
                    theEventData->veto_edep = fVetoEdep;
                }
            }

            if (detHC[iDetHC]->GetSDname() == "VetoNu1SD") {
                for (int i=0;i<kNarrowPanels;++i)
                    theEventData->narrow_veto_edep[i] = fNarrowVetoEdep[i];
            }

            if (detHC[iDetHC]->GetSDname() == "VetoNu2SD") {
                for (int i=0;i<kWidePanels;++i)
                    theEventData->wide_veto_edep[i] = fWideVetoEdep[i];
            }
        }
    }


    if (stepAct) {
        theEventData->nReflections = stepAct->GetReflectionCount();
        theEventData->totalPathLength = stepAct->GetTotalPathLength() / mm;
        theEventData->nPhotons = stepAct->GetPhotonCount();
        if (theEventData->nPhotons > 0) {
            theEventData->meanReflections = (G4double)theEventData->nReflections / theEventData->nPhotons;
            theEventData->meanPathLength = theEventData->totalPathLength / theEventData->nPhotons;
        } else {
            theEventData->meanReflections = 0.0;
            theEventData->meanPathLength = 0.0;
        }
        stepAct->ResetEventCounters();
    }

    if (eventNumber % 1000 == 0)
    {
        G4cout << "Completed Event " << eventNumber << G4endl;
        if (input->GetPGAType() == 10)
        {
            G4cout << "Special Energy = " << input->GetSpecialEnergy() << G4endl;
        }
    }

    outTree->Fill();

    hTotalPhotons->Fill(theEventData->numHits);
}

void G4d2oEventAction::ZeroEventVariables(void)
{
    eventNumber = 0;
}

double G4d2oEventAction::qe_ev(double energy_ev)
{
    if (!gBialkali_ev)
        return 0.0;
    if (energy_ev < gBialkali_ev->GetX()[0] || energy_ev > gBialkali_ev->GetX()[gBialkali_ev->GetN() - 1])
        return 0.0;
    return gBialkali_ev->Eval(energy_ev, 0, "");
}
