#include <fstream>
#include <iostream>

#include "G4d2oSensitiveDetector.hh"
#include "inputVariables.hh"
// #include "TreeMaker.hh"

#include "G4UnitsTable.hh"
#include "G4VProcess.hh"
#include "G4VParticleChange.hh"
#include "G4HadronicProcess.hh"
#include "G4Nucleus.hh"
#include "G4TrackStatus.hh"

#include "TH1.h"
#include "G4EventManager.hh"
#include "G4UserEventAction.hh"
#include "G4d2oEventAction.hh"
#include "G4SystemOfUnits.hh"

G4d2oSensitiveDetector::G4d2oSensitiveDetector(G4String sdName, G4String sdMaterial, G4int iType)
    : G4VSensitiveDetector(sdName),
      sensitiveDetectorName(sdName),
      sensitiveMaterialName(sdMaterial),
      collectionID(-1)
{
    G4cerr << "\tInstance of G4d2oSensitiveDetector Constructed." << G4endl;

    G4cerr << "\t\tSensitive Detector Name: " << sensitiveDetectorName << G4endl;
    G4cerr << "\t\tSensitive Material Name: " << sensitiveMaterialName << G4endl;

    //    hitsCollectionName = sensitiveMaterialName + "HitsCollection";
    hitsCollectionName = sensitiveDetectorName;
    G4cerr << "\t\tHits Collection Name: " << hitsCollectionName << G4endl;

    // Fills the G4CollectionNameVector collectionName with the name of the hits collection.
    G4d2oSensitiveDetector::collectionName.insert(hitsCollectionName);

    numHits = 0;

    // Get Date and Time
    dateTime = new TDatime();
    dateTime->Set();

    inputVariables *input = inputVariables::GetIVPointer();
    iprint = input->GetPrintStatus();

    if (iprint == 2)
        outfile = input->GetOutputFile();

    timeCutOff = 400 * ns; // time cut-off
    //	timeCutOff = 400e100*ns;  //time cut-off
    numEvents = 0;

    sdType = iType;

    collectionID = input->GetNewCollectionID(iType);
    G4cerr << "\t\tHits Collection ID: " << collectionID << G4endl;

} // END of private constructor

G4d2oSensitiveDetector::~G4d2oSensitiveDetector()
{
    G4cerr << "Deleting G4d2oSensitiveDetector (" << collectionID << ") ...";

    delete dateTime;

    G4cerr << "done" << G4endl;

} // END of destructor

void G4d2oSensitiveDetector::Initialize(G4HCofThisEvent *HCoE)
{

    if (collectionID < 0)
    {
        printf("ERROR: Hits collection ID not set properly.\n\n");
        exit(0);
    }

    hitsCollection = new G4d2oDetectorHitsCollection(sensitiveDetectorName, hitsCollectionName);

    HCoE->AddHitsCollection(collectionID, hitsCollection);

    numHits = 0;
    numEvents++;

} // END of Initialize()

G4bool G4d2oSensitiveDetector::ProcessHits(G4Step *theStep, G4TouchableHistory *roHist)
{
    roHist = 0;

    numHits++;

    if (numHits == 1)
        sumTime = 0;
    if (numHits == 1)
        sumED = 0;

    // Information that can be obtained from theStep
    G4Track *theTrack = theStep->GetTrack();
    const G4StepPoint *thePreStepPoint = theStep->GetPreStepPoint();
    const G4StepPoint *thePostStepPoint = theStep->GetPostStepPoint();
    theTotalEnergyDeposit = theStep->GetTotalEnergyDeposit();
    theDeltaTime = theStep->GetDeltaTime();

    const G4Material *mat0 = thePreStepPoint->GetMaterial();
    const G4Material *mat = thePostStepPoint->GetMaterial();
    G4String matname = mat0->GetName();
    if (mat)
        matname = mat->GetName();
    sprintf(theMaterialName, "%s", matname.c_str());

    // Information that can be obtained from theTrack
    const G4ParticleDefinition *theParticle = theTrack->GetDefinition();
    theTrackID = theTrack->GetTrackID();
    theTrackParentID = theTrack->GetParentID();
    theMomentumDirection = thePreStepPoint->GetMomentumDirection();
    thePostMomentumDirection = thePostStepPoint->GetMomentumDirection();

    // Information that can be obtained from theParticle
    G4String pName, pName2;
    pName = theParticle->GetParticleName();
    sprintf(theParticleName, "%s", pName.c_str());
    thePDGCharge = theParticle->GetPDGCharge();
    theBaryonNumber = theParticle->GetBaryonNumber();
    theLeptonNumber = theParticle->GetLeptonNumber();

    //  theGlobalTime = thePreStepPoint->GetGlobalTime();
    theGlobalTime = thePostStepPoint->GetGlobalTime();
    theKineticEnergy = thePreStepPoint->GetKineticEnergy();

    if (theGlobalTime > timeCutOff)
        theTrack->SetTrackStatus(fStopAndKill);

    //    if(collectionID == treeMakerControl->GetMaskCollID(0)){
    //        theTrack->SetTrackStatus(fStopAndKill);
    //        return true;
    //    }

    thePrePosition = thePreStepPoint->GetPosition();
    thePosition = thePostStepPoint->GetPosition();

    // Information that can be obtained from thePostStepPoint
    const G4VProcess *theProcessDefinedStep = (G4VProcess *)thePostStepPoint->GetProcessDefinedStep();
    if (theTrackID != 1)
    {
        const G4VProcess *theCreatorProcess = (G4VProcess *)theTrack->GetCreatorProcess();
        if (theCreatorProcess)
        {
            pName = theCreatorProcess->GetProcessName();
            sprintf(theCreatorProcessName, "%s", pName.c_str());
            pName = theCreatorProcess->GetProcessTypeName(theCreatorProcess->GetProcessType());
            sprintf(theCreatorProcessType, "%s", pName.c_str());
        }
        else
        {
            sprintf(theCreatorProcessName, "Unknown");
            sprintf(theCreatorProcessType, " ");
        }
    }
    else
    {
        sprintf(theCreatorProcessName, "Source");
        sprintf(theCreatorProcessType, " ");
    }

    theDetectorNumber = -1;
    if (theTrack->GetTouchable()->GetHistoryDepth() > 0)
        theDetectorNumber = theTrack->GetTouchable()->GetCopyNumber(2);
    thePixelNumber = 0;

    // Information that can be obtained from theProcessDefinedStep
    sprintf(theProcessName, "empty");
    sprintf(theProcessTypeName, "empty");
    if (theProcessDefinedStep)
    {
        pName = theProcessDefinedStep->GetProcessName();
        sprintf(theProcessName, "%s", pName.c_str());
        pName = theProcessDefinedStep->GetProcessTypeName(theProcessDefinedStep->GetProcessType());
        sprintf(theProcessTypeName, "%s", pName.c_str());
    }

    scatA = 0;
    scatZ = 0;

    if (theProcessDefinedStep->GetProcessType() == fHadronic && theProcessDefinedStep->GetProcessSubType() == fHadronElastic)
    { // see G4HadronicProcessType.hh

        if (theStep->GetSecondary() && theStep->GetSecondary()->size())
        {
            std::vector<G4Track *>::const_iterator it;

            // keeps a list of secondaries, so keep stepping through and get the last one
            for (it = theStep->GetSecondary()->begin(); it != theStep->GetSecondary()->end(); it++)
            {
                if ((*it)->GetDynamicParticle()->GetParticleDefinition()->GetAtomicNumber() == 0 ||
                    (*it)->GetPosition() != thePosition)
                    continue;

                scatA = (*it)->GetDynamicParticle()->GetParticleDefinition()->GetBaryonNumber();
                scatZ = (*it)->GetDynamicParticle()->GetParticleDefinition()->GetPDGCharge();
            }
        }
    }

    calcA = 0;
    if (theProcessDefinedStep->GetProcessType() == fHadronic && theProcessDefinedStep->GetProcessSubType() == fHadronElastic)
    { // see G4HadronicProcessType.hh

        G4double mN = 931.494061; // 1 amu
        G4double alpha = thePreStepPoint->GetKineticEnergy() - thePostStepPoint->GetKineticEnergy();
        G4ThreeVector p1 = thePreStepPoint->GetMomentum();
        G4ThreeVector p1p = thePostStepPoint->GetMomentum();
        G4double beta = (p1 - p1p).mag();

        calcA = (beta * beta - alpha * alpha) / (2.0 * alpha * mN);
    }

    if (iprint == 1)
    {
        G4cerr.setf(ios::fixed, ios::floatfield);
        if (numHits == 1 && sdType == 1)
        {
            G4cerr << G4endl;
            G4cerr << "Event # " << numEvents << G4endl;
        }
        G4cerr << " Hit # " << std::setw(3) << numHits << "  " << std::setw(3) << theTrackID << " " << std::setw(3) << theTrackParentID << " " << std::setw(3) << theDetectorNumber << " " << std::setw(3) << thePixelNumber;
        G4cerr << "  " << std::setw(16) << theCreatorProcessName;
        G4cerr << "  " << std::setw(10) << theParticleName;
        G4cerr.precision(4);
        G4cerr << "   KE = " << std::setw(8) << G4BestUnit(theKineticEnergy, "Energy");
        G4cerr.precision(2);
        G4cerr << "   gt = " << std::setw(6) << G4BestUnit(theGlobalTime, "Time");
        sumTime += theDeltaTime;
        G4cerr << "  " << std::setw(16) << theProcessName;
        G4cerr.precision(4);
        G4cerr << "  Ed = " << std::setw(8) << G4BestUnit(theTotalEnergyDeposit, "Energy");
        sumED += theTotalEnergyDeposit;
        G4cerr << "  Material = " << std::setw(16) << theMaterialName;
        G4cerr.precision(3);
        G4cerr << "   (px,py,pz)= " << std::setw(7) << theMomentumDirection.x() << " " << std::setw(7) << theMomentumDirection.y() << " " << std::setw(7) << theMomentumDirection.z();
        G4cerr.precision(1);
        G4cerr << "   (x,y,z)= " << std::setw(8) << thePosition.x() << " " << std::setw(8) << thePosition.y() << " " << std::setw(8) << thePosition.z();
        G4cerr << "   (x,y,z)= " << std::setw(8) << thePrePosition.x() << " " << std::setw(8) << thePrePosition.y() << " " << std::setw(8) << thePrePosition.z();
        G4cerr << G4endl;
    }
    if (iprint == 2)
    {
        outfile->setf(ios::fixed, ios::floatfield);
        if (numHits == 1 && sdType == 1)
        {
            *outfile << G4endl;
            *outfile << "Event # " << numEvents << G4endl;
        }
        *outfile << " Hit # " << std::setw(3) << numHits << "  " << std::setw(3) << theTrackID << " " << std::setw(3) << theTrackParentID << " " << std::setw(3) << theDetectorNumber << " " << std::setw(3) << thePixelNumber;
        *outfile << "  " << std::setw(16) << theCreatorProcessName;
        *outfile << "  " << std::setw(10) << theParticleName;
        outfile->precision(4);
        *outfile << "   KE = " << std::setw(8) << G4BestUnit(theKineticEnergy, "Energy");
        outfile->precision(2);
        *outfile << "   gt = " << std::setw(6) << G4BestUnit(theGlobalTime, "Time");
        sumTime += theDeltaTime;
        *outfile << "  " << std::setw(16) << theProcessName;
        outfile->precision(4);
        *outfile << "  Ed = " << std::setw(8) << G4BestUnit(theTotalEnergyDeposit, "Energy");
        sumED += theTotalEnergyDeposit;
        *outfile << "  Material = " << std::setw(16) << theMaterialName;
        outfile->precision(3);
        *outfile << "   (px,py,pz)= " << std::setw(7) << theMomentumDirection.x() << " " << std::setw(7) << theMomentumDirection.y() << " " << std::setw(7) << theMomentumDirection.z();
        outfile->precision(1);
        *outfile << "   (x,y,z)= " << std::setw(8) << thePosition.x() << " " << std::setw(8) << thePosition.y() << " " << std::setw(8) << thePosition.z();
        *outfile << "   (x,y,z)= " << std::setw(8) << thePrePosition.x() << " " << std::setw(8) << thePrePosition.y() << " " << std::setw(8) << thePrePosition.z();
        *outfile << G4endl;
    }

    // Insert a new hit into collection (keeps the boolean logic working)
    InsertNewHit();

    // kill/limits logic...
    if (matname == "Photo cathode (arb.)" && pName == "opticalphoton")
    {
        theTrack->SetTrackStatus(fStopAndKill);
    }
    else
    {
        if (numHits > 1e7)
        {
            G4cerr << " NUMHITS Exceeded " << std::endl;
            theTrack->SetTrackStatus(fStopAndKill);
        }
    }

    // ---- NEW: accumulate energy only for the top veto SD ----
    if (GetName() == "topVetoSD")
    {                                                           // guard to this SD instance
        const G4double edep = theStep->GetTotalEnergyDeposit(); // MeV
        if (edep > 0.)
        {
            auto *ua = G4EventManager::GetEventManager()->GetUserEventAction();
            auto *ea = static_cast<G4d2oEventAction *>(ua); // concrete EA

            if (ea)
            {
                ea->AddVetoEdep(edep);
            }
        }
    }

if (GetName() == "VetoNu1SD") {

    const G4double edep = theStep->GetTotalEnergyDeposit(); // MeV
    if (edep > 0.) {
      // panel index = parent's copy number (the aluminum mother)
      const auto touch = theStep->GetPreStepPoint()->GetTouchableHandle();
      const int panelIdx = touch->GetCopyNumber(1);  // <-- depth 1 = mother PV

      // Optional one-time debug to confirm depths:

      auto* ua = G4EventManager::GetEventManager()->GetUserEventAction();
      if (auto* ea = static_cast<G4d2oEventAction*>(ua)) {
        if (0 <= panelIdx && panelIdx < G4d2oEventAction::kNarrowPanels)
          ea->AddNarrowVetoEdep(panelIdx, edep);
      }
    }
  }

if (GetName() == "VetoNu2SD") {

    const G4double edep = theStep->GetTotalEnergyDeposit(); // MeV
    if (edep > 0.) {
      // panel index = parent's copy number (the aluminum mother)
      const auto touch = theStep->GetPreStepPoint()->GetTouchableHandle();
      const int panelIdx = touch->GetCopyNumber(1);  // <-- depth 1 = mother PV

      // Optional one-time debug to confirm depths:

      auto* ua = G4EventManager::GetEventManager()->GetUserEventAction();
      if (auto* ea = static_cast<G4d2oEventAction*>(ua)) {
        if (0 <= panelIdx && panelIdx < G4d2oEventAction::kWidePanels)
          ea->AddWideVetoEdep(panelIdx, edep);
      }
    }
  }


    // ---- end NEW ----
    return true;

} // END of ProcessHits()

void G4d2oSensitiveDetector::InsertNewHit()
{

    newHit = new G4d2oDetectorHit();
    hitsCollection->insert(newHit);

    newHit->SetPixelNumber(thePixelNumber);
    newHit->SetDetectorNumber(theDetectorNumber);
    newHit->SetTrackID(theTrackID);
    newHit->SetTrackParentID(theTrackParentID);
    newHit->SetScatA(scatA);
    newHit->SetScatZ(scatZ);

    newHit->SetTotalEnergyDeposit(theTotalEnergyDeposit);
    newHit->SetGlobalTime(theGlobalTime);
    newHit->SetKineticEnergy(theKineticEnergy);
    newHit->SetPDGCharge(thePDGCharge);
    newHit->SetLeptonNumber(theLeptonNumber);
    newHit->SetBaryonNumber(theBaryonNumber);
    newHit->SetDeltaTime(theDeltaTime);
    newHit->SetCalcA(calcA);

    newHit->SetParticleName(theParticleName);
    newHit->SetProcessName(theProcessName);
    newHit->SetProcessTypeName(theProcessTypeName);
    newHit->SetMaterialName(theMaterialName);
    newHit->SetCreatorProcessName(theCreatorProcessName);
    newHit->SetCreatorProcessType(theCreatorProcessType);

    newHit->SetMomentumDirection(theMomentumDirection);
    newHit->SetPostMomentumDirection(thePostMomentumDirection);
    newHit->SetPosition(thePosition);
    newHit->SetPrePosition(thePrePosition);
    newHit->SetVertex(theVertex);

} // END of InsertNewHit()

void G4d2oSensitiveDetector::EndOfEvent(G4HCofThisEvent *HCoE)
{
    HCoE = NULL;

} // END of EndOfEvent()

G4String G4d2oSensitiveDetector::GetSensitiveDetectorName(void)
{
    return sensitiveDetectorName;

} // END of GetSensitiveDetectorName()

G4String G4d2oSensitiveDetector::GetSensitiveMaterialName(void)
{
    return sensitiveMaterialName;

} // END of GetSensitiveMaterialName()

G4String G4d2oSensitiveDetector::GetHitsCollectionName(void)
{
    return hitsCollectionName;

} // END of GetHitsCollectionName()

