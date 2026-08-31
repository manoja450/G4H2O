#ifndef G4d2oSensitiveDetector_h
#define G4d2oSensitiveDetector_h 1

#include <fstream>

#include "G4d2oDetectorHit.hh"
//#include "TreeMaker.hh"

#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "globals.hh"

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TStopwatch.h"
#include "TH1.h"

using namespace std;

class G4d2oSensitiveDetector : public G4VSensitiveDetector
{
private:
    
    //static G4d2oSensitiveDetector *ptrSD;
    
    G4String sensitiveDetectorName;
    G4String sensitiveMaterialName;
    G4String hitsCollectionName;
    G4int sdType;
    
    // Hit and hits collection variables
    G4d2oDetectorHit *newHit;
    G4d2oDetectorHitsCollection *hitsCollection;
    G4int collectionID;
  	
    // Hit information variables
    G4ThreeVector theMomentumDirection;
    G4ThreeVector thePostMomentumDirection;
    G4ThreeVector thePosition, thePrePosition;
    G4ThreeVector theVertex;
    
    G4int thePixelNumber;
    G4int theDetectorNumber;
    G4int theTrackID;
    G4int theTrackParentID;
    G4int scatA, scatZ;

    G4double theTotalEnergyDeposit;
    G4double theGlobalTime;
    G4double theKineticEnergy;
    G4double thePDGCharge;
    G4double theLeptonNumber;
    G4double theBaryonNumber;
    G4double theDeltaTime;
    G4double calcA;
    
    char theParticleName[100];
    char theProcessName[100];
    char theProcessTypeName[100];
    char theMaterialName[100];
    char theCreatorProcessName[100];
    char theCreatorProcessType[100];
    
    // Counters
    G4int numEvents;
    G4int numHits;
    G4double sumTime;
    G4double sumED;
    G4double timeCutOff;
    
    //I/O variables
    ofstream * outfile;
    G4int iprint;
    
    //Date & Time variables
    TDatime *dateTime;
    
    //  TH1F *hSourceSpectrum[2];
//    TreeMaker *treeMakerControl;
    
protected:
    void InsertNewHit( void );
    
public:
	G4d2oSensitiveDetector(G4String sdName, G4String sdMaterial, G4int iType);
    ~G4d2oSensitiveDetector();
    //static G4d2oSensitiveDetector* GetSDPointer( G4String sdName, G4String sdMaterial );
    //static G4d2oSensitiveDetector* GetSDPointerIfExists( void );
    void Initialize(G4HCofThisEvent *HCoE);
    G4bool ProcessHits(G4Step *theStep,G4TouchableHistory *roHist); //roHist stands for "read-out history"
    void EndOfEvent(G4HCofThisEvent *HCoE);
    G4String GetSensitiveMaterialName( void );
    G4String GetSensitiveDetectorName( void );
    G4String GetHitsCollectionName( void );
    
}; //END of class G4d2oSensitiveDetector

#endif
