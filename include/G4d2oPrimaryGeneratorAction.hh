#ifndef G4d2oPrimaryGeneratorAction_h
#define G4d2oPrimaryGeneratorAction_h 1

#include "globals.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"

#include "TTimer.h"

#include "G4d2oNeutrinoAlley.hh"
#include "G4d2oDetector.hh"

#include "ReplayTools.h"

class G4ParticleGun;
class G4Event;
class TFile;
class TH2D;

class G4d2oPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    G4d2oPrimaryGeneratorAction()
    : fluxFile(nullptr),
      fluxW(nullptr),
      angleOffset(0.0)
    {;}

    virtual ~G4d2oPrimaryGeneratorAction() {;}
    
public:
    virtual void GeneratePrimaries(G4Event*) {;}
	
    virtual G4double      GetSourceEnergy()     { return sourceEnergy; }
    virtual G4ThreeVector GetOriginalPosition() { return SourcePosition; }
    virtual G4ThreeVector GetInitialDirection() { return initDir; }
    
    ReplayTools *GetReplayTool() { return etfTimer; }
    
protected:
    
    G4ParticleGun* particleGun;
    
    G4double      sourceEnergy;
    G4ThreeVector SourcePosition;    
    G4ThreeVector initDir;

    G4d2oNeutrinoAlley *detCon;
    G4d2oDetector      *theDet;
    
    inputVariables *input;
    
    ReplayTools *etfTimer;

    G4ThreeVector tankSize;
    
    G4double zPlaneOffset;

    // --- NEW STUFF ---
    TFile*   fluxFile;     // will point to fluxWeight.root
    TH2D*    fluxW;        // the "fluxW" histogram
    G4double angleOffset;  // same offset as in analysis macro, turns out to be useless but whatever
    // -----------------------------------
 
}; // END of class G4d2oPrimaryGeneratorAction

#endif
