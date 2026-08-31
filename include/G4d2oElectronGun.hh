#ifndef G4d2oElectronGun_h
#define G4d2oElectronGun_h 1

#include "globals.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "G4DataVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4d2oPrimaryGeneratorAction.hh"

#include "inputVariables.hh"
#include "G4d2oNeutrinoAlley.hh"

#include "G4d2oDetector.hh"

#include "TF1.h"
#include <RQ_OBJECT.h>

class G4ParticleGun;
class G4Event;

class G4d2oElectronGun : public G4d2oPrimaryGeneratorAction
{
    RQ_OBJECT("G4d2oElectronGun")

public:
    G4d2oElectronGun();
    ~G4d2oElectronGun();
    
public:
    virtual void GeneratePrimaries(G4Event* anEvent);
	
    virtual G4double GetSourceEnergy() {return sourceEnergy;}
    virtual G4ThreeVector GetOriginalPosition() {return SourcePosition;}
    
protected:
    
private:
    
    G4int totalEvents;
    
    G4double cthetaRange1;
    G4double cthetaRange2;
    
    G4int theEventNum;
        
}; //END of class G4d2oElectronGun

#endif
