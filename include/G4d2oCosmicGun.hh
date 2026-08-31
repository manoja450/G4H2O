#ifndef G4d2oCosmicGun_h
#define G4d2oCosmicGun_h 1

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

#ifdef USE_CRY
#include "CRYSetup.h"
#include "CRYGenerator.h"
#include "CRYParticle.h"
#include "CRYUtils.h"
#endif

#include "RNGWrapper.hh"

#include <vector>

#include "TF1.h"
#include <RQ_OBJECT.h>



class G4ParticleGun;
class G4Event;

class G4d2oCosmicGun : public G4d2oPrimaryGeneratorAction
{
    RQ_OBJECT("G4d2oCosmicGun")

public:
    G4d2oCosmicGun();
    ~G4d2oCosmicGun();
    
public:
    virtual void GeneratePrimaries(G4Event* anEvent);
    void InputCRY();
    void UpdateCRY(std::string* MessInput);
    void CRYFromFile(G4String newValue);
    
    virtual G4double GetSourceEnergy() {return sourceEnergy;}
    virtual G4ThreeVector GetOriginalPosition() {return SourcePosition;}
    
protected:
    
private:
    
    G4int totalEvents;
    
    G4int theEventNum;
#ifdef USE_CRY
    std::vector<CRYParticle*> vect; // vector of generated particles
    CRYGenerator* gen;
#endif
    G4int InputState;
    
    G4double zPlaneOfCosmics;
    
}; //END of class G4d2oCosmicGun

#endif
