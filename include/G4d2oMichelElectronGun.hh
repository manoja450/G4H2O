#ifndef G4d2oMichelElectronGun_h
#define G4d2oMichelElectronGun_h 1

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

class G4d2oMichelElectronGun : public G4d2oPrimaryGeneratorAction
{
    RQ_OBJECT("G4d2oMichelElectronGun")

public:
    G4d2oMichelElectronGun();
    ~G4d2oMichelElectronGun();
    
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

    G4double EMMU; // muon mass
    G4double EMASS;
    G4double michel_rho;
    G4double michel_delta;
    G4double michel_xsi;
    G4double michel_eta;

    G4double F_c(G4double x, G4double x0);
    G4double F_theta(G4double x, G4double x0);
    G4double R_c(G4double x);
        
}; //END of class G4d2oMichelElectronGun

#endif

