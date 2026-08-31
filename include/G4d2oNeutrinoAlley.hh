#ifndef G4d2oNeutrinoAlley_H
#define G4d2oNeutrinoAlley_H 1

class G4LogicalVolume;
class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"
#include "G4d2oMaterialsDefinition.hh"

#include "G4d2oDetector.hh"

#include "inputVariables.hh"

#include "G4Tubs.hh"

class G4d2oNeutrinoAlley : public G4VUserDetectorConstruction
{
public:
	
	G4d2oNeutrinoAlley();
	~G4d2oNeutrinoAlley();
	
	virtual G4VPhysicalVolume* Construct();
    virtual G4d2oDetector * GetDetectorPtr() {return theDet;}
    virtual G4double GetTopOfConcreteCeiling() const {return topOfConcreteCeiling;}
protected:
	
private:
    
    const G4double PI;
    
    // Dimensions
    const G4double in;
    const G4double ft;
    
    //pointers
    G4d2oDetector *theDet;
    G4d2oMaterialsDefinition *matPtr;
    inputVariables *input;
    
    // Experimental Hall
    G4VSolid *worldVol_solid;
    G4LogicalVolume* worldVol_logV;
    G4VPhysicalVolume* worldVol_physV;
    
    G4VSolid *hallConcreteVol_solid;
    G4LogicalVolume* hallConcreteVol_logV;
    G4VPhysicalVolume* hallConcreteVol_physV;

    G4VSolid *hallOpening_solid;
    G4LogicalVolume* hallOpening_logV;
    G4VPhysicalVolume* hallOpening_physV;
    
    G4double worldVol_x;
    G4double worldVol_y;
    G4double worldVol_z;
    G4double topOfConcreteCeiling;
    G4double wallThickness;
    G4double hallway_height;
    G4double hallway_width;
    G4double hallway_length;
    G4double outerWall_y;
    G4double outerWall_z;
    G4double detectorcenter_offsetfromInnerWall;
    G4double walloffset_x;
    G4double walloffset_y;
    G4double walloffset_z;
    
};//END of class G4d2oNeutrinoAlley

#endif
