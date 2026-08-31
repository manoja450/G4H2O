#include "G4d2oNeutrinoAlley.hh"
#include "G4d2oDetector.hh"
#include "G4d2oCylindricalDetector.hh"
#include "G4d2oRunAction.hh"
#include "G4d2oSensitiveDetector.hh"
#include "inputVariables.hh"

#include "TMath.h"
#include "G4PhysicalVolumeStore.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4Tet.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "globals.hh"
#include "G4VisAttributes.hh"
#include "G4NistManager.hh"
#include "G4PVParameterised.hh"
#include "G4VisExtent.hh"
#include "G4SubtractionSolid.hh"
#include "G4TessellatedSolid.hh"
#include "G4QuadrangularFacet.hh"
#include "G4TriangularFacet.hh"
#include "G4VFacet.hh"
#include "G4GeometryManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"

G4d2oNeutrinoAlley::G4d2oNeutrinoAlley()
:  PI( CLHEP::pi), in(2.54*cm), ft(12.0*2.54*cm)
{ 
    input = inputVariables::GetIVPointer();
    
    worldVol_x = 11.0*m;    // Full length along X axis
    worldVol_y = worldVol_x; // Full length along Y axis
    worldVol_z = worldVol_x; // Full length along Z axis
    
    wallThickness = 8*ft;
    hallway_height = 12*ft; // was 3.7*m
    hallway_width = 11*ft; // was 3.3*m
    hallway_length = 8.0*m;
    outerWall_y = hallway_width + 2.0 * wallThickness;
    outerWall_z = hallway_height + 2.0 * wallThickness;
    detectorcenter_offsetfromInnerWall = 0.5*m;

    //create a flag for this eventually
//    theDet = new G4d2oDetector();
//    walloffset_x = 0.0;
//    walloffset_y = hallway_width/2.0-detectorcenter_offsetfromInnerWall;
//    walloffset_z = 0.0;
    theDet = new G4d2oCylindricalDetector();
    walloffset_x = 0.0;
    walloffset_y = hallway_width/2.0-(theDet->OuterVetoX()/2.0+1*cm);
    walloffset_z = 0.0;
    if(theDet){
        walloffset_z = hallway_height/2.0-(theDet->OuterVetoZ()/2.0+1*cm);
    }
    topOfConcreteCeiling=outerWall_z/2.0 + walloffset_z;
    std::cout<<"Setting Neutrino Alley TopOfCeiling to " <<topOfConcreteCeiling<<std::endl;

}//END of constructor

G4d2oNeutrinoAlley::~G4d2oNeutrinoAlley()
{
	
    G4cout<<"Deleting G4d2oNeutrinoAlley...";
    
    G4cout<<"done."<<G4endl;
    
}//END of destructor

G4VPhysicalVolume* G4d2oNeutrinoAlley::Construct()
{
	//Materials pointer
    matPtr = G4d2oRunAction::GetMaterialsPointer();
        
    ////////////////////////////////////////////////////////////////////
    //  World volume                                                  //
    ////////////////////////////////////////////////////////////////////
    
	worldVol_solid = new G4Box("worldVol_solid", worldVol_x/2.0, worldVol_y/2.0, worldVol_z/2.0);

	// Create a logical volume
	worldVol_logV = new G4LogicalVolume(worldVol_solid, matPtr->GetMaterial( AIR ), "worldVol_logV");
	worldVol_logV->SetVisAttributes(G4VisAttributes(false));
    
	// Create physical volume
	worldVol_physV = new G4PVPlacement(0, G4ThreeVector(0,0,0), worldVol_logV, "worldVol_physV", 0, false,	0);
    
    ////////////////////////////////////////////////////////////////////
    //  Concrete Hallway volume                                                  //
    ////////////////////////////////////////////////////////////////////
    // let's create the outer concrete such that detector stays centered in global coordinates even if offset in the hallway
    
    walloffset_x=0.0; walloffset_y = 0.0; walloffset_z = 0.0; //!! Commented this in to center detector in owrld volume!!
    
    hallConcreteVol_solid = new G4Box("hallConcreteVol_solid", hallway_length/2.0, outerWall_y/2.0, outerWall_z/2.0);
    
    // Create a logical volume
    hallConcreteVol_logV = new G4LogicalVolume(hallConcreteVol_solid, matPtr->GetMaterial( CONCRETE ), "hallConcreteVol_logV");
    //hallConcreteVol_logV->SetVisAttributes(G4VisAttributes(false));

    // Create physical volume
    hallConcreteVol_physV = new G4PVPlacement(0, G4ThreeVector(walloffset_x,walloffset_y,walloffset_z), hallConcreteVol_logV, "hallConcreteVol_physV", worldVol_logV, false, 0, true);
    G4VisAttributes *visAttHallO = new G4VisAttributes();
    visAttHallO->SetColour(G4Color::Blue());
    //    visAttHallO->SetForceSolid(true);
    hallConcreteVol_logV->SetVisAttributes( visAttHallO );

    G4VisAttributes * visWorld = new G4VisAttributes();
    visWorld->SetColor(G4Color::Red());
    //    visWorld->SetForceSolid(true);
    worldVol_logV->SetVisAttributes( visWorld );
    
    hallOpening_solid = new G4Box("hallOpening_solid", hallway_length/2.0, hallway_width/2.0, hallway_height/2.0);
    
    // Create a logical volume
    hallOpening_logV = new G4LogicalVolume(hallOpening_solid, matPtr->GetMaterial( AIR ), "hallOpening_logV");
    //hallConcreteVol_logV->SetVisAttributes(G4VisAttributes(false));
    G4VisAttributes *visAttHallI = new G4VisAttributes();
    visAttHallI->SetColour(G4Color::Blue());
    //    visAttHallI->SetForceSolid(true);
    hallOpening_logV->SetVisAttributes( visAttHallI );
    
    // Create physical volume
    hallOpening_physV = new G4PVPlacement(0, G4ThreeVector(0,0,0), hallOpening_logV, "hallOpening_physV", hallConcreteVol_logV, false, 0, true);

    
    if(theDet){
        G4LogicalVolume *detLogV = (G4LogicalVolume*)theDet->GetDetector();

        if(detLogV){
            //detector into worldVol
            //new G4PVPlacement(0,G4ThreeVector(-walloffset_x,-walloffset_y,-walloffset_z),detLogV,"detPhysV",hallOpening_logV,false,0,true);
            
            //detector into worldVol. ADDED 5*in IN Z-COORDINATE for WATER TC!!!! 
            G4RotationMatrix *rotDet = new G4RotationMatrix(); //Accomodating Detector with PMTs on the sides:
            rotDet->rotateZ(-90.0*deg);

	  
	    //	    vetoOuterHeight = detLogV->vetoOuterHeight;
	    G4VPhysicalVolume *totalDetPhysV = new G4PVPlacement(rotDet,G4ThreeVector(-walloffset_x,-44.808/2.0*in+hallway_width/2.0,-hallway_height/2.0+88.67*in/2.0+0.38*in/2.0),detLogV,"detPhysV",hallOpening_logV,false,0,true);
	    G4d2oMaterialsDefinition::SetReflector(totalDetPhysV,G4PhysicalVolumeStore::GetInstance()->GetVolume("h2oPhysV"),0.02,0.0);
	    //0.5inches ir airthickness of veto panels
	    //            new G4PVPlacement(rotDet,G4ThreeVector(-walloffset_x,-walloffset_y,-(hallway_height)/2+1.0*(191.854)/2),detLogV,"detPhysV",hallOpening_logV,false,0,true);
	    //191.854*cm was outerVetoHeight on October 03, 2024
	    // was	    new G4PVPlacement(rotDet,G4ThreeVector(-walloffset_x,-walloffset_y,-(walloffset_z/2.0) + 5.0*in),detLogV,"detPhysV",hallOpening_logV,false,0,true);
             
        }
    }
    
    return worldVol_physV;
    
}//END of Construct()

