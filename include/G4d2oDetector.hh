#ifndef G4d2oDetector_H
#define G4d2oDetector_H 1

class G4LogicalVolume;
class G4VPhysicalVolume;

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4SDManager.hh"

#include "G4VUserDetectorConstruction.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPVParameterisation.hh"
#include "G4RotationMatrix.hh"

#include "TGraph.h"
#include "TH2F.h"

#include "G4d2oMaterialsDefinition.hh"



class G4d2oDetector
{
public:
    
    G4d2oDetector();
    ~G4d2oDetector();
    
    virtual G4LogicalVolume * GetDetector();
    G4ThreeVector GetPMTPosition(Int_t iPMT);
    TH1D * GetPMTPositionHistogram(Int_t iDim) {return hPanel[iDim];}
    G4int GetTotalPMTS() { return totPMT; }
    G4int GetPMTRows() { return numRows; }
    G4int GetPMTsInRow() { return numInRows; }
    
    G4double TankLength() {return d2oLength;} //x
    G4double TankWidth() {return d2oWidth;} //y
    G4double TankHeight() {return d2oHeight;} //z

    G4double OuterTankLength() {return h2oInnerLength;} //x
    G4double OuterTankWidth() {return h2oInnerWidth;} //y
    G4double OuterTankHeight() {return h2oInnerHeight;} //z
  
    G4double TankX() {return d2oLength;} //x
    G4double TankY() {return d2oWidth;} //y
    G4double TankZ() {return d2oHeight;} //z

    G4double OuterTankX() {return h2oInnerLength;} //x
    G4double OuterTankY() {return h2oInnerWidth;} //y
    G4double OuterTankZ() {return h2oInnerHeight;} //z

    G4double OuterReflectorX() {return h2oInnerLength;} //x
    G4double OuterReflectorY() {return h2oInnerWidth;} //y
    G4double OuterReflectorZ() {return h2oInnerHeight;} //z
    
    G4double OuterContX() {return contOuterLength;} //x
    G4double OuterContY() {return contOuterWidth;} //y
    G4double OuterContZ() {return contOuterHeight;} //z

    G4double OuterVetoX() {return vetoOuterLength;} //x
    G4double OuterVetoY() {return vetoOuterWidth;} //y
    G4double OuterVetoZ() {return vetoOuterHeight;} //z

    G4double SidePMTX() {return sideLiningX;}
    G4double SidePMTZ() {return sideLiningZ;}

    G4double EvalPMTQE(G4double energy_ev);
    
protected:

    // Dimensions
    const G4double in;
    
    inputVariables *input;
    G4d2oMaterialsDefinition *matPtr;
    G4SDManager* sdManager;

    G4LogicalVolume *GetTyvek();
    G4LogicalVolume *GetSphericalPMT();
    G4LogicalVolume *GetEllipsoidPMT();
    G4LogicalVolume *GetAreaPMT();
    virtual void Initialize();
    virtual void DetermineSpacing();

    TH1D *hPanel[3];
    G4int totPMT;
    G4int numRows;
    G4int numInRows;
    
    G4double bubbleRadius;
    G4double bubbleHeight;
    
    G4double d2oLength; //x
    G4double d2oWidth; //y
    G4double d2oHeight; //z

    G4double h2oTCthickness;

    G4double acrylicThickness;
    G4double pmtDepth;
    G4double tyvekThickness;
    G4double tyvekReflectivity;
    G4double tyvekSigmaAlpha; //roughness parameter. 0=polished
    G4double acrylicReflectivity;
    G4double acrylicSigmaAlpha;
    G4double airh2oReflectivity;
    G4double airh2oSigmaAlpha;
    G4double pmtDiameter; //controlled from beamOn.dat or command line now
    G4double pmtMinorAxis; // Minor axis for ellipsoidal pmts
    G4double pmtMajorAxis; // Major axis for ellipsoidal pmts
    G4double pmtLegLength;
    G4double muMetalThickness;
    G4double pmtWindowThickness;  
    G4double minGapBetweenPMTs;
    G4double outerContainerThickness;

    G4double shieldThickness;
    G4double h2oRefl;
    G4double shieldLength;
    G4double sink;
    G4double aluminumthickness;
    G4double airthickness;
    G4double shieldWidth;
    G4double shieldHeight;

    G4double h2oLength; //x
    G4double h2oWidth; //y
    G4double h2oHeight; //z
    G4double pmtReflectorLength; //x
    G4double pmtReflectorWidth; //y
    G4double pmtReflectorHeight; //z
    
    G4double h2oInnerLength; //x
    G4double h2oInnerWidth; //y
    G4double h2oInnerHeight; //z
    
    G4double sideLiningX, sideLiningZ;
    G4double theSideLiningThickness;
    
    G4double contOuterLength; //x
    G4double contOuterWidth; //y
    G4double contOuterHeight; //z

    G4double muonVetoThickness;
    G4int muonVetoLayers;

    G4double vetoOuterLength; //x
    G4double vetoOuterWidth; //y
    G4double vetoOuterHeight; //z

    
    G4int numInX, numInY, numInZ;
    G4double spacingInX, spacingInY, spacingInZ;

    TGraph * fPMTQE;
    G4double fminEnergyQE;
    G4double fmaxEnergyQE;
    G4int InitializeQEArray();

    
};//END of class G4d2oDetector

#endif
