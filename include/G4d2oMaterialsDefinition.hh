#ifndef G4d2oMaterialsDefinition_H
#define G4d2oMaterialsDefinition_H 1

#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4VPhysicalVolume.hh"
#include "G4SurfaceProperty.hh"

#include "inputVariables.hh"

class TFile;

enum materialName
{
    AIR, VINYLTOLUENE, ALUMINUM, LEAD,  
    POLY, STEEL,
    PMMA, MUMETAL,
    COPPER,
    PLASTIC,
    VACUUM, H2O, D2O,
    FUSEDSILICA, BOROSILICATE,
    PHOTOCATHODE,
    TEFLON,
    TYVEK,
    DELRIN,
    EPOXY,
    CONCRETE
};

class G4d2oMaterialsDefinition
{
public:
    
    G4d2oMaterialsDefinition();
    ~G4d2oMaterialsDefinition();
    
    G4Material * GetMaterial( materialName matName );
    static void SetReflector(G4VPhysicalVolume *theExitingVolume, G4VPhysicalVolume *theEnteringVolume,
		      G4double theReflectivity, G4double theSigmaAlpha=0.0, G4SurfaceType=dielectric_dielectric);
    
    static void SetDataDrivenReflector(G4VPhysicalVolume*, G4VPhysicalVolume*, G4double, const G4String&);

    static void AttachReflectionTree(TFile* file);
    static void CloseReflectionTree();
    static void SetCurrentEventNumber(G4int eventNum);
    
protected:
    inputVariables *input;
    G4double h2oRefl;
    
private:
    
    G4NistManager *manager;
    
    G4Material *matAir, *matVinylToluene, *matAl, *matSteel;
    G4Material *matPb, *matPoly;
    G4Material *matPMMA;
    G4Material *matMuMetal;
    G4Material *matCopper;
    G4Material *matPlastic;
    G4Material *matVacuum;
    G4Material *matH2O, *matD2O;
    G4Material *matFusedSilica;
    G4Material *matBorosilicate;
    G4Material *matPhotoCathode;
    G4Material *matTeflon;
    G4Material *matTyvek;
    G4Material *matDelrin;
    G4Material *matEpoxy;
    G4Material *matConcrete;
    
    void SetUniformOpticalProperties(G4Material *theMat, G4double theIndex, G4double theAbsLength);
    void SetOpticalProperties(G4Material *theMat, G4double theIndex, G4double h2oReflF, G4String absLengthFile);
    
};

#endif
