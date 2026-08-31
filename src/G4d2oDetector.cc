#include "G4d2oDetector.hh"
#include "G4d2oRunAction.hh"
#include "G4d2oSensitiveDetector.hh"

#include "TMath.h"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Ellipsoid.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4VisAttributes.hh"
#include "G4PVParameterised.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

#include "G4OpticalSurface.hh" //NEW
#include "G4LogicalBorderSurface.hh" //NEW

#include "G4VisAttributes.hh" //NEW

using namespace std;

G4d2oDetector::G4d2oDetector()
: in(2.54*cm)
{
    G4cerr << "\tConstructing G4d2oDetector..." ;

    input = inputVariables::GetIVPointer();
    
    for(G4int i=0; i<3; i++) hPanel[i] = 0;
    totPMT = 0;
    numRows = 0;
    numInRows = 0;
    
    Initialize();
    
    sdManager = G4SDManager::GetSDMpointer();

}

void G4d2oDetector::Initialize(){
    
    //relevant parameters
    d2oLength = 100.0*cm;//140.0*cm; //x
    d2oWidth = 60.0*cm; //y
    d2oHeight = 100.0*cm;//140.0*cm; //z
    acrylicThickness = 1.0*in;
    pmtDepth = 10.0*in;
    tyvekThickness = 0.25*in;
    tyvekReflectivity = 0.99; //need to verify
    tyvekSigmaAlpha = 0.1; //roughness parameter. use 0 for polished
    muMetalThickness = 0.03*in;
    pmtWindowThickness = 3.0*mm;  //making this up
    minGapBetweenPMTs = 1.0*mm;
    //going for 10cm of H2O around entire D2O tank
    h2oTCthickness = input->GetTailCatcherThick()*cm;

    //pmtDiameter controlled from beamOn.dat or command line now
    // if pmtDiameter is negative we will use the ellipsoidal pmts!
    pmtDiameter = input->GetPMTDiameter()*in;
    pmtMinorAxis = pmtDiameter;
    pmtMajorAxis = pmtDiameter;
    pmtLegLength = 62.*mm; // was 62.*mm
    if (pmtDiameter < 0.0) {
      pmtMinorAxis = 76.72*mm * 2.0;
      pmtMajorAxis = 103.6*mm * 2.0;
    }

    h2oLength = d2oLength + 2*(acrylicThickness + h2oTCthickness + pmtDepth); //x
    h2oWidth = d2oWidth + 2*(acrylicThickness + h2oTCthickness); //y
    //h2oHeight = d2oHeight + 2*(acrylicThickness + h2oTCthickness + pmtDepth); //z
    if (pmtDiameter > 0.0) {
      h2oHeight = d2oHeight + 2*(acrylicThickness + h2oTCthickness + pmtMinorAxis/2.0); //z
    } else {
      h2oHeight = d2oHeight + 2*(acrylicThickness + h2oTCthickness + pmtMinorAxis + pmtLegLength); //z
    }
    outerContainerThickness = 0.25*in;

    shieldThickness = 0.0*in;
    shieldLength = shieldWidth = shieldHeight = 0.0*in;

    ///// Full detector volume /////
    contOuterLength = h2oLength + 2*outerContainerThickness; //x
    contOuterWidth = h2oWidth + 2*theSideLiningThickness + 2*outerContainerThickness; //y
    contOuterHeight = h2oHeight + 2*outerContainerThickness; //z
    
    //Reflector in PMT arrays
    pmtReflectorLength = d2oLength + 2*(acrylicThickness + h2oTCthickness + (pmtMinorAxis/2.0) + tyvekThickness);
    pmtReflectorWidth = d2oWidth + 2*(acrylicThickness + h2oTCthickness);
    pmtReflectorHeight = d2oHeight + 2*(acrylicThickness + h2oTCthickness + (pmtMinorAxis/2.0) + tyvekThickness);
    
    h2oInnerLength = d2oLength + 2*(acrylicThickness + h2oTCthickness + (pmtMinorAxis/2.0) ); //x
    h2oInnerWidth = d2oWidth + 2*(acrylicThickness + h2oTCthickness); //y
    h2oInnerHeight = d2oHeight + 2*(acrylicThickness + h2oTCthickness + (pmtMinorAxis/2.0)); //z
    
    muonVetoThickness = 2.54*cm;
    muonVetoLayers = 2;
    
    vetoOuterLength = contOuterLength + 2.0*muonVetoLayers*muonVetoThickness; //x
    vetoOuterWidth = contOuterWidth + 2.0*muonVetoLayers*muonVetoThickness; //y
    vetoOuterHeight = contOuterHeight + 2.0*muonVetoLayers*muonVetoThickness; //z

    //DetermineSpacing must be called after all the parameters above are defined
    DetermineSpacing();

    fPMTQE = 0;
    InitializeQEArray();
        
    G4cerr << "done." << G4endl;
    
}//END of constructor

G4d2oDetector::~G4d2oDetector()
{
    if (fPMTQE) delete fPMTQE;

    G4cout << "Instance of G4d2oDetector Destructed!" << G4endl;
    
}//END of destructor

void G4d2oDetector::DetermineSpacing(){
        
    //minimum gap between PMTs and between PMTs and edge of volume
    G4double nominalSpacing = (pmtMajorAxis) + minGapBetweenPMTs;
    
    //spacing along height (z)
    //-0.5*nominalSpacing b/c we stagger, so need to account for not putting them out over the edge
    G4double totalHeightInZ = h2oHeight-2*pmtDepth-0.5*nominalSpacing-2*minGapBetweenPMTs;
    numInZ = TMath::FloorNint(totalHeightInZ/nominalSpacing);
    spacingInZ = totalHeightInZ/double(numInZ);
 
    //spacing along length (x)
    G4double totalLengthInX = h2oLength-2*pmtDepth-0.5*nominalSpacing-2*minGapBetweenPMTs;
    numInX = TMath::FloorNint(totalLengthInX/nominalSpacing);
    spacingInX = totalLengthInX/double(numInX);

    G4double nominalRowSpacing = TMath::Sqrt( pow(pmtMajorAxis,2.0) - pow(0.5*spacingInZ,2.0) ) + minGapBetweenPMTs;
    
    //here we fix total width to use in calculation b/c we don't want to overlap the edge,
    //thus -2*(pmtDiameter-nominalRowSpacing)
    G4double totalWidthInY = h2oWidth  - 2*minGapBetweenPMTs - 2*((pmtMajorAxis/2.0)-nominalRowSpacing);
    numInY = TMath::FloorNint(totalWidthInY/nominalRowSpacing);
    spacingInY = totalWidthInY/double(numInY);
    
    //some totals that get passed to other classes
    numRows = numInY;
    numInRows = 2*(numInZ + numInX);
    totPMT = numRows*numInRows;
    
    //sideLiningX = h2oLength-2*pmtDepth;
    //sideLiningZ = h2oHeight-2*pmtDepth;

    sideLiningX = h2oInnerLength;
    sideLiningZ = h2oInnerHeight;


}

G4LogicalVolume * G4d2oDetector::GetDetector(){

    //Materials pointer
    matPtr = G4d2oRunAction::GetMaterialsPointer();
    
//     Geometry hierarchy
//     totalDetLogV (mother volume)
//        -> H2O tank
//            -> PMTs
//            -> acrylic tank
//                -> D2O
//        -> Outer tyvek tyvek layers
    
    //define the thickness of the side lining
    if(input->GetSideLining()==0) theSideLiningThickness = tyvekThickness;
    else if(input->GetSideLining()==1) theSideLiningThickness = pmtWindowThickness;
    else{
        printf("\tWarning in G4d2oDetector::GetDetector() - side lining specification unknown (%d)\n",input->GetSideLining());
    }
    
    G4Box* outerSolid = new G4Box("outerSolid",vetoOuterLength/2.0,vetoOuterWidth/2.0,vetoOuterHeight/2.0);
    G4LogicalVolume *totalDetLogV = new G4LogicalVolume(outerSolid,
                                                        matPtr->GetMaterial( AIR ),
                                                        "totalDetLogV");
    /// Muon Veto Botton
    G4Box* muVetoBISolid = new G4Box("muVetoBISolid",
                                     (vetoOuterLength-2.0*muonVetoThickness*muonVetoLayers)/2.0,
                                     (vetoOuterWidth-2.0*muonVetoThickness*muonVetoLayers)/2.0,
                                     muonVetoThickness/2.0);
    G4LogicalVolume *muVetoBILogV = new G4LogicalVolume(muVetoBISolid,
                                                        matPtr->GetMaterial( PLASTIC ),
                                                        "muVetoBILogV");
    G4Box* muVetoBOSolid = new G4Box("muVetoBOSolid",
                                     (vetoOuterLength-2.0*muonVetoThickness*muonVetoLayers)/2.0,
                                     (vetoOuterWidth-2.0*muonVetoThickness*muonVetoLayers)/2.0,
                                     muonVetoThickness/2.0);
    G4LogicalVolume *muVetoBOLogV = new G4LogicalVolume(muVetoBOSolid,
                                                        matPtr->GetMaterial( PLASTIC ),
                                                        "muVetoBOLogV");
    /// Muon Veto Top
    G4Box* muVetoTISolid = new G4Box("muVetoTISolid",vetoOuterLength/2.0,
                                                     vetoOuterWidth/2.0,
                                                     muonVetoThickness/2.0);
    G4LogicalVolume *muVetoTILogV = new G4LogicalVolume(muVetoTISolid,
                                                        matPtr->GetMaterial( PLASTIC ),
                                                        "muVetoTILogV");
    G4Box* muVetoTOSolid = new G4Box("muVetoTOSolid",vetoOuterLength/2.0,vetoOuterWidth/2.0,muonVetoThickness/2.0);
    G4LogicalVolume *muVetoTOLogV = new G4LogicalVolume(muVetoTOSolid,
                                                        matPtr->GetMaterial( PLASTIC ),
                                                        "muVetoTOLogV");
    /// Muon Veto Left Side
    G4Box* muVetoLISolid = new G4Box("muVetoLISolid",
                                     muonVetoThickness/2.0,
                                     (vetoOuterWidth-2.0*muonVetoThickness*muonVetoLayers)/2.0,
                                     (vetoOuterHeight-muonVetoThickness*muonVetoLayers)/2.0);
    G4LogicalVolume *muVetoLILogV = new G4LogicalVolume(muVetoLISolid,
                                                        matPtr->GetMaterial( PLASTIC ),
                                                        "muVetoLILogV");
    G4Box* muVetoLOSolid = new G4Box("muVetoLOSolid",
                                     muonVetoThickness/2.0,
                                     (vetoOuterWidth-2.0*muonVetoThickness*muonVetoLayers)/2.0,
                                     (vetoOuterHeight-muonVetoThickness*muonVetoLayers)/2.0);
    G4LogicalVolume *muVetoLOLogV = new G4LogicalVolume(muVetoLOSolid,
                                                        matPtr->GetMaterial( PLASTIC ),
                                                        "muVetoLOLogV");
    /// Muon Veto Right Side
    G4Box* muVetoRISolid = new G4Box("muVetoRISolid",
                                     muonVetoThickness/2.0,
                                     (vetoOuterWidth-2.0*muonVetoThickness*muonVetoLayers)/2.0,
                                     (vetoOuterHeight-muonVetoThickness*muonVetoLayers)/2.0);
    G4LogicalVolume *muVetoRILogV = new G4LogicalVolume(muVetoRISolid,
                                                        matPtr->GetMaterial( PLASTIC ),
                                                        "muVetoRILogV");
    G4Box* muVetoROSolid = new G4Box("muVetoROSolid",
                                     muonVetoThickness/2.0,
                                     (vetoOuterWidth-2.0*muonVetoThickness*muonVetoLayers)/2.0,
                                     (vetoOuterHeight-muonVetoThickness*muonVetoLayers)/2.0);
    G4LogicalVolume *muVetoROLogV = new G4LogicalVolume(muVetoROSolid,
                                                        matPtr->GetMaterial( PLASTIC ),
                                                        "muVetoROLogV");
    /// Muon Veto Far Side
    G4Box* muVetoFISolid = new G4Box("muVetoFISolid",
                                     vetoOuterLength/2.0,
                                     muonVetoThickness/2.0,
                                     (vetoOuterHeight-muonVetoThickness*muonVetoLayers)/2.0);
    G4LogicalVolume *muVetoFILogV = new G4LogicalVolume(muVetoFISolid,
                                                        matPtr->GetMaterial( PLASTIC ),
                                                        "muVetoFILogV");
    G4Box* muVetoFOSolid = new G4Box("muVetoFOSolid",
                                     vetoOuterLength/2.0,
                                     muonVetoThickness/2.0,
                                     (vetoOuterHeight-muonVetoThickness*muonVetoLayers)/2.0);
    G4LogicalVolume *muVetoFOLogV = new G4LogicalVolume(muVetoFOSolid,
                                                        matPtr->GetMaterial( PLASTIC ),
                                                        "muVetoFOLogV");
    /// Muon Veto Near Side
    G4Box* muVetoNISolid = new G4Box("muVetoNISolid",
                                     vetoOuterLength/2.0,
                                     muonVetoThickness/2.0,
                                     (vetoOuterHeight-muonVetoThickness*muonVetoLayers)/2.0);
    G4LogicalVolume *muVetoNILogV = new G4LogicalVolume(muVetoNISolid,
                                                        matPtr->GetMaterial( PLASTIC ),
                                                        "muVetoNILogV");
    G4Box* muVetoNOSolid = new G4Box("muVetoNOSolid",
                                     vetoOuterLength/2.0,
                                     muonVetoThickness/2.0,
                                     (vetoOuterHeight-muonVetoThickness*muonVetoLayers)/2.0);
    G4LogicalVolume *muVetoNOLogV = new G4LogicalVolume(muVetoNOSolid,
                                                        matPtr->GetMaterial( PLASTIC ),
                                                        "muVetoNOLogV");

    //Full detector volume will be steel for now. Steel will be displaced
    //as we add materials inside it. Steel has no optical parameters, so any photons that pass from
    //the H2O to the steel will be killed
    G4Box *outerVessel = new G4Box("outerVessel",contOuterLength/2.0,contOuterWidth/2.0,contOuterHeight/2.0);
    G4LogicalVolume *outerVesselLogV = new G4LogicalVolume(outerVessel,
                                                        matPtr->GetMaterial( STEEL ),
                                                        "outerVesselLogV");

   // outerVesselLogV->SetVisAttributes(G4VisAttributes::Invisible);
    outerVesselLogV->SetVisAttributes(G4VisAttributes::GetInvisible()); //new for use in topaz!!

 
    
    
    ///// D2O /////
    G4Box *d2oSolid = new G4Box("d2oSolid",d2oLength/2.0,d2oWidth/2.0,d2oHeight/2.0);
    G4LogicalVolume *d2oLogV = new G4LogicalVolume(d2oSolid,
                                                   matPtr->GetMaterial( D2O ),
                                                   "d2oLogV");
    G4VisAttributes *visAttD2O = new G4VisAttributes();
    visAttD2O->SetColour(G4Color::Blue());
    //    visAttD2O->SetForceSolid(true);
    d2oLogV->SetVisAttributes( visAttD2O );

    ///// PMMA /////
    G4double acrylicLength = d2oLength + 2.0*acrylicThickness; //x
    G4double acrylicWidth = d2oWidth + 2.0*acrylicThickness; //y
    G4double acrylicHeight = d2oHeight + 2.0*acrylicThickness; //z
    
    G4Box *acrylicSolid = new G4Box("acrylicSolid",acrylicLength/2.0,acrylicWidth/2.0,acrylicHeight/2.0);
    G4LogicalVolume *acrylicLogV = new G4LogicalVolume(acrylicSolid,
                                                   matPtr->GetMaterial( PMMA ),
                                                   "acrylicLogV");

    ///// H2O Inner /////
    G4Box *h2oSolid = new G4Box("h2oSolid",h2oLength/2.0,h2oWidth/2.0,h2oHeight/2.0);
    G4LogicalVolume *h2oLogV = new G4LogicalVolume(h2oSolid,
                                                   matPtr->GetMaterial( H2O ),
                                                   "h2oLogV");
    
    ///// TYVEK Outer /////
    G4Box *tyvekPMTReflectorO = new G4Box("tyvekPMTReflectorO",pmtReflectorLength/2.0,pmtReflectorWidth/2.0,pmtReflectorHeight/2.0);
    G4Box *tyvekPMTReflectorI = new G4Box("tyvekPMTReflectorI",h2oInnerLength/2.0,h2oInnerWidth/2.0,h2oInnerHeight/2.0);

    G4cout<<"pmtReflector("<<pmtReflectorLength
          <<","<<pmtReflectorWidth
          <<","<<pmtReflectorHeight<<")"<<G4endl;
    G4cout<<"h2oInner("<<h2oInnerLength
          <<","<<h2oInnerWidth
          <<","<<h2oInnerHeight<<")"<<G4endl;
    
    G4SubtractionSolid* tyvekPMTReflectorP = new G4SubtractionSolid("tyvekPMTReflectorP",tyvekPMTReflectorO,tyvekPMTReflectorI,0,G4ThreeVector(0,0,0));
    int pmtSubI = 1;
    
    G4VisAttributes *visLightWater = new G4VisAttributes();
    visLightWater->SetColour(G4Color::Cyan());
    //    visLightWater->SetForceSolid(true);
    acrylicLogV->SetVisAttributes( visLightWater );
    h2oLogV->SetVisAttributes( visLightWater );
    
    
    ///// PMTs /////
    G4LogicalVolume *pmtLogV = 0;
    if (pmtDiameter < 1.0) {
      pmtLogV = GetEllipsoidPMT();
    } else {
      pmtLogV = GetSphericalPMT();
    }
    
    //Mu-metal housing
    //at present, this is included to ensure that photons don't enter PMT volume from behind

    G4double pmtHousingDepth = pmtDepth-(pmtMinorAxis/2.0);
    G4Tubs *pmtHousingSolid = new G4Tubs("pmtHousingSolid",0.0,pmtMajorAxis/2.0,
                                         pmtHousingDepth/2.0,0.0,360.0*deg);
    G4LogicalVolume *pmtHousingLogV = new G4LogicalVolume(pmtHousingSolid,
                                                          matPtr->GetMaterial( VACUUM ),
                                                          "pmtHousingLogV");
    G4Tubs *pmtMuMetalSolid = new G4Tubs("pmtMuMetalSolid",(pmtMajorAxis/2.0)-muMetalThickness,pmtMinorAxis/2.0,
                                         pmtHousingDepth/2.0,0.0,360.0*deg);
    G4LogicalVolume *pmtMuMetalLogV = new G4LogicalVolume(pmtMuMetalSolid,
                                                          matPtr->GetMaterial( MUMETAL ),
                                                          "pmtMuMetalLogV");
    new G4PVPlacement(0,G4ThreeVector(0,0,0),pmtMuMetalLogV,"pmtMuMetalPhysV",pmtHousingLogV,false,0,true);


    ///// Place the PMTs /////
    
    //histograms that hold coordinates of the front faces of the PMTs
    for(G4int iDim=0; iDim<3; iDim++)
        hPanel[iDim] = new TH1D(Form("hPanel[%d]",iDim),Form("Dimension %d",iDim),totPMT,-0.5,totPMT-0.5);

    G4int iCopy = 0;
    //y-z plane, +x
    G4RotationMatrix *rotA = new G4RotationMatrix();
    rotA->rotateY(90.0*deg);
    for(G4int iCopyZ = 0; iCopyZ<numInZ; iCopyZ++){
        for(G4int iCopyY = 0; iCopyY<numInY; iCopyY++){
            
            G4double theY = -spacingInY*numInY/2.0 + double(iCopyY+0.5)*spacingInY;
            G4double theShift = 0.25;
            if(iCopyY%2==0) theShift = -0.25;
	    // messing with "theZ"
	    G4double theZ = -spacingInZ*numInZ/2.0 + double(iCopyZ+0.5)*spacingInZ + theShift*spacingInZ + 500.0*mm;
	    //            G4double theZ = -spacingInZ*numInZ/2.0 + double(iCopyZ+0.5)*spacingInZ + theShift*spacingInZ;
	    // done messing with "theZ"
            G4ThreeVector thisPMTPlacement(h2oLength/2.0-pmtHousingDepth,theY,theZ);
            new G4PVPlacement(rotA,thisPMTPlacement,pmtLogV,Form("pmtPhysV_%d",iCopy),h2oLogV,true,iCopy,true);

            G4ThreeVector thisPMTHousingPlacement((h2oLength-pmtHousingDepth)/2.0,theY,theZ);
            new G4PVPlacement(rotA,thisPMTHousingPlacement,pmtHousingLogV,Form("pmtHousingPhysV_%d",iCopy),h2oLogV,true,iCopy,true);
           
            //place holes in reflector for PMT
            tyvekPMTReflectorP = new G4SubtractionSolid(Form("tyvekPMTReflectorP%d",pmtSubI),tyvekPMTReflectorP,
                                                         new G4Tubs(Form("pmtMuMetalSolid_%d",pmtSubI),0.0,pmtMinorAxis/2.0,
                                                                                                                  pmtHousingDepth,0.0,360.0*deg),rotA,thisPMTHousingPlacement);
            
            pmtSubI++;
            iCopy++;

            //save coordinates of the front faces of the PMTs
            hPanel[0]->SetBinContent(iCopy,thisPMTPlacement.x()-(pmtMinorAxis/2.0));
            hPanel[1]->SetBinContent(iCopy,thisPMTPlacement.y());
            hPanel[2]->SetBinContent(iCopy,thisPMTPlacement.z());
        }
    }

    //x-y plane, +z
    G4RotationMatrix *rotB = new G4RotationMatrix();
    rotB->rotateY(180.0*deg);
    //if(false)
    for(G4int iCopyX = 0; iCopyX<numInX; iCopyX++){
        for(G4int iCopyY = 0; iCopyY<numInY; iCopyY++){
            
            G4double theY = -spacingInY*numInY/2.0 + double(iCopyY+0.5)*spacingInY;
            G4double theShift = 0.25;
            if(iCopyY%2==0) theShift = -0.25;
            
            G4double theX = spacingInX*numInX/2.0 - (double(iCopyX+0.5)*spacingInX + theShift*spacingInX);
            G4ThreeVector thisPMTPlacement(theX,theY,h2oHeight/2.0-pmtHousingDepth+500.0*mm); // -50*mm
            new G4PVPlacement(rotB,thisPMTPlacement,pmtLogV,Form("pmtPhysV_%d",iCopy),h2oLogV,true,iCopy,true);
            
            G4ThreeVector thisPMTHousingPlacement(theX,theY,(h2oHeight-pmtHousingDepth)/2.0);
            new G4PVPlacement(rotB,thisPMTHousingPlacement,pmtHousingLogV,Form("pmtHousingPhysV_%d",iCopy),h2oLogV,true,iCopy,true);
            
            //place holes in reflector for PMT
            tyvekPMTReflectorP = new G4SubtractionSolid(Form("tyvekPMTReflectorP%d",pmtSubI),tyvekPMTReflectorP,
                                                         new G4Tubs(Form("pmtMuMetalSolid_%d",pmtSubI),0.0,pmtMinorAxis/2.0,
                                                                    pmtHousingDepth,0.0,360.0*deg),rotB,thisPMTHousingPlacement);
            
            pmtSubI++;
            iCopy++;
            
            //save coordinates of the front faces of the PMTs
            hPanel[0]->SetBinContent(iCopy,thisPMTPlacement.x());
            hPanel[1]->SetBinContent(iCopy,thisPMTPlacement.y());
            hPanel[2]->SetBinContent(iCopy,thisPMTPlacement.z()-pmtMinorAxis/2.0);
        }
    }

    //y-z plane, -x
    G4RotationMatrix *rotC = new G4RotationMatrix();
    rotC->rotateY(-90.0*deg);
    for(G4int iCopyZ = 0; iCopyZ<numInZ; iCopyZ++){
        for(G4int iCopyY = 0; iCopyY<numInY; iCopyY++){
            
            G4double theY = -spacingInY*numInY/2.0 + double(iCopyY+0.5)*spacingInY;
            G4double theShift = 0.25;
            if(iCopyY%2==0) theShift = -0.25;
            
            G4double theZ = spacingInZ*numInZ/2.0 - (double(iCopyZ+0.5)*spacingInZ + theShift*spacingInZ) + 500.0*mm; // - 500.0*mm
            G4ThreeVector thisPMTPlacement(-(h2oLength/2.0-pmtHousingDepth),theY,theZ);
            new G4PVPlacement(rotC,thisPMTPlacement,pmtLogV,Form("pmtPhysV_%d",iCopy),h2oLogV,true,iCopy,true);
            
            G4ThreeVector thisPMTHousingPlacement(-(h2oLength-pmtHousingDepth)/2.0,theY,theZ);
            new G4PVPlacement(rotC,thisPMTHousingPlacement,pmtHousingLogV,Form("pmtHousingPhysV_%d",iCopy),h2oLogV,true,iCopy,true);
            
            //place holes in reflector for PMT
            tyvekPMTReflectorP = new G4SubtractionSolid(Form("tyvekPMTReflectorP%d",pmtSubI),tyvekPMTReflectorP,
                                                         new G4Tubs(Form("pmtMuMetalSolid_%d",pmtSubI),0.0,pmtMinorAxis/2.0,
                                                                    pmtHousingDepth,0.0,360.0*deg),rotC,thisPMTHousingPlacement);
            
            pmtSubI++;
            iCopy++;

            //save coordinates of the front faces of the PMTs
            hPanel[0]->SetBinContent(iCopy,thisPMTPlacement.x()+pmtMinorAxis/2.0);
            hPanel[1]->SetBinContent(iCopy,thisPMTPlacement.y());
            hPanel[2]->SetBinContent(iCopy,thisPMTPlacement.z());
        }
    }

    //x-y plane, -z
    G4RotationMatrix *rotD = new G4RotationMatrix();
//    rotD->rotateY(0.0*deg);
    //if(false)
    for(G4int iCopyX = 0; iCopyX<numInX; iCopyX++){
        for(G4int iCopyY = 0; iCopyY<numInY; iCopyY++){
            
            G4double theY = -spacingInY*numInY/2.0 + double(iCopyY+0.5)*spacingInY;
            G4double theShift = 0.25;
            if(iCopyY%2==0) theShift = -0.25;
            
            G4double theX = -spacingInX*numInX/2.0 + double(iCopyX+0.5)*spacingInX + theShift*spacingInX;
            G4ThreeVector thisPMTPlacement(theX,theY,-(h2oHeight/2.0-pmtHousingDepth)+500.0*mm); // +500.0*mm
            new G4PVPlacement(rotD,thisPMTPlacement,pmtLogV,Form("pmtPhysV_%d",iCopy),h2oLogV,true,iCopy,true);
            
            G4ThreeVector thisPMTHousingPlacement(theX,theY,-(h2oHeight-pmtHousingDepth)/2.0);
            new G4PVPlacement(rotD,thisPMTHousingPlacement,pmtHousingLogV,Form("pmtHousingPhysV_%d",iCopy),h2oLogV,true,iCopy,true);
            
            //place holes in reflector for PMT
            tyvekPMTReflectorP = new G4SubtractionSolid(Form("tyvekPMTReflectorP%d",pmtSubI),tyvekPMTReflectorP,
                                                         new G4Tubs(Form("pmtMuMetalSolid_%d",pmtSubI),0.0,pmtMinorAxis/2.0,
                                                                    pmtHousingDepth,0.0,360.0*deg),rotD,thisPMTHousingPlacement);
            
            pmtSubI++;
            iCopy++;

            //save coordinates of the front faces of the PMTs
            hPanel[0]->SetBinContent(iCopy,thisPMTPlacement.x());
            hPanel[1]->SetBinContent(iCopy,thisPMTPlacement.y());
            hPanel[2]->SetBinContent(iCopy,thisPMTPlacement.z()+pmtMinorAxis/2.0);
        }
    }
    // Create tyvek reflector with holes for PMTs
    G4LogicalVolume *tyvekPMTReflectorLogV = new G4LogicalVolume(tyvekPMTReflectorP,
                                                                  matPtr->GetMaterial( TYVEK ),
                                                                  "tyvekPMTReflectorLogV");

    G4VPhysicalVolume* tyvekPMTReflectorPhysV = new G4PVPlacement(0,G4ThreeVector(0,0,0),tyvekPMTReflectorLogV,"tyvekPMTReflectorPhysV",h2oLogV,false,0,true);
    
    //D2O into acrylic
    new G4PVPlacement(0,G4ThreeVector(0,0,0),d2oLogV,"d2oPhysV",acrylicLogV,false,0,true);
    //acrylic into H2O
    new G4PVPlacement(0,G4ThreeVector(0,0,0),acrylicLogV,"acrylicPhysV",h2oLogV,false,0,true);
    
    //H2O into full detector volume
    G4VPhysicalVolume *outerVesselPhysV = new G4PVPlacement(0,G4ThreeVector(0,0,0),
                                                            outerVesselLogV,"outerVesselPhysV",
                                                            totalDetLogV,false,0,true);
    
    // Place Muon Vetos
    // Muon Vetos Top
    new G4PVPlacement(0,G4ThreeVector(0,0,+(contOuterHeight+muonVetoThickness)/2.0),
                      muVetoTILogV,"muVetoTIPhysV",
                      totalDetLogV,false,0,true);
    new G4PVPlacement(0,G4ThreeVector(0,0,+(contOuterHeight + muonVetoLayers*muonVetoThickness + muonVetoThickness)/2.0),
                      muVetoTOLogV,"muVetoTOPhysV",
                      totalDetLogV,false,0,true);
    // Muon Vetos Bottom
    new G4PVPlacement(0,G4ThreeVector(0,0,-(contOuterHeight+muonVetoThickness)/2.0),
                      muVetoBILogV,"muVetoBIPhysV",
                      totalDetLogV,false,0,true);
    new G4PVPlacement(0,G4ThreeVector(0,0,-(contOuterHeight + muonVetoLayers*muonVetoThickness + muonVetoThickness)/2.0),
                      muVetoBOLogV,"muVetoBOPhysV",
                      totalDetLogV,false,0,true);
    // Muon Vetos Left
    new G4PVPlacement(0,G4ThreeVector((contOuterLength + muonVetoThickness)/2.0,
                                      0,
                                      -muonVetoThickness),
                      muVetoLILogV,"muVetoLIPhysV",
                      totalDetLogV,false,0,true);
    new G4PVPlacement(0,G4ThreeVector((contOuterLength + muonVetoLayers*muonVetoThickness + muonVetoThickness)/2.0,
                                      0,
                                      -muonVetoThickness),
                      muVetoLOLogV,"muVetoLOPhysV",
                      totalDetLogV,false,0,true);
    // Muon Vetos Right
    new G4PVPlacement(0,G4ThreeVector(-(contOuterLength + muonVetoThickness)/2.0,
                                      0,
                                      -muonVetoThickness),
                      muVetoRILogV,"muVetoRIPhysV",
                      totalDetLogV,false,0,true);
    new G4PVPlacement(0,G4ThreeVector(-(contOuterLength + muonVetoLayers*muonVetoThickness + muonVetoThickness)/2.0,
                                      0,
                                      -muonVetoThickness),
                      muVetoROLogV,"muVetoROPhysV",
                      totalDetLogV,false,0,true);
    // Muon Vetos Far
    new G4PVPlacement(0,G4ThreeVector(0,
                                      -(contOuterWidth + muonVetoThickness)/2.0,
                                      -muonVetoThickness),
                      muVetoFILogV,"muVetoFIPhysV",
                      totalDetLogV,false,0,true);
    new G4PVPlacement(0,G4ThreeVector(0,
                                      -(contOuterWidth + muonVetoLayers*muonVetoThickness + muonVetoThickness)/2.0,
                                      -muonVetoThickness),
                      muVetoFOLogV,"muVetoFOPhysV",
                      totalDetLogV,false,0,true);
    // Muon Vetos Near
    new G4PVPlacement(0,G4ThreeVector(0,
                                      +(contOuterWidth + muonVetoThickness)/2.0,
                                      -muonVetoThickness),
    muVetoNILogV,"muVetoNIPhysV",
    totalDetLogV,false,0,true);
    new G4PVPlacement(0,G4ThreeVector(0,
                                      +(contOuterWidth + muonVetoLayers*muonVetoThickness + muonVetoThickness)/2.0,
                                      -muonVetoThickness),
                      muVetoNOLogV,"muVetoNOPhysV",
                      totalDetLogV,false,0,true);

    // Add Event Handlers for Muon Vetos
    // Muon Vetos Top
    G4d2oSensitiveDetector *senDetmuVetoTO = new G4d2oSensitiveDetector("muVetoTO", muVetoTOLogV->GetMaterial()->GetName(), 0);
    sdManager->AddNewDetector(senDetmuVetoTO);
    muVetoTOLogV->SetSensitiveDetector(senDetmuVetoTO);
    G4d2oSensitiveDetector *senDetmuVetoTI = new G4d2oSensitiveDetector("muVetoTI", muVetoTILogV->GetMaterial()->GetName(), 0);
    sdManager->AddNewDetector(senDetmuVetoTI);
    muVetoTILogV->SetSensitiveDetector(senDetmuVetoTI);
    // Muon Vetos Bottom
    G4d2oSensitiveDetector *senDetmuVetoBO = new G4d2oSensitiveDetector("muVetoBO", muVetoBOLogV->GetMaterial()->GetName(), 0);
    sdManager->AddNewDetector(senDetmuVetoBO);
    muVetoBOLogV->SetSensitiveDetector(senDetmuVetoBO);
    G4d2oSensitiveDetector *senDetmuVetoBI = new G4d2oSensitiveDetector("muVetoBI", muVetoBILogV->GetMaterial()->GetName(), 0);
    sdManager->AddNewDetector(senDetmuVetoBI);
    muVetoBILogV->SetSensitiveDetector(senDetmuVetoBI);
    // Muon Vetos Left
    G4d2oSensitiveDetector *senDetmuVetoLO = new G4d2oSensitiveDetector("muVetoLO", muVetoLOLogV->GetMaterial()->GetName(), 0);
    sdManager->AddNewDetector(senDetmuVetoLO);
    muVetoLOLogV->SetSensitiveDetector(senDetmuVetoLO);
    G4d2oSensitiveDetector *senDetmuVetoLI = new G4d2oSensitiveDetector("muVetoLI", muVetoLILogV->GetMaterial()->GetName(), 0);
    sdManager->AddNewDetector(senDetmuVetoLI);
    muVetoLILogV->SetSensitiveDetector(senDetmuVetoLI);
    // Muon Vetos Right
    G4d2oSensitiveDetector *senDetmuVetoRO = new G4d2oSensitiveDetector("muVetoRO", muVetoROLogV->GetMaterial()->GetName(), 0);
    sdManager->AddNewDetector(senDetmuVetoRO);
    muVetoROLogV->SetSensitiveDetector(senDetmuVetoRO);
    G4d2oSensitiveDetector *senDetmuVetoRI = new G4d2oSensitiveDetector("muVetoRI", muVetoRILogV->GetMaterial()->GetName(), 0);
    sdManager->AddNewDetector(senDetmuVetoRI);
    muVetoRILogV->SetSensitiveDetector(senDetmuVetoRI);
    // Muon Vetos Far
    G4d2oSensitiveDetector *senDetmuVetoFO = new G4d2oSensitiveDetector("muVetoFO", muVetoFOLogV->GetMaterial()->GetName(), 0);
    sdManager->AddNewDetector(senDetmuVetoFO);
    muVetoFOLogV->SetSensitiveDetector(senDetmuVetoFO);
    G4d2oSensitiveDetector *senDetmuVetoFI = new G4d2oSensitiveDetector("muVetoFI", muVetoFILogV->GetMaterial()->GetName(), 0);
    sdManager->AddNewDetector(senDetmuVetoFI);
    muVetoFILogV->SetSensitiveDetector(senDetmuVetoFI);
    // Muon Vetos Near
    G4d2oSensitiveDetector *senDetmuVetoNO = new G4d2oSensitiveDetector("muVetoNO", muVetoNOLogV->GetMaterial()->GetName(), 0);
    sdManager->AddNewDetector(senDetmuVetoNO);
    muVetoNOLogV->SetSensitiveDetector(senDetmuVetoNO);
    G4d2oSensitiveDetector *senDetmuVetoNI = new G4d2oSensitiveDetector("muVetoNI", muVetoNILogV->GetMaterial()->GetName(), 0);
    sdManager->AddNewDetector(senDetmuVetoNI);
    muVetoNILogV->SetSensitiveDetector(senDetmuVetoNI);

    
    G4VPhysicalVolume *h2oPhysV = new G4PVPlacement(0,G4ThreeVector(0,0,0),h2oLogV,"h2oPhysV",outerVesselLogV,false,0,true);

    //Set reflectivity of tyvek Reflector between PMTs
    matPtr->SetReflector(h2oPhysV, tyvekPMTReflectorPhysV, tyvekReflectivity, tyvekSigmaAlpha);
    
    //Now line the sides with either tyvek or area PMT
    G4VPhysicalVolume *sideLiningPhysV_0 = 0;
    G4VPhysicalVolume *sideLiningPhysV_1 = 0;
    
    G4LogicalVolume *theSideLiningLogV = 0;
    if(input->GetSideLining()==0) theSideLiningLogV = GetTyvek();
    else if(input->GetSideLining()==1) theSideLiningLogV = GetAreaPMT();

    if(theSideLiningLogV){
        
        G4RotationMatrix *rotSide = new G4RotationMatrix();
        rotSide->rotateX(180.0*deg);
 
        sideLiningPhysV_0 = new G4PVPlacement(0,G4ThreeVector(0, contOuterWidth/2.0 - outerContainerThickness/2.0 -theSideLiningThickness,0.0),
                                              theSideLiningLogV,"sideLiningPhysV_0",h2oLogV,false,0,true);
        sideLiningPhysV_1 = new G4PVPlacement(rotSide,G4ThreeVector(0,-(contOuterWidth/2.0 - outerContainerThickness/2.0 -theSideLiningThickness),0.0),
                                              theSideLiningLogV,"sideLiningPhysV_1",h2oLogV,false,1,true);
    
        if(input->GetSideLining()==0){ //Set reflectivity of tyvek
            matPtr->SetReflector(h2oPhysV, sideLiningPhysV_0, tyvekReflectivity, tyvekSigmaAlpha);
            matPtr->SetReflector(h2oPhysV, sideLiningPhysV_1, tyvekReflectivity, tyvekSigmaAlpha);
        }
    }
    
    return totalDetLogV;
    
    
}//END of Construct()

G4LogicalVolume *G4d2oDetector::GetAreaPMT(){
    
    ///// Area PMT /////
    G4Box *areaPMTSolid = new G4Box("areaPMTSolid",sideLiningX/2.0,theSideLiningThickness/2.0,sideLiningZ/2.0);
    G4LogicalVolume *areaPMTLogV = new G4LogicalVolume(areaPMTSolid,
                                                       matPtr->GetMaterial( BOROSILICATE ),
                                                       "areaPMTLogV");
    
    G4Box *areaPhotoCSolid = new G4Box("areaPhotoCSolid",sideLiningX/2.0,1.0*mm/2.0,sideLiningZ/2.0);
    G4LogicalVolume *areaPhotoCLogV = new G4LogicalVolume(areaPhotoCSolid,
                                                          matPtr->GetMaterial( PHOTOCATHODE ),
                                                          "areaPhotoCLogV");
    
    G4VisAttributes *visAttPhotoC = new G4VisAttributes();
    visAttPhotoC->SetColour(G4Color::Yellow());
    //    visAttPhotoC->SetForceSolid(true);
    areaPhotoCLogV->SetVisAttributes( visAttPhotoC );
    
    G4cout<<G4endl;
    G4d2oSensitiveDetector *senDetSide = new G4d2oSensitiveDetector("sidePMT", areaPhotoCLogV->GetMaterial()->GetName(), 0);
    sdManager->AddNewDetector(senDetSide);
    areaPhotoCLogV->SetSensitiveDetector(senDetSide);
    G4cout<<G4endl;
    
    new G4PVPlacement(0,G4ThreeVector(0,(theSideLiningThickness-1.0*mm)/2.0,0),areaPhotoCLogV,"areaPhotoCPhysV",areaPMTLogV,false,0,true);
    
    return areaPMTLogV;
    
}
/*
// Taken from cenns10geant4 code
G4LogicalVolume *G4d2oDetector::GetEllipsoidPMT()
{
  G4cout<<"Creating Ellipsoidal PMT\n";
  G4LogicalVolume *pmtLogV = 0;

  // Top ellipsoid of the PMT
  // 76.72mm long, 103.6mm wide
  // Nominal PMT axes are 103.6mm, 76.72mm. Keep that ratio even for alternative pmt diameters
  G4double OPMTTopESemiX      = 103.6*mm;
  G4double OPMTTopESemiY      = OPMTTopESemiX;
  G4double OPMTTopESemiZ      = 76.72*mm;
  G4double OPMTTopEBottomCut  = 0;
  G4double OPMTTopETopCut     = OPMTTopESemiZ;
  G4double PMTFrameThickness  = pmtWindowThickness;

  G4String OPMTTopEVolName    = "PMTTopEVol";
  G4Ellipsoid *OPMTTopESolid = new G4Ellipsoid(OPMTTopEVolName,OPMTTopESemiX, OPMTTopESemiY, OPMTTopESemiZ, OPMTTopEBottomCut, OPMTTopETopCut);

  // Bottom ellipsoid of the PMT
  G4String OPMTBottomEName      = "OPMTBottomEVol";
  G4double OPMTBottomEBottomCut = -OPMTTopETopCut;
  G4double OPMTBottomETopCut    = 0;
 
  G4Ellipsoid *OPMTBottomESolid = new G4Ellipsoid(OPMTBottomEName, OPMTTopESemiX, OPMTTopESemiY, OPMTTopESemiZ, OPMTBottomEBottomCut, OPMTBottomETopCut);

  // Outer Spherical PMT Part
  G4String OPMTSPartName = "OuterPMTSphericalVol";
  G4Transform3D Transform = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,0));
  G4UnionSolid *OPMTSPartSolid  = new G4UnionSolid(OPMTSPartName, OPMTTopESolid, OPMTBottomESolid, Transform);

  // Outer PMT leg
  // Should probably have a check here that the leg doesn't extend out both sides ofthe PMT...
  G4double PMTLegPenetration = OPMTTopESemiZ - 60.84*mm;  // 60.84 - radius value in leg's begining poing
  G4double PMTLegProtruding  = 62*mm;
  G4double PMTLegDiameter    = 84.5*mm;

  G4String OPMTLegName        = "OPMTLegVol";
  G4double OPMTLegInnerRad    = 0;
  G4double OPMTLegOuterRad    = PMTLegDiameter/2;
  G4double OPMTLegZsize       = PMTLegPenetration + PMTLegProtruding;
  G4double OPMTLegSphi        = 0;
  G4double OPMTLegPphi        = 2*M_PI*rad;
  G4double OPMTZposition      = OPMTLegZsize/2 + 60.84*mm;             // Center Position
  G4Tubs *OPMTLegSolid = new G4Tubs(OPMTLegName, OPMTLegInnerRad, OPMTLegOuterRad, OPMTLegZsize/2, OPMTLegSphi, OPMTLegPphi);

  
  // Outer PMT Full Shape
  Transform = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,-OPMTZposition));
  G4String OPMTVolName = "pmtLogV";
  G4UnionSolid      *OPMTSolid  = new G4UnionSolid(OPMTVolName, OPMTSPartSolid, OPMTLegSolid, Transform);
  pmtLogV    = new G4LogicalVolume(OPMTSolid, matPtr->GetMaterial( BOROSILICATE ), OPMTVolName); 

  // Inner PMT shape

  // FIrst the top ellipsoid
  G4String IPMTTopEName         = "IPMTtopEVol";
  G4double IPMTTopESemiX        = OPMTTopESemiX - PMTFrameThickness;
  G4double IPMTTopESemiY        = OPMTTopESemiY - PMTFrameThickness;
  G4double IPMTTopESemiZ        = OPMTTopESemiZ - PMTFrameThickness;
  G4double IPMTTopEBottomCut    = 0;
  G4double IPMTTopECut          = IPMTTopESemiZ;
  G4Ellipsoid *IPMTTopESolid = new G4Ellipsoid(IPMTTopEName,IPMTTopESemiX, IPMTTopESemiY, IPMTTopESemiZ, IPMTTopEBottomCut, IPMTTopECut);

  // THen Inner PMT Bottom Ellipsoid
  G4String IPMTBottomEName      = "IPMTbottomEVol";
  G4double IPMTBottomESemiX     = IPMTTopESemiX;
  G4double IPMTBottomESemiY     = IPMTTopESemiY;
  G4double IPMTBottomESemiZ     = IPMTTopESemiZ;
  G4double IPMTBottomEBottomCut = -IPMTBottomESemiZ;
  G4double IPMTBottomETopCut    = 0;
  G4Ellipsoid *IPMTBottomESolid = new G4Ellipsoid(IPMTBottomEName,IPMTBottomESemiX, IPMTBottomESemiY, IPMTBottomESemiZ, IPMTBottomEBottomCut, IPMTBottomETopCut);

  // Inner PMT Spherical part solid
  G4String  IPMTSPartName = "InnerPMTSphericalVol";
  Transform = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,0));
  G4UnionSolid *IPMTSPartSolid = new G4UnionSolid(IPMTSPartName, IPMTTopESolid,IPMTBottomESolid, Transform);

  // Inner PMT Leg
  G4String IPMTLegName      = "IPMTLegVol";
  G4double IPMTLegInnerRad  = 0;
  G4double IPMTLegOuterRad  = OPMTLegOuterRad - PMTFrameThickness;
  G4double IPMTLegZsize     = OPMTLegZsize - 2*PMTFrameThickness;
  G4double IPMTLegSphi      = 0;
  G4double IPMTLegPphi      = 2*M_PI*rad;
  G4double IPMTZposition    = 60.84*mm  + IPMTLegZsize/2;
  G4Tubs *IPMTLegSolid = new G4Tubs(IPMTLegName, IPMTLegInnerRad, IPMTLegOuterRad, IPMTLegZsize/2, IPMTLegSphi, IPMTLegPphi);
  
  // Inner PMT Full Shape
  G4String IPMTVolName = "InnerPMTVol";
  Transform = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0, 0, -IPMTZposition));
  G4LogicalVolume   *IPMTLog    = new G4LogicalVolume(IPMTSPartSolid, matPtr->GetMaterial( VACUUM ), IPMTVolName);
  IPMTLog->SetVisAttributes(G4VisAttributes::Invisible);
  
  // Then place inner PMT volume inside outer PMT volume
  G4RotationMatrix *IPMTRotation = new G4RotationMatrix();
  IPMTRotation->rotateX(M_PI*rad);
  G4VPhysicalVolume *IPMTPhys      = new G4PVPlacement(0, G4ThreeVector(0,0,0), IPMTLog, IPMTVolName, pmtLogV, false, 0, true); //top

  // Now place the photocathode
  G4String PMTPhtVolName    = "PhotocathodeVol";
  G4double PMTPhtThickness  =  2.6E-5*mm;
  G4double PMTPhtSRad       =  131*mm;                                                // Photocathode Cons Surface Side
  G4double PMTPhtSurDiam    =  190*mm;                                                // Photocathode surface diameter
  G4double PMTPhtAngle      =  acos(PMTPhtSurDiam/(2*PMTPhtSRad))*rad;
  G4double PMTPhtOuterRad   =  PMTPhtSRad-PMTFrameThickness;
  G4double PMTPhtInnerRad   =  PMTPhtSRad - PMTFrameThickness - PMTPhtThickness;
//  G4double PMTPhtZpozition  =  PMTPhtSRad - OPMTTopESemiZ + PMTFrameThickness + PMTPhtThickness/2;
  G4double PMTPhtZpozition  =  - PMTPhtSRad + OPMTTopESemiZ - PMTFrameThickness;
  
  G4Sphere          *PMTPhtSolid  = new G4Sphere(PMTPhtVolName, PMTPhtInnerRad, PMTPhtOuterRad, 0, 2*M_PI*rad, 0, PMTPhtAngle);
  G4LogicalVolume   *PMTPhtLog    = new G4LogicalVolume(PMTPhtSolid, matPtr->GetMaterial( PHOTOCATHODE ), PMTPhtVolName);
  
  G4String PMTPhtTopVolName = "PhotocathodeVol";
  G4VPhysicalVolume *PMTPhtPhys      = new G4PVPlacement(0, G4ThreeVector(0, 0, PMTPhtZpozition), PMTPhtLog, PMTPhtTopVolName, IPMTLog, false, 0,  true);  //top

  ////// Make it sensitive ////////
  G4cout<<G4endl;
  G4d2oSensitiveDetector *senDet = new G4d2oSensitiveDetector("pmt", PMTPhtLog->GetMaterial()->GetName(), 0);
  sdManager->AddNewDetector(senDet);
  PMTPhtLog->SetSensitiveDetector(senDet);
  G4cout<<G4endl;

  G4VisAttributes *visAttPhotoC = new G4VisAttributes();
  visAttPhotoC->SetColour(G4Color::Yellow());
  //  visAttPhotoC->SetForceSolid(true);
  PMTPhtLog->SetVisAttributes( visAttPhotoC );
  pmtLogV->SetVisAttributes( visAttPhotoC );

  return pmtLogV;
}*/

  G4LogicalVolume *G4d2oDetector::GetEllipsoidPMT(){
     ///// PMTs /////
    //--------Ellipsoids--------
    
    // Taken from cenns10geant4 code
  G4cout<<"Creating Ellipsoidal PMT\n";
  G4LogicalVolume *pmtLogV = 0;

  // Top ellipsoid of the PMT
  // 76.72mm long, 103.6mm wide
  // Nominal PMT axes are 103.6mm, 76.72mm. Keep that ratio even for alternative pmt diameters
  G4double OPMTTopESemiX      = 101.0*mm; // was 103.6*mm
  G4double OPMTTopESemiY      = OPMTTopESemiX;
  G4double OPMTTopESemiZ      = 76.72*mm;
  G4double OPMTTopEBottomCut  = 0;
  G4double OPMTTopETopCut     = OPMTTopESemiZ;
  G4double PMTFrameThickness  = pmtWindowThickness;

  G4String OPMTTopEVolName    = "PMTTopEVol";
  G4Ellipsoid *OPMTTopESolid = new G4Ellipsoid(OPMTTopEVolName,OPMTTopESemiX, OPMTTopESemiY, OPMTTopESemiZ, OPMTTopEBottomCut, OPMTTopETopCut);

  // Bottom ellipsoid of the PMT
  G4String OPMTBottomEName      = "OPMTBottomEVol";
  G4double OPMTBottomEBottomCut = -OPMTTopETopCut;
  G4double OPMTBottomETopCut    = 0;
  
 
  G4Ellipsoid *OPMTBottomESolid = new G4Ellipsoid(OPMTBottomEName, OPMTTopESemiX, OPMTTopESemiY, OPMTTopESemiZ, OPMTBottomEBottomCut, OPMTBottomETopCut);

  // Outer Spherical PMT Part
  G4String OPMTSPartName = "OuterPMTSphericalVol";
  G4Transform3D Transform = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,0));
  
  G4UnionSolid *OPMTSPartSolid  = new G4UnionSolid(OPMTSPartName, OPMTTopESolid, OPMTBottomESolid, Transform); //joining top with bottom OUTER ellipsoids.

  // Outer PMT leg
  // Should probably have a check here that the leg doesn't extend out both sides ofthe PMT...
  G4double PMTLegPenetration = OPMTTopESemiZ - 60.84*mm;  // 60.84 - radius value in leg's begining poing
  G4double PMTLegProtruding  = 62*mm; // was 62*mm
  G4double PMTLegDiameter    = 84.5*mm;

  G4String OPMTLegName        = "OPMTLegVol";
  G4double OPMTLegInnerRad    = 0;
  G4double OPMTLegOuterRad    = PMTLegDiameter/2;
  //  G4double OPMTLegZsize       = PMTLegPenetration + PMTLegProtruding;
  G4double OPMTLegZsize       = 147.0*mm;
  G4double OPMTLegSphi        = 0;
  G4double OPMTLegPphi        = 2*M_PI*rad;
  G4double OPMTZposition      = OPMTLegZsize/2 + PMTLegProtruding;             // Center Position
  //  G4double OPMTZposition      = OPMTLegZsize/2 + 60.84*mm;             // Center Position
  // messing with the leg
  
  G4Tubs *OPMTLegSolid = new G4Tubs(OPMTLegName, OPMTLegInnerRad, OPMTLegOuterRad, OPMTLegZsize/2, OPMTLegSphi, OPMTLegPphi);
  //  G4Tubs *OPMTNewSolid = new G4Tubs("OPMTNewVol", 0*mm, 0.1*mm, (147.0*mm)/2, OPMTLegSphi, OPMTLegPphi);
  //  G4Tubs *OPMTLegSolid = new G4Tubs(OPMTLegName, OPMTLegInnerRad, OPMTLegOuterRad, OPMTLegZsize/2, OPMTLegSphi, OPMTLegPphi);

  // end messing with the leg
  // Outer PMT Full Shape
  
  Transform = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,-OPMTZposition));
  
  G4String OPMTVolName = "pmtLogV";
  //  G4UnionSolid      *OPMTIntermediate  = new G4UnionSolid("Intermediate", OPMTSPartSolid, OPMTLegSolid, Transform); //Joining OUTER PMT head with leg
  G4UnionSolid      *OPMTSolid  = new G4UnionSolid(OPMTVolName, OPMTSPartSolid, OPMTLegSolid, Transform); //Joining OUTER PMT head with leg
  pmtLogV    = new G4LogicalVolume(OPMTSolid, matPtr->GetMaterial( BOROSILICATE ), OPMTVolName); //Creating full PMT logical volume, made of Borosilicate.  

  // Inner PMT shape

  // FIrst the TOP ellipsoid
  G4String IPMTTopEName         = "IPMTtopEVol";
  G4double IPMTTopESemiX        = OPMTTopESemiX - PMTFrameThickness;
  G4double IPMTTopESemiY        = OPMTTopESemiY - PMTFrameThickness;
  G4double IPMTTopESemiZ        = OPMTTopESemiZ - PMTFrameThickness;
  G4double IPMTTopEBottomCut    = 0;
  G4double IPMTTopECut          = IPMTTopESemiZ;
  G4Ellipsoid *IPMTTopESolid = new G4Ellipsoid(IPMTTopEName,IPMTTopESemiX, IPMTTopESemiY, IPMTTopESemiZ, IPMTTopEBottomCut, IPMTTopECut); //Top Inner ellipsoid (VACUUM)
  
  G4LogicalVolume *ITopPMTLog= new G4LogicalVolume(IPMTTopESolid, matPtr->GetMaterial( VACUUM ), "ITopPMTLog"); //This is where the photocathode will be placed into.
  
//----------
  
  
  //---Inner BOTTOM PMT Spherical part solid (silver coating to prevent photons to escape)
  
  // Inner PMT Bottom Ellipsoid
  G4String IPMTBottomEName      = "IPMTbottomSilverEVol";
  G4double IPMTBottomESemiX     = IPMTTopESemiX;
  G4double IPMTBottomESemiY     = IPMTTopESemiY;
  G4double IPMTBottomESemiZ     = IPMTTopESemiZ;
  G4double IPMTBottomEBottomCut = -IPMTBottomESemiZ;
  G4double IPMTBottomETopCut    = 0;
  G4Ellipsoid *IPMTBottomESolid = new G4Ellipsoid(IPMTBottomEName,IPMTBottomESemiX, IPMTBottomESemiY, IPMTBottomESemiZ, IPMTBottomEBottomCut, IPMTBottomETopCut); //Bottom Inner ellipsoid
  
  // Inner PMT Leg
  G4String IPMTLegName      = "IPMTSilverLegVol";
  G4double IPMTLegInnerRad  = 0;
  G4double IPMTLegOuterRad  = OPMTLegOuterRad - PMTFrameThickness;
  G4double IPMTLegZsize     = OPMTLegZsize - PMTFrameThickness; //****** 2*
  G4double IPMTLegSphi      = 0;
  G4double IPMTLegPphi      = 2*M_PI*rad;
  G4double IPMTZposition    = 60.84*mm  + IPMTLegZsize/2;
  G4Tubs *IPMTLegSolid = new G4Tubs(IPMTLegName, IPMTLegInnerRad, IPMTLegOuterRad, IPMTLegZsize/2, IPMTLegSphi, IPMTLegPphi); //Inner leg

  
  
  //Joining Inner Bottom and Leg 
  G4String  IPMTSilverPartName = "InnerPMTSilverVol";
  Transform = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,-IPMTZposition));
  
  G4UnionSolid *IPMTSPartSolid = new G4UnionSolid(IPMTSilverPartName, IPMTBottomESolid, IPMTLegSolid, Transform); 


  
  // Inner Bottom SILVER shape:
  
  G4String IPMTSilverVolName = "InnerSilverPMTVol";
  G4LogicalVolume *IPMTSilverLog = new G4LogicalVolume(IPMTSPartSolid, matPtr->GetMaterial( ALUMINUM ), IPMTSilverVolName); //Creating logical volume of INNER HEAD and make it vacuum
  
  
  //IPMTSilverLog->SetVisAttributes(G4VisAttributes::Invisible);
  
  //!!new:
  G4VisAttributes *visAttIPMT = new G4VisAttributes();
  visAttIPMT->SetColour(G4Color::Green());
  //  visAttIPMT->SetForceSolid(true);
  IPMTSilverLog->SetVisAttributes(visAttIPMT);
  
  
  //---BOTTOM Inner made of VACUUM:
  
    // Inner PMT Bottom Ellipsoid
  G4String IPMTBottomEVName      = "IPMTbottomVacuumEVol";
  G4double IPMTBottomEVSemiX     = IPMTTopESemiX - PMTFrameThickness/2.0;
  G4double IPMTBottomEVSemiY     = IPMTTopESemiY - PMTFrameThickness/2.0;
  G4double IPMTBottomEVSemiZ     = IPMTTopESemiZ - PMTFrameThickness/2.0;
  G4double IPMTBottomEVBottomCut = -IPMTBottomEVSemiZ;
  G4double IPMTBottomEVTopCut    = 0;
  G4Ellipsoid *IPMTBottomVESolid = new G4Ellipsoid(IPMTBottomEVName,IPMTBottomEVSemiX, IPMTBottomEVSemiY, IPMTBottomEVSemiZ, IPMTBottomEVBottomCut, IPMTBottomEVTopCut); //Bottom Inner ellipsoid
  
  // Inner PMT Leg
  G4String IPMTLegVName      = "IPMTVacuumLegVol";
  G4double IPMTLegVInnerRad  = 0;
  G4double IPMTLegVOuterRad  = OPMTLegOuterRad - 2.0*PMTFrameThickness;
  G4double IPMTLegVZsize     = OPMTLegZsize - 2.0*PMTFrameThickness;
  G4double IPMTLegSVphi      = 0;
  G4double IPMTLegPVphi      = 2*M_PI*rad;
  G4double IPMTVZposition    = 60.84*mm  + IPMTLegVZsize/2;
  G4Tubs *IPMTVLegSolid = new G4Tubs(IPMTLegVName, IPMTLegVInnerRad, IPMTLegVOuterRad, IPMTLegVZsize/2, IPMTLegSVphi, IPMTLegPVphi); //Inner leg

  
  //Joining Inner Bottom and Leg 
  G4String  IPMTSVPartName = "InnerPMTVacuumVol";
  Transform = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,-IPMTVZposition));
  //Transform = G4Transform3D(*PMTRotation, G4ThreeVector(0,0,0));//*New
  G4UnionSolid *IPMTSVPartSolid = new G4UnionSolid(IPMTSVPartName, IPMTBottomVESolid, IPMTVLegSolid, Transform); 

  
  // Inner BOTTOM VACUUM Shape
  
  G4String IPMTVolName = "InnerVacuumPMTVol";
  
  G4LogicalVolume   *IPMTVacuumLog    = new G4LogicalVolume(IPMTSVPartSolid, matPtr->GetMaterial( VACUUM ), IPMTVolName); //Creating logical volume of INNER HEAD and make it vacuum
  
  //IPMTVacuumLog->SetVisAttributes(G4VisAttributes::Invisible);
  
  //!!new:
  G4VisAttributes *visAttIVPMT = new G4VisAttributes();
  visAttIVPMT->SetColour(G4Color::Magenta());
  //  visAttIVPMT->SetForceSolid(true);
  IPMTVacuumLog->SetVisAttributes(visAttIVPMT);
  
  
  
  
  //--- PLACE VACUUM BOTTOM PART INTO SILVER BOTTOM PART 
  
  G4String  BottomPMTPartName = "BottomPart";
  
  G4VPhysicalVolume *BottomPMT  = 
  new G4PVPlacement(0, G4ThreeVector(0,0,0), IPMTVacuumLog, BottomPMTPartName, IPMTSilverLog, false, 0, true); // "Vacuum in Silver"
  
  // SILVER IN OUTER pmtLogV
  G4VPhysicalVolume *AgCoating = 
  new G4PVPlacement(0, G4ThreeVector(0,0,0), IPMTSilverLog, "BottomTotal", pmtLogV, false, 0, true);
  
  G4Tubs *ThickRingSolid = new G4Tubs("ThickRingSolid", 85.0*mm/2.0, 116.0*mm/2.0, 85*mm/2.0, IPMTLegSVphi, IPMTLegPVphi); //Inner legq
  G4LogicalVolume *ThickRingLogV = new G4LogicalVolume(ThickRingSolid, matPtr->GetMaterial( EPOXY ), "ThickRingLogV"); //Creating logical volume of INNER HEAD and make it vacuum
  G4VisAttributes *visThickRing = new G4VisAttributes();
  visThickRing->SetColour(G4Color::Blue());
  ThickRingLogV->SetVisAttributes(visThickRing);
  new G4PVPlacement(0, G4ThreeVector(0,0,-158.0*mm/2.0-147.0*mm+85.0*mm/2.0), ThickRingLogV, "ThickRing", pmtLogV, false, 0, true);
  
  G4Tubs *ThinRingSolid = new G4Tubs("ThinRingSolid", 85.0*mm/2.0, 95.0*mm/2.0, 62.0*mm/2.0, IPMTLegSVphi, IPMTLegPVphi); //Inner legq
  G4LogicalVolume *ThinRingLogV = new G4LogicalVolume(ThinRingSolid, matPtr->GetMaterial( EPOXY ), "ThinRingLogV"); //Creating logical volume of INNER HEAD and make it vacuum
  G4VisAttributes *visThinRing = new G4VisAttributes();
  visThinRing->SetColour(G4Color::Red());
  ThinRingLogV->SetVisAttributes(visThinRing);
  new G4PVPlacement(0, G4ThreeVector(0,0,-158.0*mm/2.0-62.0*mm/2.0), ThinRingLogV, "ThinRing", pmtLogV, false, 0, true);

  G4Tubs *HolderRingSolid = new G4Tubs("HolderRingSolid", 96.0*mm/2.0, 116.0*mm/2.0, 20.0*mm/2.0, IPMTLegSVphi, IPMTLegPVphi); //Inner legq
  G4LogicalVolume *HolderRingLogV = new G4LogicalVolume(HolderRingSolid, matPtr->GetMaterial( DELRIN ), "HolderRingLogV"); //Creating logical volume of INNER HEAD and make it vacuum
  G4VisAttributes *visHolderRing = new G4VisAttributes();
  visHolderRing->SetColour(G4Color::Yellow());
  HolderRingLogV->SetVisAttributes(visHolderRing);
  new G4PVPlacement(0, G4ThreeVector(0,0,-158.0*mm/2.0-20.0*mm/2.0), HolderRingLogV, "HolderRing", pmtLogV, false, 0, true);  
  
  G4Tubs *HolderRodSolid = new G4Tubs("HolderRodSolid", 0.0*mm/2.0, 20.0*mm/2.0, 85.0*mm/2.0, IPMTLegSVphi, IPMTLegPVphi); //Inner legq
  G4LogicalVolume *HolderRodLogV = new G4LogicalVolume(HolderRodSolid, matPtr->GetMaterial( DELRIN ), "HolderRodLogV"); //Creating logical volume of INNER HEAD and make it vacuum
  G4VisAttributes *visHolderRod = new G4VisAttributes();
  visHolderRod->SetColour(G4Color::Brown());
  HolderRodLogV->SetVisAttributes(visHolderRod);
  new G4PVPlacement(0, G4ThreeVector((116.0+20.0)*mm/2.0*sin(45*deg),(116.0+20.0)*mm/2.0*cos(45*deg),-158.0*mm/2.0-147.0*mm+85.0*mm/2.0), HolderRodLogV, "HolderRod", pmtLogV, false, 0, true);
  new G4PVPlacement(0, G4ThreeVector((116.0+20.0)*mm/2.0*sin(45+90*deg),(116.0+20.0)*mm/2.0*cos(45+90*deg),-158.0*mm/2.0-147.0*mm+85.0*mm/2.0), HolderRodLogV, "HolderRod", pmtLogV, false, 0, true);
  new G4PVPlacement(0, G4ThreeVector((116.0+20.0)*mm/2.0*sin(45+180*deg),(116.0+20.0)*mm/2.0*cos(45+180*deg),-158.0*mm/2.0-147.0*mm+85.0*mm/2.0), HolderRodLogV, "HolderRod", pmtLogV, false, 0, true);
  new G4PVPlacement(0, G4ThreeVector((116.0+20.0)*mm/2.0*sin(45+270*deg),(116.0+20.0)*mm/2.0*cos(45+270*deg),-158.0*mm/2.0-147.0*mm+85.0*mm/2.0), HolderRodLogV, "HolderRod", pmtLogV, false, 0, true);  

  
//////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////
//   Optical properties of the interface between the Air and Reflective Surface
//   For Mirror, reflectivity is set at 100% and specular reflection is assumed.
  
     
    /*G4OpticalSurface* OpticalAirMirror = new G4OpticalSurface("OpticalAirMirror");
    
    new G4LogicalBorderSurface("OpticalAirMirror", AgCoating , BottomPMT, OpticalAirMirror );
    
    OpticalAirMirror->SetModel(unified);
    OpticalAirMirror->SetType(dielectric_dielectric);
    OpticalAirMirror->SetFinish(polished);
    
    const G4int NUM = 3;
    G4double XX[NUM] = {1.5*eV, 4.25*eV, 7.0*eV} ; 
    G4double ICEREFLECTIVITY[NUM]  = { 1.0, 1.0, 1.0 };
    
     G4MaterialPropertiesTable *AirMirrorMPT = new G4MaterialPropertiesTable();
     AirMirrorMPT->AddProperty("REFLECTIVITY", XX, ICEREFLECTIVITY,NUM);
     OpticalAirMirror->SetMaterialPropertiesTable(AirMirrorMPT);*/
  
    matPtr->SetReflector(BottomPMT, AgCoating, 1.0, 0.0);
     
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////     
     
  
  // TOP PART IN OUTER pmtLogV
  
  //G4VPhysicalVolume *TopPMT = 
  new G4PVPlacement(0, G4ThreeVector(0,0,0), ITopPMTLog, "TopTotal", pmtLogV, false, 0, true);
  
  
  // Then place inner PMT volume inside outer PMT volume
  //G4VPhysicalVolume *IPMTPhys  = new G4PVPlacement(0, G4ThreeVector(0,0,0), IPMTLog, IPMTVolName, pmtLogV, false, 0, true); //top
 

  // Now place the photocathode
  
  G4String PMTPhtVolName    = "PhotocathodeVol";
  G4double PMTPhtThickness  =  2.6E-5*mm;
  G4double PMTPhtSRad       =  131*mm;                                                // Photocathode Cons Surface Side
  G4double PMTPhtSurDiam    =  190*mm;                                                // Photocathode surface diameter
  G4double PMTPhtAngle      =  acos(PMTPhtSurDiam/(2*PMTPhtSRad))*rad;
  G4double PMTPhtOuterRad   =  PMTPhtSRad-PMTFrameThickness;
  G4double PMTPhtInnerRad   =  PMTPhtSRad - PMTFrameThickness - PMTPhtThickness;
//  G4double PMTPhtZpozition  =  PMTPhtSRad - OPMTTopESemiZ + PMTFrameThickness + PMTPhtThickness/2;
  G4double PMTPhtZpozition  =  - PMTPhtSRad + OPMTTopESemiZ - PMTFrameThickness;
  
  G4Sphere          *PMTPhtSolid  = new G4Sphere(PMTPhtVolName, PMTPhtInnerRad, PMTPhtOuterRad, 0, 2*M_PI*rad, 0, PMTPhtAngle);
  G4LogicalVolume   *PMTPhtLog    = new G4LogicalVolume(PMTPhtSolid, matPtr->GetMaterial( PHOTOCATHODE ), PMTPhtVolName);
  
  G4String PMTPhtTopVolName = "PhotocathodeVol";
  G4VPhysicalVolume *PMTPhtPhys  = new G4PVPlacement(0, G4ThreeVector(0, 0, PMTPhtZpozition), PMTPhtLog, PMTPhtTopVolName, ITopPMTLog, false, 0,  true);  //top
  

  ////// Make it sensitive ////////
  G4cout<<G4endl;
  G4d2oSensitiveDetector *senDet = new G4d2oSensitiveDetector("pmt", PMTPhtLog->GetMaterial()->GetName(), 0);
  G4SDManager *sdManager = G4SDManager::GetSDMpointer(); //New
  sdManager->AddNewDetector(senDet);
  PMTPhtLog->SetSensitiveDetector(senDet);
  G4cout<<G4endl;

  G4VisAttributes *visAttPhotoC = new G4VisAttributes();
  visAttPhotoC->SetColour(G4Color::Yellow());
  //  visAttPhotoC->SetForceSolid(true);
  
  G4VisAttributes *visAttPhotoM = new G4VisAttributes();
  visAttPhotoM->SetColour(G4Color::White());
  //  visAttPhotoM->SetForceSolid(true);
  
  PMTPhtLog->SetVisAttributes( visAttPhotoC );
  pmtLogV->SetVisAttributes( visAttPhotoM );

  return pmtLogV;
}

G4LogicalVolume *G4d2oDetector::GetSphericalPMT(){
    
    ///// PMTs /////
    //Half-spheres for now
    G4Sphere *pmtSolid = new G4Sphere("pmtSolid",0.0,pmtDiameter/2.0,0.0,360.0*deg,0.0,90.0*deg);
    G4LogicalVolume *pmtLogV = new G4LogicalVolume(pmtSolid,
                                                   matPtr->GetMaterial( VACUUM ),
                                                   "pmtLogV");
   // pmtLogV->SetVisAttributes(G4VisAttributes::Invisible);
    pmtLogV->SetVisAttributes(G4VisAttributes::GetInvisible()); //new for running in topaz!!
    
    //use a borosilicate glass window
    G4Sphere *pmtGlassSolid = new G4Sphere("pmtGlassSolid",pmtDiameter/2.0-pmtWindowThickness,pmtDiameter/2.0,
                                           0.0,360.0*deg,0.0,90.0*deg);
    G4LogicalVolume *pmtGlassLogV = new G4LogicalVolume(pmtGlassSolid,
                                                        matPtr->GetMaterial( BOROSILICATE ),
                                                        "pmtGlassLogV");
    new G4PVPlacement(0,G4ThreeVector(0,0,0),pmtGlassLogV,"pmtGlassPhysV",pmtLogV,false,0,true);
    //for now, photocathode layer that just absorbs whatever passes through glass
    // photocathod thickness from Nuclear Instruments and Methods in Physics Research A 539 (2005) 217–235
    // TODO: Adjust photocathode  be less than 90 degrees
    G4Sphere *pmtPhotoCathodeSolid = new G4Sphere("pmtPhotoCathodeSolid",
                                                  pmtDiameter/2.0-pmtWindowThickness-20.0e-9*m,
                                                  pmtDiameter/2.0-pmtWindowThickness,
                                                  0.0,360.0*deg,0.0,90.0*deg);
    G4LogicalVolume *pmtPhotoCathodeLogV = new G4LogicalVolume(pmtPhotoCathodeSolid,
                                                               matPtr->GetMaterial( PHOTOCATHODE ),
                                                               "pmtPhotoCathodeLogV");
    new G4PVPlacement(0,G4ThreeVector(0,0,0),pmtPhotoCathodeLogV,"pmtPhotoCathodePhysV",pmtLogV,false,0,true);
    G4VisAttributes *visAttPhotoC = new G4VisAttributes();
    visAttPhotoC->SetColour(G4Color::Yellow());
    //    visAttPhotoC->SetForceSolid(true);
    pmtPhotoCathodeLogV->SetVisAttributes( visAttPhotoC );
    pmtGlassLogV->SetVisAttributes( visAttPhotoC );
    
    ////// Make it sensitive ////////
    G4cout<<G4endl;
    G4d2oSensitiveDetector *senDet = new G4d2oSensitiveDetector("pmt", pmtPhotoCathodeLogV->GetMaterial()->GetName(), 0);
    sdManager->AddNewDetector(senDet);
    pmtPhotoCathodeLogV->SetSensitiveDetector(senDet);
    G4cout<<G4endl;

    return pmtLogV; // not this one
    
}

G4LogicalVolume *G4d2oDetector::GetTyvek(){
    
    ///// Tyvek reflector /////
    G4Box *tyvekSolid = new G4Box("d2oSolid",sideLiningX/2.0,theSideLiningThickness/2.0,sideLiningZ/2.0);
    G4LogicalVolume *tyvekLogV = new G4LogicalVolume(tyvekSolid,
                                                      matPtr->GetMaterial( TYVEK ),
                                                      "tyvekLogV");
    
    return tyvekLogV;
    
}

G4ThreeVector G4d2oDetector::GetPMTPosition(Int_t iPMT){
    
    if(!hPanel[0] || !hPanel[1] || !hPanel[2]){
        printf("Error in G4d2oDetector::GetPMTPosition() - pixel centers not defined.\n\n");
        exit(0);
    }
    
    G4ThreeVector thePosition(hPanel[0]->GetBinContent(iPMT+1),
                              hPanel[1]->GetBinContent(iPMT+1),
                              hPanel[2]->GetBinContent(iPMT+1));
    
    return thePosition;
  
}

G4int G4d2oDetector::InitializeQEArray() {
    const double hc_evnm = 1.23984193 *1e3;
    std::vector<double> qe_bialkali_y;
    std::vector<double> qe_bialkali_xgev;
    std::vector<double> qe_bialkali_xnm;

    int pmt_qe_mapping = input->GetPMTQE();
  
    if (!pmt_qe_mapping) {
      qe_bialkali_y = {0.000,0.000,0.000,0.000,
        0.001,0.003,0.008,0.017,
        0.033,0.052,0.076,0.106,
        0.139,0.164,0.193,0.204,
        0.228,0.240,0.254,0.254,
        0.240,0.213,0.164,0.032,
        0.003,0.0};

      qe_bialkali_xgev = {1.63e-9,1.68e-9,1.72e-9,1.77e-9,
        1.82e-9,1.88e-9,1.94e-9,2.00e-9,
        2.07e-9,2.14e-9,2.21e-9,2.30e-9,
        2.38e-9,2.48e-9,2.58e-9,2.70e-9,
        2.82e-9,2.95e-9,3.10e-9,3.26e-9,
        3.44e-9,3.65e-9,3.88e-9,4.13e-9,
        4.43e-9,4.7e-9};
    } else {
      qe_bialkali_xgev = {
          1.71510330e-09, 1.72505427e-09, 1.73613456e-09, 1.74735811e-09,
          1.75872773e-09, 1.75872773e-09, 1.77101314e-09, 1.78298526e-09,
          1.79590892e-09, 1.79504152e-09, 1.80822117e-09, 1.82232657e-09,
          1.83696324e-09, 1.85299105e-09, 1.86812650e-09, 1.87165448e-09,
          1.88829612e-09, 1.88973633e-09, 1.90912889e-09, 1.93247960e-09,
          1.95922137e-09, 1.98804205e-09, 2.02183955e-09, 2.06251326e-09,
          2.10933689e-09, 2.14278094e-09, 2.16147404e-09, 2.21625389e-09,
          2.26867743e-09, 2.31459658e-09, 2.36805769e-09, 2.42800588e-09,
          2.45405844e-09, 2.52277310e-09, 2.55091084e-09, 2.62061070e-09,
          2.65098636e-09, 2.73135221e-09, 2.76436551e-09, 2.85734752e-09,
          2.89349699e-09, 3.00157732e-09, 3.04149380e-09, 3.16114168e-09,
          3.20544622e-09, 3.33862347e-09, 3.38808152e-09, 3.54565660e-09,
          3.60149007e-09, 3.75618664e-09, 3.81400721e-09, 3.97195740e-09,
          4.07641626e-09, 4.14520577e-09, 4.19263608e-09, 4.23805531e-09,
          4.27052054e-09, 4.29313044e-09, 4.32435077e-09, 4.42825911e-09,
          4.48977939e-09, 4.54931759e-09, 4.57029657e-09, 4.59857138e-09
      };

      qe_bialkali_y = {
          1.21916616e-04, 1.64916271e-04, 2.25148378e-04, 3.03070686e-04,
          3.89609228e-04, 4.19063429e-04, 5.60020355e-04, 7.57130441e-04,
          9.66787262e-04, 1.02733496e-03, 1.38190893e-03, 1.87533797e-03,
          2.53999850e-03, 3.42553496e-03, 4.35274392e-03, 4.65981070e-03,
          5.86334223e-03, 6.28450848e-03, 8.35604877e-03, 1.11955399e-02,
          1.50356257e-02, 1.99577242e-02, 2.60925624e-02, 3.38645791e-02,
          4.33569287e-02, 5.27538911e-02, 5.46512894e-02, 6.84281357e-02,
          8.84422738e-02, 1.15409924e-01, 1.49302391e-01, 1.74484430e-01,
          1.78010685e-01, 1.98627112e-01, 2.03661354e-01, 2.31722450e-01,
          2.37462978e-01, 2.65408769e-01, 2.70018275e-01, 2.96629739e-01,
          3.00940803e-01, 3.22502973e-01, 3.26551781e-01, 3.38995211e-01,
          3.39631673e-01, 3.43448449e-01, 3.43517826e-01, 3.39532867e-01,
          3.37523668e-01, 3.11806132e-01, 3.02439900e-01, 2.39203198e-01,
          1.77953871e-01, 1.30590694e-01, 9.59773034e-02, 6.89963758e-02,
          5.22438279e-02, 4.00873446e-02, 2.95875022e-02, 2.00769923e-02,
          1.08495722e-02, 4.16906330e-03, 3.15345936e-03, 2.30471385e-03
      };
    }

    for(auto xgev : qe_bialkali_xgev) qe_bialkali_xnm.push_back(hc_evnm/(xgev*1e9));
  
    if (!fPMTQE) fPMTQE = new TGraph(qe_bialkali_xgev.size());
  
    for(size_t ipoint = 0; ipoint<qe_bialkali_xnm.size(); ipoint++ ){
        
       fPMTQE->SetPoint(ipoint,qe_bialkali_xgev[ipoint]*1e9,(1.033*1.01/0.8368)*qe_bialkali_y[ipoint]); //Scaled to match DATA and MC! 02/2023

    }
    fminEnergyQE = qe_bialkali_xgev[0]*1.e9;
    fmaxEnergyQE = qe_bialkali_xgev[qe_bialkali_xgev.size() - 1]*1.e9;

    return 0;
}

G4double G4d2oDetector::EvalPMTQE(G4double energy_ev) 
{
  if (energy_ev > fmaxEnergyQE || energy_ev < fminEnergyQE)
    return 0.0;

  return fPMTQE->Eval(energy_ev); // Evaluate the pmt QE for a given wavelength
}
