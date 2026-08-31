// PMTs ON TOP GEOMETRY!!

#include "G4d2oCylindricalDetector.hh"
#include "G4d2oRunAction.hh"
#include "G4d2oSensitiveDetector.hh"
#include "inputVariables.hh"
#include "G4d2oDataDrivenReflector.hh"

#include "TMath.h"

#include "G4PhysicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4VisAttributes.hh"
#include "G4PVParameterised.hh"
#include "G4SubtractionSolid.hh"
#include <iostream>

using namespace std;

G4d2oCylindricalDetector::G4d2oCylindricalDetector()
{
    G4cerr << "\n\tUsing cylindrical geometry" << G4endl ;
    Initialize();
}

void G4d2oCylindricalDetector::Initialize(){

    //relevant parameters
    sink = (2.0*in + 2*0.06*in + 0.5*in)/2.0;
    airthickness = 0.38*in;
    aluminumthickness = 0.06; // in inches
    //D2O
    d2oLength = 26.5*in;//radius // was 70.0*cm
    d2oWidth = d2oLength; //radius
    d2oHeight = 55.12*in - 2*in;//z // WAS 140.0*cm

    //acrylic
    acrylicThickness = 0.5*in; // was 0.25*in
    acrylicEndCapThickness = 1.0*in;

    //pmt
    pmtDepth = 10.0*in;
    muMetalThickness = 0.03*in;
    pmtWindowThickness = 3.0*mm;  //making this up
    minGapBetweenPMTs = 3.0*mm;
    pmtDiameter = input->GetPMTDiameter()*in; //pmtDiameter controlled from beamOn.dat or command line
    pmtMinorAxis = pmtMajorAxis = pmtDiameter;
    pmtLegLength = 62.*mm; // was 62.*mm
    if (pmtDiameter < 0.0) {
      pmtMinorAxis = 76.72*mm * 2.0;
      pmtMajorAxis = 103.6*mm * 2.0;
    }

    //tyvek
    tyvekThickness = 0.25*in;
    tyvekReflectivity = input->GetReflectivity();
    G4cout << "tyvekReflectivity = " << G4endl;
    G4cout << tyvekReflectivity << G4endl;
    tyvekSigmaAlpha = 0.1; //roughness parameter. 0.1 for tyvek, use 0 for polished

    //H2O
    h2oTCthickness = 7.5*in/2;
    h2oLength = d2oLength + 2*(acrylicThickness + h2oTCthickness); //x
    h2oWidth = h2oLength; //y
    if (pmtDiameter > 0.0) {
      //      h2oHeight = d2oHeight + 2*(acrylicEndCapThickness + h2oTCthickness + pmtMinorAxis/2.0); //z // comenting because never happens
    } else {
      bottomThickness = 4.13*in;
      h2oHeight = 75.55*in - 0.5*in;
    }

    //not using these for now, but set them anyway so that other parts of code don't have problems
    h2oInnerLength = h2oLength;
    h2oInnerWidth = h2oWidth;
    h2oInnerHeight = h2oHeight;

    //outer vessel
    outerContainerThickness = 0.25*in;

    ///// Full detector volume /////
    contOuterLength = h2oLength + 2*outerContainerThickness; //x
    contOuterWidth = contOuterLength; //y
    contOuterHeight = h2oHeight + 2*outerContainerThickness; //z

    ///// Shielding volume /////
    shieldThickness = 2*in;// input->GetShieldThickness()*in;
    h2oRefl = input->GetH2oRefl();
    shieldLength = contOuterLength + 2.0*shieldThickness; // x
    shieldWidth = shieldLength;
    shieldHeight = 80.856*in; // 342/9*2*in+.5*in+2*in+2*in+.356*in from engineering drawing

    //Veto Layers
    muonVetoThickness = 2.54*cm;
    muonVetoLayers = 2;

    vetoOuterLength = (44.688+2*0.06)*in; //x
    vetoOuterWidth = vetoOuterLength; //y
    if (!input->GetBottomVeto() ) {
      vetoOuterHeight = shieldHeight + muonVetoLayers*muonVetoThickness;
    }

    iUseBottomPMTs = input->GetBottomPMTs();

    DetermineSpacing();

    G4cerr << "done." << G4endl;

}//END of constructor

G4d2oCylindricalDetector::~G4d2oCylindricalDetector()
{
    G4cout << "Instance of G4d2oCylindricalDetector Destructed!" << G4endl;
}

void G4d2oCylindricalDetector::DetermineSpacing(){
  spacingInX = pmtMajorAxis + (2.0*minGapBetweenPMTs);
  numInX = TMath::FloorNint( (h2oLength / spacingInX) );

  spacingInY = TMath::Sqrt(3.0) * ( (pmtMajorAxis/2.0) + 2.0*minGapBetweenPMTs);
  numInY = TMath::FloorNint( (h2oLength / spacingInY) );
}

G4LogicalVolume * G4d2oCylindricalDetector::GetDetector(){

    //Materials pointer
    matPtr = G4d2oRunAction::GetMaterialsPointer();

    ///// create total detector volume /////
    G4LogicalVolume *totalDetLogV = GetTotalDetectorLogV();

    ///// create lead shielding volume /////
    G4LogicalVolume *shieldingLogV = GetShieldingLogV();

    ///// create steel vessel /////
    G4LogicalVolume *outerVesselLogV = GetOuterVesselLogV();

    ///// create airel vessel /////
    G4LogicalVolume *outerAirselLogV = GetOuterAirselLogV();

    ///// Create Tyvek Lining /////
    G4LogicalVolume *tyvekLining = GetTyvekLiningLogV();

    ///// H2O Inner /////
    G4LogicalVolume *h2oLogV = GetH2OLogV();

    ///// PMMA /////
    G4LogicalVolume *acrylicLogV = GetAcrylicLogV();

    // ============================================================
    // ============================================================
    G4LogicalVolume *targetLogV = GetH2OModuleLogV();
    G4cout << "\n>>> USING MODULE 2 (both H2O) <<<" << G4endl;

    ///// PMTs with mu-metal shield /////
    G4LogicalVolume *pmtLogV = 0;
    if (pmtDiameter < 1.0) {
      pmtLogV = GetEllipsoidPMT();
    } else {
      pmtLogV = GetSphericalPMT();
    }

    ///// InnerVeto Layer
    G4LogicalVolume *innerVetoLogV = GetInnerVetoLogV();

    /////OuterVeto Layer

    G4LogicalVolume *topAluminumLogV = GetTopAluminumLogV();
    G4LogicalVolume *topAirLogV = GetTopAirLogV();
    G4LogicalVolume *basePlateLogV = GetBasePlateLogV();
    G4LogicalVolume *ribsLogV = GetRibsLogV();
    G4LogicalVolume *strawLogV = GetStrawLogV();
    G4LogicalVolume *airStrawLogV = GetAirStrawLogV();
    G4LogicalVolume *leadCapLogV = GetLeadCapLogV();
    G4LogicalVolume *feetLogV = GetFeetLogV();
    G4LogicalVolume *narrowAluminumLogV = GetNarrowAluminumLogV();
    G4LogicalVolume *narrowAirLogV = GetNarrowAirLogV();
    G4LogicalVolume *wideAluminumLogV = GetWideAluminumLogV();
    G4LogicalVolume *wideAirLogV = GetWideAirLogV();
    G4LogicalVolume *topVetoLogV = GetTopVetoLogV();
    G4LogicalVolume *narrowVetoLogV = GetNarrowVetoLogV();
    G4LogicalVolume *wideVetoLogV = GetWideVetoLogV();

    //     Geometry hierarchy
    //     totalDetLogV (mother volume)
    //       -> veto layers
    //        -> outer steel
    //           -> tyvek lining
    //              -> H2O tank
    //                 -(1)-> PMTs
    //                 -(2)-> acrylic tank
    //                        -> target (H2O, Module 2)

    ///// Place the target in acrylic tank /////
    G4VPhysicalVolume *targetPhysV = new G4PVPlacement(0,G4ThreeVector(0,0,0),targetLogV,"targetPhysV",acrylicLogV,false,0,true);

    ///// Place the acrylic tank in H2O /////
    G4VPhysicalVolume *acrylicPhysV = new G4PVPlacement(0,G4ThreeVector(0,0,(+55.12*in-(h2oHeight-0.5*pmtLegLength-0.5*pmtMinorAxis))/2.0+4.13*in+tyvekThickness),acrylicLogV,"acrylicPhysV",h2oLogV,false,0,true);

    ///// Place PMTs into H2O /////
    PlacePMTs(pmtLogV,h2oLogV);

    //// Create tyvek end caps ////
    G4LogicalVolume *tyvekCapTop = GetTyvekLiningCapLogV("topCap", true);
    G4LogicalVolume *tyvekCapBot = GetTyvekLiningCapLogV("botCap", iUseBottomPMTs);

    ///// Place the H2O in outer vessel /////
    G4VPhysicalVolume *h2oPhysV = new G4PVPlacement(0,G4ThreeVector(0,0,-(0.5*pmtLegLength+0.5*pmtMinorAxis)/2.0),h2oLogV,"h2oPhysV",outerAirselLogV,false,0,true);

    ///// Place the tyvek lining in H2O/////
    G4VPhysicalVolume *tyvekPhysV = new G4PVPlacement(0,G4ThreeVector(0,0,contOuterHeight/2.0-(h2oHeight-0.5*pmtLegLength-0.5*pmtMinorAxis)/2.0-outerContainerThickness/2.0),tyvekLining,"tyvekPhysV",h2oLogV,false,0,true);
    auto g4pvs = G4PhysicalVolumeStore::GetInstance();

    // Place the tyvek caps
    G4VPhysicalVolume *tyvekCapBotPhysV = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, -(h2oHeight-tyvekThickness-0.5*pmtLegLength-0.5*pmtMinorAxis)/2.0), tyvekCapBot, "tyvekCapBotPhysV", h2oLogV, false, 0, true);

    // ============================================================
    // REFLECTOR MODEL: Data-Driven 
    // ============================================================
    G4cout << "\n>>> USING DATA-DRIVEN MODEL <<<" << G4endl;
    G4String dataDir = "angular_data";
    G4d2oDataDrivenReflector::CreateInstance(tyvekReflectivity, dataDir);
    matPtr->SetDataDrivenReflector(tyvekPhysV, h2oPhysV, tyvekReflectivity, dataDir);
    matPtr->SetDataDrivenReflector(h2oPhysV, tyvekPhysV, tyvekReflectivity, dataDir);
    matPtr->SetDataDrivenReflector(h2oPhysV, tyvekCapBotPhysV, tyvekReflectivity, dataDir);
    matPtr->SetDataDrivenReflector(tyvekCapBotPhysV, h2oPhysV, tyvekReflectivity, dataDir);

    // setting other reflectivities
    acrylicReflectivity = 0.0263;
    acrylicSigmaAlpha = 0.0;
    matPtr->SetReflector(h2oPhysV, acrylicPhysV, acrylicReflectivity, acrylicSigmaAlpha);
    matPtr->SetReflector(acrylicPhysV, targetPhysV, acrylicReflectivity, acrylicSigmaAlpha);

    ///// Place the outer vessel inside the shielding
    G4double offset = -0.5*shieldThickness;
    std::cout << "the offset is=" << offset << " ";
    G4VPhysicalVolume *outerVesselPhysV = new G4PVPlacement(0,G4ThreeVector(0,0,-contOuterHeight/2.0+shieldHeight/2.0-1.75*in),outerVesselLogV,"outerVesselPhysV",shieldingLogV,false,0,true);
    G4VPhysicalVolume *outerAirselPhysV = new G4PVPlacement(0,G4ThreeVector(0,0,0),outerAirselLogV,"outerAirselPhysV",outerVesselLogV,false,0,true);

    ///// Place the outer vessel in total detector /////

    aluminumthickness = 0.06;
    new G4PVPlacement(0,G4ThreeVector(0,0,86.430*in/2.0+(2.*in+4*0.06*in)/2.0-sink+airthickness/2.0),topAluminumLogV,"topAluminumPhysV",totalDetLogV,false,0,true);
    new G4PVPlacement(0,G4ThreeVector(0,0,-airthickness/2.0),topVetoLogV,"topVetoPhysV",topAluminumLogV,false,0,true);
    new G4PVPlacement(0,G4ThreeVector(0,0,+2.0*in/2.0),topAirLogV,"topAirPhysV",topAluminumLogV,false,0,true);

    new G4PVPlacement(0,G4ThreeVector(0,0,shieldHeight/2.0-contOuterHeight-1.75*in-3.556*in/2.0),basePlateLogV,"basePlatePhysV",shieldingLogV,false,0,true);

    for (int pos=0; pos<4; pos++){
      G4RotationMatrix* Rotationribs = new G4RotationMatrix();
      Rotationribs->rotateX(0*deg);
      Rotationribs->rotateY(0*deg);
      Rotationribs->rotateZ(pos*90*deg);
      std::string ribs_pv_name = "ribsPhysV_" + std::to_string(pos);
      new G4PVPlacement(Rotationribs,G4ThreeVector((27.5/2.0*in+3.375/2.0*in)*sin(pos*90*deg),(27.5/2.0*in+3.375/2.0*in)*cos(pos*90*deg),(-3.75*in-55.12/2.0*in+48.0*in)),ribsLogV,ribs_pv_name,h2oLogV,false,0,true);
      matPtr->SetReflector(h2oPhysV, g4pvs->GetVolume(ribs_pv_name), acrylicReflectivity, acrylicSigmaAlpha);
    }

    for (int pos=0; pos<10; pos++){
      G4RotationMatrix* Rotationfeet = new G4RotationMatrix();
      Rotationfeet->rotateX(0*deg);
      Rotationfeet->rotateY(0*deg);
      Rotationfeet->rotateZ(0*deg);
      std::string feet_pv_name = "feetPhysV_" + std::to_string(pos);
      new G4PVPlacement(0,G4ThreeVector(24.825/2.0*in*sin(pos*36*deg),24.825/2.0*in*cos(pos*36*deg),-((h2oHeight-0.5*pmtLegLength-0.5*pmtMinorAxis)/2.0)+4.13*in/2.0+tyvekThickness),feetLogV,feet_pv_name,h2oLogV,false,0,true);
      matPtr->SetReflector(h2oPhysV, g4pvs->GetVolume(feet_pv_name), acrylicReflectivity, acrylicSigmaAlpha);
    }

    new G4PVPlacement(0,G4ThreeVector(0,0,(+55.12*in-(h2oHeight-0.5*pmtLegLength-0.5*pmtMinorAxis))/2.0+4.13*in+tyvekThickness+(d2oHeight + 2.0*acrylicEndCapThickness)/2.0+17.40*in/2.0-1.0*in),strawLogV,"strawPhysV",h2oLogV,false,0,true);
    new G4PVPlacement(0,G4ThreeVector(0,0,0),airStrawLogV,"airStrawPhysV",strawLogV,false,0,true);

    for (int pos=0; pos<4; pos++){
      G4RotationMatrix* Rotationnarrow = new G4RotationMatrix();
      Rotationnarrow->rotateX(0*deg);
      Rotationnarrow->rotateY(0*deg);
      Rotationnarrow->rotateZ(45*deg+pos*90*deg);
      new G4PVPlacement(Rotationnarrow,G4ThreeVector((shieldLength+2.*in+2*0.06*in+airthickness)/2.0*sin(45*deg+pos*90*deg),(shieldLength+2.*in+2*0.06*in+airthickness)/2.0*cos(45*deg+pos*90*deg),0-sink),narrowAluminumLogV,"narrowAluminumPhysV",totalDetLogV,false,pos,true);
    }

    new G4PVPlacement(0,G4ThreeVector(0,-airthickness/2.0,0),narrowVetoLogV,"narrowVetoPhysV",narrowAluminumLogV,false,0,true);
    new G4PVPlacement(0,G4ThreeVector(0,+2.0*in/2.0,0),narrowAirLogV,"narrowAirPhysV",narrowAluminumLogV,false,0,true);

    for (int pos=0; pos<4; pos++){
      G4RotationMatrix* Rotationwide = new G4RotationMatrix();
      Rotationwide->rotateX(0*deg);
      Rotationwide->rotateY(0*deg);
      Rotationwide->rotateZ(pos*90*deg);
      if (pos==3){
        Rotationwide->rotateZ(180*deg);
      }
      new G4PVPlacement(Rotationwide,G4ThreeVector((shieldLength+2.*in+7*0.06*in+airthickness)/2.0*sin(pos*90*deg),(shieldLength+2.*in+7*0.06*in+airthickness)/2.0*cos(pos*90*deg),0-sink),wideAluminumLogV,"wideAluminumPhysV",totalDetLogV,false,pos,true);
    }

    new G4PVPlacement(0,G4ThreeVector(0,-airthickness/2.0,0),wideVetoLogV,"wideVetoPhysV",wideAluminumLogV,false,0,true);
    new G4PVPlacement(0,G4ThreeVector(0,+2.0*in/2.0,0),wideAirLogV,"wideAirPhysV",wideAluminumLogV,false,0,true);

    new G4PVPlacement(0,G4ThreeVector(0,0,+shieldHeight/2.0-88.67*in/2.0-0.38*in/2.0),shieldingLogV,"shieldingPhysV",totalDetLogV,false,0,true);

    return totalDetLogV;
}

// ============================================================
// Helper methods
// ============================================================

G4LogicalVolume *G4d2oCylindricalDetector::GetH2OModuleLogV(){
    G4Tubs* h2oModuleSolid = new G4Tubs("h2oModuleSolid",0.0,d2oLength/2.0,d2oHeight/2.0,0.0,360.0*deg);
    G4LogicalVolume *h2oModuleLogV = new G4LogicalVolume(h2oModuleSolid,
                                                         matPtr->GetMaterial( H2O ),
                                                         "h2oModuleLogV");
    G4VisAttributes *visAttH2OModule = new G4VisAttributes();
    visAttH2OModule->SetColour(G4Colour(0.0, 0.8, 0.8));
    h2oModuleLogV->SetVisAttributes(visAttH2OModule);
    return h2oModuleLogV;
}

G4LogicalVolume *G4d2oCylindricalDetector::GetTotalDetectorLogV(){
    G4Box* outerSolid = new G4Box("outerSolid",1.0*vetoOuterLength/2.0,1.0*vetoOuterLength/2.0,(86.43*in+2.*in+4*0.06*in+airthickness)/2.0);
    G4LogicalVolume *totalDetLogV = new G4LogicalVolume(outerSolid,
                                                        matPtr->GetMaterial( AIR ),
                                                        "totalDetLogV");
    std::cout << "\n vetoOuterHeight = \n\n";
    std::cout <<(vetoOuterHeight);
    G4VisAttributes *visTotal = new G4VisAttributes();
    visTotal->SetColour(G4Color::White());
    totalDetLogV->SetVisAttributes( visTotal );
    return totalDetLogV;
}

G4LogicalVolume *G4d2oCylindricalDetector::GetShieldingLogV() {
    G4LogicalVolume * shieldingLogV = 0;
    G4Tubs * shieldVessel = new G4Tubs("shieldVessel", 0.0, shieldLength/2.0, shieldHeight/2.0, 0.0, 360.0*deg);
    shieldingLogV = new G4LogicalVolume(shieldVessel,
                                        matPtr->GetMaterial( LEAD ),
                                        "shieldingLogV");
    G4VisAttributes *visShielding = new G4VisAttributes();
    visShielding->SetColour(G4Color::Brown());
    shieldingLogV->SetVisAttributes( visShielding );
    return shieldingLogV;
}

G4LogicalVolume *G4d2oCylindricalDetector::GetOuterVesselLogV(){
    G4Tubs* outerVessel = new G4Tubs("outerVessel",0.0,contOuterLength/2.0,contOuterHeight/2.0,0.0,360.0*deg);
    G4LogicalVolume *outerVesselLogV = new G4LogicalVolume(outerVessel,
                                                           matPtr->GetMaterial( STEEL ),
                                                           "outerVesselLogV");
    G4VisAttributes *visouterVessel = new G4VisAttributes();
    visouterVessel->SetColour(G4Colour(0.75, 0.75, 0.75));
    outerVesselLogV->SetVisAttributes(visouterVessel);
    return outerVesselLogV;
}

G4LogicalVolume *G4d2oCylindricalDetector::GetOuterAirselLogV(){
    G4Tubs* outerAirsel = new G4Tubs("outerAirsel",0.0,contOuterLength/2.0-0.5/2.0*in,contOuterHeight/2.0-0.5/2.0*in,0.0,360.0*deg);
    G4LogicalVolume *outerAirselLogV = new G4LogicalVolume(outerAirsel,
                                                           matPtr->GetMaterial( AIR ),
                                                           "outerAirselLogV");
    G4VisAttributes *visouterAirsel = new G4VisAttributes();
    visouterAirsel->SetColour(G4Colour(0.929, 0.91, 0.816));
    outerAirselLogV->SetVisAttributes(visouterAirsel);
    return outerAirselLogV;
}

G4LogicalVolume *G4d2oCylindricalDetector::GetInnerVetoLogV(){
    G4double halfheight = shieldHeight/2.0+1.0*muonVetoThickness;
    if (!input->GetBottomVeto() )
      halfheight = (shieldHeight + muonVetoThickness) / 2.0;
    G4Tubs* innerVeto = new G4Tubs("v",0.0,shieldLength/2.0+1.0*muonVetoThickness,
                                   halfheight,0.0,360.0*deg);
    G4LogicalVolume *innerVetoLogV = new G4LogicalVolume(innerVeto,
                                                           matPtr->GetMaterial( PLASTIC ),
                                                           "innerVetoLogV");
    G4VisAttributes *visinnerVeto = new G4VisAttributes();
    visinnerVeto->SetColour(G4Color::Cyan());
    innerVetoLogV->SetVisAttributes(visinnerVeto);
    return innerVetoLogV;
}

G4LogicalVolume *G4d2oCylindricalDetector::GetBasePlateLogV(){
    G4double basethickness = 3.556*in;
    G4Tubs* basePlate = new G4Tubs("basePlate",0.0,contOuterLength/2.0,(shieldHeight-contOuterHeight-1.75*in)/2.0,0.0,360*deg);
    G4LogicalVolume *basePlateLogV = new G4LogicalVolume(basePlate,
                                                         matPtr->GetMaterial( ALUMINUM ),
                                                         "basePlateLogV");
    G4VisAttributes *visbasePlate = new G4VisAttributes();
    visbasePlate->SetColour(G4Color::Gray());
    basePlateLogV->SetVisAttributes(visbasePlate);
    return basePlateLogV;
}

G4LogicalVolume *G4d2oCylindricalDetector::GetRibsLogV(){
    G4Box* ribs = new G4Box("ribs",1.0/2.0*in,3.375/2.0*in,2.0/2.0*in);
    G4LogicalVolume *ribsLogV = new G4LogicalVolume(ribs,
                                                    matPtr->GetMaterial( PMMA ),
                                                    "ribsLogV");
    G4VisAttributes *visribs = new G4VisAttributes();
    visribs->SetColour(G4Color::Magenta());
    ribsLogV->SetVisAttributes(visribs);
    return ribsLogV;
}

G4LogicalVolume *G4d2oCylindricalDetector::GetAirStrawLogV(){
    auto* airStraw = new G4Tubs("airStraw",0.0*in/2.0,0.5*in/2.0,17.4*in/2.0-1.0*in,0,360*deg);
    G4LogicalVolume *airStrawLogV = new G4LogicalVolume(airStraw,
                                                        matPtr->GetMaterial( AIR ),
                                                        "airStrawLogV");
    G4VisAttributes *visairStraw = new G4VisAttributes();
    visairStraw->SetColour(G4Color::White());
    airStrawLogV->SetVisAttributes(visairStraw);
    return airStrawLogV;
}

G4LogicalVolume *G4d2oCylindricalDetector::GetStrawLogV(){
    auto* straw = new G4Tubs("straw",0,1.0*in/2.0,17.4*in/2.0-1.0*in,0,360*deg);
    G4LogicalVolume *strawLogV = new G4LogicalVolume(straw,
                                                     matPtr->GetMaterial( PMMA ),
                                                     "strawLogV");
    G4VisAttributes *visstraw = new G4VisAttributes();
    visstraw->SetColour(G4Color::Magenta());
    strawLogV->SetVisAttributes(visstraw);
    return strawLogV;
}

G4LogicalVolume *G4d2oCylindricalDetector::GetLeadCapLogV(){
    auto* leadCap = new G4Tubs("leadCap",0,2.0*in/2.0,2.0*in/2.0,0,360*deg);
    G4LogicalVolume *leadCapLogV = new G4LogicalVolume(leadCap,
                                                       matPtr->GetMaterial( LEAD ),
                                                       "leadCapLogV");
    G4VisAttributes *visleadCap = new G4VisAttributes();
    visleadCap->SetColour(G4Color::Brown());
    leadCapLogV->SetVisAttributes(visleadCap);
    return leadCapLogV;
}

G4LogicalVolume *G4d2oCylindricalDetector::GetFeetLogV(){
    G4Box* feet = new G4Box("feet",1.38/2.0*in,1.38/2.0*in,4.13/2.0*in);
    G4LogicalVolume *feetLogV = new G4LogicalVolume(feet,
                                                    matPtr->GetMaterial( PMMA ),
                                                    "feetLogV");
    G4VisAttributes *visfeet = new G4VisAttributes();
    visfeet->SetColour(G4Color::Magenta());
    feetLogV->SetVisAttributes(visfeet);
    return feetLogV;
}

G4LogicalVolume *G4d2oCylindricalDetector::GetTopAluminumLogV(){
    G4double aluminumthickness = 0.06;
    G4Box* topAluminum = new G4Box("v",(44.688+2.0*aluminumthickness)/2.0*in,(44.688+2.0*aluminumthickness)/2.0*in,(2.000+2.0*aluminumthickness)/2.0*in+airthickness/2.0);
    G4LogicalVolume *topAluminumLogV = new G4LogicalVolume(topAluminum,
                                                           matPtr->GetMaterial( ALUMINUM ),
                                                           "topAluminumLogV");
    G4VisAttributes *vistopAluminum = new G4VisAttributes();
    vistopAluminum->SetColour(G4Color::Gray());
    topAluminumLogV->SetVisAttributes(vistopAluminum);
    return topAluminumLogV;
}

G4LogicalVolume *G4d2oCylindricalDetector::GetTopAirLogV(){
    G4double aluminumthickness = 0.06;
    G4Box* topAir = new G4Box("v",(44.688)/2.0*in,(44.688)/2.0*in,airthickness/2.0);
    G4LogicalVolume *topAirLogV = new G4LogicalVolume(topAir,
                                                      matPtr->GetMaterial( AIR ),
                                                      "topAirLogV");
    G4VisAttributes *vistopAir = new G4VisAttributes();
    vistopAir->SetColour(G4Color::White());
    topAirLogV->SetVisAttributes(vistopAir);
    return topAirLogV;
}

G4LogicalVolume *G4d2oCylindricalDetector::GetNarrowAirLogV(){
    G4double aluminumthickness = 0.06;
    G4Box* narrowAir = new G4Box("v",11.688/2.0*in,airthickness/2.0,86.430/2.0*in);
    G4LogicalVolume *narrowAirLogV = new G4LogicalVolume(narrowAir,
                                                         matPtr->GetMaterial( AIR ),
                                                         "narrowAirLogV");
    G4VisAttributes *visnarrowAir = new G4VisAttributes();
    visnarrowAir->SetColour(G4Color::White());
    narrowAirLogV->SetVisAttributes(visnarrowAir);
    return narrowAirLogV;
}

G4LogicalVolume *G4d2oCylindricalDetector::GetWideAirLogV(){
    G4double aluminumthickness = 0.06;
    G4Box* wideAir = new G4Box("v",25.188/2.0*in,airthickness/2.0,86.430/2.0*in);
    G4LogicalVolume *wideAirLogV = new G4LogicalVolume(wideAir,
                                                       matPtr->GetMaterial( AIR ),
                                                       "wideAirLogV");
    G4VisAttributes *viswideAir = new G4VisAttributes();
    viswideAir->SetColour(G4Color::White());
    wideAirLogV->SetVisAttributes(viswideAir);
    return wideAirLogV;
}

G4LogicalVolume *G4d2oCylindricalDetector::GetNarrowAluminumLogV(){
    G4double aluminumthickness = 0.06;
    G4Box* narrowAluminum = new G4Box("v",(11.688+2.0*aluminumthickness)/2.0*in,(2.000+2.0*aluminumthickness)/2.0*in+airthickness/2.0,(86.430+2.0*aluminumthickness)/2.0*in);
    G4LogicalVolume *narrowAluminumLogV = new G4LogicalVolume(narrowAluminum,
                                                              matPtr->GetMaterial( ALUMINUM ),
                                                              "narrowAluminiumLogV");
    G4VisAttributes *visnarrowAluminum = new G4VisAttributes();
    visnarrowAluminum->SetColour(G4Color::Gray());
    narrowAluminumLogV->SetVisAttributes(visnarrowAluminum);
    return narrowAluminumLogV;
}

G4LogicalVolume *G4d2oCylindricalDetector::GetWideAluminumLogV(){
    G4double aluminumthickness = 0.06;
    G4Box* wideAluminum = new G4Box("v",(25.188+2.0*aluminumthickness)/2.0*in,(2.000+2.0*aluminumthickness)/2.0*in+airthickness/2.0,(86.430+2.0*aluminumthickness)/2.0*in);
    G4LogicalVolume *wideAluminumLogV = new G4LogicalVolume(wideAluminum,
                                                            matPtr->GetMaterial( ALUMINUM ),
                                                            "wideAluminiumLogV");
    G4VisAttributes *viswideAluminum = new G4VisAttributes();
    viswideAluminum->SetColour(G4Color::Gray());
    wideAluminumLogV->SetVisAttributes(viswideAluminum);
    return wideAluminumLogV;
}

G4LogicalVolume *G4d2oCylindricalDetector::GetTopVetoLogV(){
    G4double halfheight = shieldHeight/2.0+2.0*muonVetoThickness;
    if (!input->GetBottomVeto() )
      halfheight = (shieldHeight + 2.0*muonVetoThickness)/2.0;
    G4Box* topVeto = new G4Box("v",44.688/2.0*in,44.688/2.0*in,2.000/2.0*in);
    G4LogicalVolume *topVetoLogV = new G4LogicalVolume(topVeto,
                                                       matPtr->GetMaterial( VINYLTOLUENE ),
                                                       "topVetoLogV");
    G4VisAttributes *vistopVeto = new G4VisAttributes();
    vistopVeto->SetColour(G4Color::Yellow());
    topVetoLogV->SetVisAttributes(vistopVeto);
    auto* sdManager = G4SDManager::GetSDMpointer();
    auto* vetoSd = new G4d2oSensitiveDetector("topVetoSD", topVetoLogV->GetMaterial()->GetName(), 0);
    sdManager->AddNewDetector(vetoSd);
    topVetoLogV->SetSensitiveDetector(vetoSd);
    return topVetoLogV;
}

G4LogicalVolume *G4d2oCylindricalDetector::GetNarrowVetoLogV(){
    G4Box* narrowVeto = new G4Box("v",11.688/2.0*in,2.000/2.0*in,86.430/2.0*in);
    G4LogicalVolume *narrowVetoLogV = new G4LogicalVolume(narrowVeto,
                                                          matPtr->GetMaterial( VINYLTOLUENE ),
                                                          "narrowVetoLogV");
    G4VisAttributes *visnarrowVeto = new G4VisAttributes();
    visnarrowVeto->SetColour(G4Color::Yellow());
    narrowVetoLogV->SetVisAttributes(visnarrowVeto);
    auto* sdManager = G4SDManager::GetSDMpointer();
    auto* vetoNo1Sd = new G4d2oSensitiveDetector("VetoNu1SD", narrowVetoLogV->GetMaterial()->GetName(), 0);
    sdManager->AddNewDetector(vetoNo1Sd);
    narrowVetoLogV->SetSensitiveDetector(vetoNo1Sd);
    return narrowVetoLogV;
}

G4LogicalVolume *G4d2oCylindricalDetector::GetWideVetoLogV(){
    G4Box* wideVeto = new G4Box("v",25.188/2.0*in,2.000/2.0*in,86.430/2.0*in);
    G4LogicalVolume *wideVetoLogV = new G4LogicalVolume(wideVeto,
                                                        matPtr->GetMaterial( VINYLTOLUENE ),
                                                        "wideVetoLogV");
    G4VisAttributes *viswideVeto = new G4VisAttributes();
    viswideVeto->SetColour(G4Color::Yellow());
    wideVetoLogV->SetVisAttributes(viswideVeto);
    auto* sdManager = G4SDManager::GetSDMpointer();
    auto* vetoNo2Sd = new G4d2oSensitiveDetector("VetoNu2SD", wideVetoLogV->GetMaterial()->GetName(), 0);
    sdManager->AddNewDetector(vetoNo2Sd);
    wideVetoLogV->SetSensitiveDetector(vetoNo2Sd);
    return wideVetoLogV;
}

G4LogicalVolume *G4d2oCylindricalDetector::GetH2OLogV(){
    G4Tubs* h2oSolid = new G4Tubs("h2oSolid",0.0,(h2oLength)/2.0,(h2oHeight-0.5*pmtLegLength-0.5*pmtMinorAxis)/2.0,0.0,360.0*deg);
    G4LogicalVolume *h2oLogV = new G4LogicalVolume(h2oSolid,
                                                   matPtr->GetMaterial( H2O ),
                                                   "h2oLogV");
    G4VisAttributes *visLightWater = new G4VisAttributes();
    visLightWater->SetColour(G4Color::Cyan());
    h2oLogV->SetVisAttributes( visLightWater );
    return h2oLogV;
}

G4LogicalVolume *G4d2oCylindricalDetector::GetAcrylicLogV(){
    G4double acrylicLength = d2oLength + 2.0*acrylicThickness;
    G4double acrylicHeight = d2oHeight + 2.0*acrylicEndCapThickness;
    G4Tubs* acrylicSolid = new G4Tubs("acrylicSolid",0.0,acrylicLength/2.0,acrylicHeight/2.0,0.0,360.0*deg);
    G4LogicalVolume *acrylicLogV = new G4LogicalVolume(acrylicSolid,
                                                       matPtr->GetMaterial( PMMA ),
                                                       "acrylicLogV");
    G4VisAttributes *visAcrylic = new G4VisAttributes();
    visAcrylic->SetColour(G4Color::Magenta());
    acrylicLogV->SetVisAttributes( visAcrylic );
    return acrylicLogV;
}

G4LogicalVolume *G4d2oCylindricalDetector::GetTyvekLiningCapLogV(const char * name, bool boundPMTs=true){
    G4LogicalVolume *capLogV = 0;
    G4double width = (h2oLength-tyvekThickness)/2.0;
    G4double height = (tyvekThickness)/2.0;
    G4VSolid * capSolid = new G4Tubs(Form("%sSolid", name), 0.0, width, height, 0.0, 360.0*deg);
    if (boundPMTs) {
        G4Tubs *capSub = new G4Tubs("subCap", 0.0, pmtMajorAxis/2.0 + 1.0*mm, height+1.0*cm, 0.0 ,360.0);
        G4Transform3D transform;
        int counter = 0;
        for (int bin=1; bin<=hPanel[2]->GetNbinsX(); ++bin) {
            if (hPanel[2]->GetBinContent(bin) <= 0.0) continue;
            transform = G4Transform3D(G4RotationMatrix(), G4ThreeVector(hPanel[0]->GetBinContent(bin), hPanel[1]->GetBinContent(bin), -1.0*mm));
            capSolid = new G4SubtractionSolid(Form("%sSolid", name), capSolid, capSub, transform);
            ++counter;
        }
        capLogV = new G4LogicalVolume(capSolid,
                                      matPtr->GetMaterial( TYVEK ),
                                      Form("%sLogV", name) );
        G4VisAttributes *visCap = new G4VisAttributes();
        visCap->SetColour(G4Color::Red());
        capLogV->SetVisAttributes( visCap );
    } else {
        capLogV = new G4LogicalVolume(capSolid,
                                      matPtr->GetMaterial( TYVEK ),
                                      Form("%sLogV", name) );
        G4VisAttributes *visCap = new G4VisAttributes();
        visCap->SetColour(G4Color::Gray());
        capLogV->SetVisAttributes( visCap );
    }
    return capLogV;
}

G4LogicalVolume *G4d2oCylindricalDetector::GetTyvekLiningLogV(){
    G4double tyvekheight = (d2oHeight/2.0) + acrylicEndCapThickness + h2oTCthickness + (pmtMinorAxis/2.0) - (pmtLegLength+3.0*tyvekThickness)/2.0 + pmtLegLength;
    G4Tubs * tyvekSolid = new G4Tubs("tyvekSolid",(h2oLength-tyvekThickness)/2.0, h2oLength/2.0, contOuterHeight/2.0-outerContainerThickness, 0.0, 360.0*deg);
    G4LogicalVolume *theSideLiningLogV = new G4LogicalVolume(tyvekSolid,
                                                             matPtr->GetMaterial( TYVEK ),
                                                             "theSideLiningLogV");
    G4VisAttributes *visTyvekLining = new G4VisAttributes();
    visTyvekLining->SetColour(G4Color::Gray());
    theSideLiningLogV->SetVisAttributes( visTyvekLining );
    return theSideLiningLogV;
}

void G4d2oCylindricalDetector::PlacePMTs(G4LogicalVolume *thePMTLogV, G4LogicalVolume *theMotherLogV){
    G4int iCopy = 0;
    //histograms that hold coordinates of the front faces of the PMTs
    //x-y plane, +z
    G4RotationMatrix *rotTop = new G4RotationMatrix();
    rotTop->rotateY(180.0*deg);

    if (!(numInY%2) && !(numInX%2)) {
      G4double yshift = 0.0;
      if (!(numInY%2)) {
        yshift = (pmtMajorAxis / 2.0) / TMath::Sqrt(3.0);
      }
      for (int iY=0; iY<numInY; ++iY) {
        G4int ycount = iY - (numInY/2);
        G4double ycoord = static_cast<G4double>(ycount);
        ycoord *= TMath::Sqrt(3.0) * ( (pmtMajorAxis/2.0) + minGapBetweenPMTs);
        ycoord += yshift;
        G4int nx = numInX - TMath::Abs(ycount);
        G4double xshift = 0.0;
        if (!(nx%2)) {
          xshift = pmtMajorAxis/2.0;
        }
        for (int iX=0; iX<nx; ++iX) {
          G4int xcount = iX - (nx/2);
          G4double xcoord = static_cast<G4double>(xcount);
          xcoord *= pmtMajorAxis + minGapBetweenPMTs;
          xcoord += xshift;
          G4double radius = TMath::Sqrt(TMath::Power(ycoord, 2.0) + TMath::Power(xcoord, 2.0) );
          radius += (pmtMajorAxis / 2.0);
          if (radius >= (h2oLength-tyvekThickness)/2.0)
            continue;
          G4ThreeVector thisPMTPlacement;
          if (pmtDiameter > 0.0) {
            thisPMTPlacement = G4ThreeVector(xcoord,ycoord,h2oHeight/2.0-tyvekThickness+pmtLegLength);
          } else {
            thisPMTPlacement = G4ThreeVector(xcoord,ycoord,(h2oHeight-tyvekThickness)/2.0 - pmtMinorAxis/2.0 - pmtLegLength+pmtLegLength - 2*50.0*mm);
          }
          new G4PVPlacement(rotTop,thisPMTPlacement,thePMTLogV,Form("pmtPhysV_%d",iCopy),theMotherLogV,true,iCopy,true);
          ++iCopy;
        }
      }
    } else {
      G4double yshift = 0.0;
      if (!(numInY%2)) {
        yshift = pmtMajorAxis / 2.0;
      }
      for (int iY=0; iY<numInY; ++iY) {
        G4int ycount = iY - (numInY/2);
        G4double ycoord = static_cast<G4double>(ycount);
        ycoord *= TMath::Sqrt(3.0) * ( (pmtMajorAxis/2.0) + minGapBetweenPMTs);
        ycoord += yshift;
        G4int nx = numInX - TMath::Abs(ycount);
        G4double xshift = 0.0;
        if (!(nx%2)) {
          xshift = pmtMajorAxis/2.0;
        }
        for (int iX=0; iX<nx; ++iX) {
          G4int xcount = iX - (nx/2);
          G4double xcoord = static_cast<G4double>(xcount);
          xcoord *= pmtMajorAxis + minGapBetweenPMTs;
          xcoord += xshift;
          G4double radius = TMath::Sqrt(TMath::Power(ycoord, 2.0) + TMath::Power(xcoord, 2.0) );
          radius += (pmtMajorAxis / 2.0);
          if (radius >= (h2oLength+tyvekThickness)/2.0)
            continue;
          G4ThreeVector thisPMTPlacement;
          if (pmtDiameter > 0.0) {
            thisPMTPlacement = G4ThreeVector(xcoord,ycoord,h2oHeight/2.0-tyvekThickness+pmtLegLength);
          } else {
            thisPMTPlacement = G4ThreeVector(xcoord,ycoord,h2oHeight/2.0-tyvekThickness - pmtMinorAxis/2.0 - pmtLegLength+pmtLegLength);
          }
          new G4PVPlacement(rotTop,thisPMTPlacement,thePMTLogV,Form("pmtPhysV_%d",iCopy),theMotherLogV,true,iCopy,true);
          ++iCopy;
        }
      }
    }

    //x-y plane, -z
    if(iUseBottomPMTs) {
      G4RotationMatrix *rotBottom = new G4RotationMatrix();
      if (!(numInY%2) && !(numInX%2)) {
        G4double yshift = 0.0;
        if (!(numInY%2)) {
          yshift = (pmtMajorAxis / 2.0) / TMath::Sqrt(3.0);
        }
        for (int iY=0; iY<numInY; ++iY) {
          G4int ycount = iY - (numInY/2);
          G4double ycoord = static_cast<G4double>(ycount);
          ycoord *= TMath::Sqrt(3.0) * ( (pmtMajorAxis/2.0) + minGapBetweenPMTs);
          ycoord += yshift;
          G4int nx = numInX - TMath::Abs(ycount);
          G4double xshift = 0.0;
          if (!(nx%2)) {
            xshift = pmtMajorAxis/2.0;
          }
          for (int iX=0; iX<nx; ++iX) {
            G4int xcount = iX - (nx/2);
            G4double xcoord = static_cast<G4double>(xcount);
            xcoord *= pmtMajorAxis + minGapBetweenPMTs;
            xcoord += xshift;
            G4double radius = TMath::Sqrt(TMath::Power(ycoord, 2.0) + TMath::Power(xcoord, 2.0) );
            radius += (pmtMajorAxis / 2.0);
            if (radius >= h2oLength/2.0)
              continue;
            G4ThreeVector thisPMTPlacement;
            if (pmtDiameter > 0.0) {
              thisPMTPlacement = G4ThreeVector(xcoord,ycoord,-(h2oHeight/2.0-tyvekThickness));
            } else {
              thisPMTPlacement = G4ThreeVector(xcoord,ycoord,-(h2oHeight/2.0-tyvekThickness - pmtMinorAxis/2.0 - pmtLegLength));
            }
            new G4PVPlacement(rotBottom,thisPMTPlacement,thePMTLogV,Form("pmtPhysV_%d",iCopy),theMotherLogV,true,iCopy,true);
            ++iCopy;
          }
        }
      } else {
        G4double yshift = 0.0;
        if (!(numInY%2)) {
          yshift = pmtMajorAxis / 2.0;
        }
        for (int iY=0; iY<numInY; ++iY) {
          G4int ycount = iY - (numInY/2);
          G4double ycoord = static_cast<G4double>(ycount);
          ycoord *= TMath::Sqrt(3.0) * ( (pmtMajorAxis/2.0) + minGapBetweenPMTs);
          ycoord += yshift;
          G4int nx = numInX - TMath::Abs(ycount);
          G4double xshift = 0.0;
          if (!(nx%2)) {
            xshift = pmtMajorAxis/2.0;
          }
          for (int iX=0; iX<nx; ++iX) {
            G4int xcount = iX - (nx/2);
            G4double xcoord = static_cast<G4double>(xcount);
            xcoord *= pmtMajorAxis + minGapBetweenPMTs;
            xcoord += xshift;
            G4double radius = TMath::Sqrt(TMath::Power(ycoord, 2.0) + TMath::Power(xcoord, 2.0) );
            radius += (pmtMajorAxis / 2.0);
            if (radius >= h2oLength/2.0)
              continue;
            G4ThreeVector thisPMTPlacement;
            if (pmtDiameter > 0.0) {
              thisPMTPlacement = G4ThreeVector(xcoord,ycoord,-(h2oHeight/2.0-tyvekThickness));
            } else {
              thisPMTPlacement = G4ThreeVector(xcoord,ycoord,-(h2oHeight/2.0-tyvekThickness - pmtMinorAxis/2.0 - pmtLegLength));
            }
            new G4PVPlacement(rotBottom,thisPMTPlacement,thePMTLogV,Form("pmtPhysV_%d",iCopy),theMotherLogV,true,iCopy,true);
            ++iCopy;
          }
        }
      }
    }

    totPMT = iCopy;
    numInRows = totPMT;
    numRows = totPMT;

    // Now store PMT locations...
    iCopy = 1;
    for(G4int iDim=0; iDim<3; iDim++)
        hPanel[iDim] = new TH1D(Form("hPanel[%d]",iDim),Form("Dimension %d",iDim),totPMT,-0.5,totPMT-0.5);

    if (!(numInY%2) && !(numInX%2)) {
      G4double yshift = 0.0;
      if (!(numInY%2)) {
        yshift = (pmtMajorAxis / 2.0) / TMath::Sqrt(3.0);
      }
      for (int iY=0; iY<numInY; ++iY) {
        G4int ycount = iY - (numInY/2);
        G4double ycoord = static_cast<G4double>(ycount);
        ycoord *= TMath::Sqrt(3.0) * ( (pmtMajorAxis/2.0) + minGapBetweenPMTs);
        ycoord += yshift;
        G4int nx = numInX - TMath::Abs(ycount);
        G4double xshift = 0.0;
        if (!(nx%2)) {
          xshift = pmtMajorAxis/2.0;
        }
        for (int iX=0; iX<nx; ++iX) {
          G4int xcount = iX - (nx/2);
          G4double xcoord = static_cast<G4double>(xcount);
          xcoord *= pmtMajorAxis + minGapBetweenPMTs;
          xcoord += xshift;
          G4double radius = TMath::Sqrt(TMath::Power(ycoord, 2.0) + TMath::Power(xcoord, 2.0) );
          radius += (pmtMajorAxis / 2.0);
          if (radius >= h2oLength/2.0)
            continue;
          G4ThreeVector thisPMTPlacement;
          if (pmtDiameter > 0.0) {
            thisPMTPlacement = G4ThreeVector(xcoord,ycoord,h2oHeight/2.0-tyvekThickness+pmtLegLength);
          } else {
            thisPMTPlacement = G4ThreeVector(xcoord,ycoord,h2oHeight/2.0-tyvekThickness - pmtMinorAxis/2.0 - pmtLegLength+pmtLegLength);
          }
          hPanel[0]->SetBinContent(iCopy,thisPMTPlacement.x());
          hPanel[1]->SetBinContent(iCopy,thisPMTPlacement.y());
          if (pmtDiameter > 0.0) {
            hPanel[2]->SetBinContent(iCopy,thisPMTPlacement.z()-pmtMinorAxis/2.0);
          } else {
            hPanel[2]->SetBinContent(iCopy,thisPMTPlacement.z()-pmtMinorAxis/2.0 - pmtLegLength);
          }
          ++iCopy;
        }
      }
    } else {
      G4double yshift = 0.0;
      if (!(numInY%2)) {
        yshift = pmtMajorAxis / 2.0;
      }
      for (int iY=0; iY<numInY; ++iY) {
        G4int ycount = iY - (numInY/2);
        G4double ycoord = static_cast<G4double>(ycount);
        ycoord *= TMath::Sqrt(3.0) * ( (pmtMajorAxis/2.0) + minGapBetweenPMTs);
        ycoord += yshift;
        G4int nx = numInX - TMath::Abs(ycount);
        G4double xshift = 0.0;
        if (!(nx%2)) {
          xshift = pmtMajorAxis/2.0;
        }
        for (int iX=0; iX<nx; ++iX) {
          G4int xcount = iX - (nx/2);
          G4double xcoord = static_cast<G4double>(xcount);
          xcoord *= pmtMajorAxis + minGapBetweenPMTs;
          xcoord += xshift;
          G4double radius = TMath::Sqrt(TMath::Power(ycoord, 2.0) + TMath::Power(xcoord, 2.0) );
          radius += (pmtMajorAxis / 2.0);
          if (radius >= h2oLength/2.0)
            continue;
          G4ThreeVector thisPMTPlacement;
          if (pmtDiameter > 0.0) {
            thisPMTPlacement = G4ThreeVector(xcoord,ycoord,h2oHeight/2.0-tyvekThickness+pmtLegLength);
          } else {
            thisPMTPlacement = G4ThreeVector(xcoord,ycoord,h2oHeight/2.0-tyvekThickness - pmtMinorAxis/2.0 - pmtLegLength+pmtLegLength);
          }
          hPanel[0]->SetBinContent(iCopy,thisPMTPlacement.x());
          hPanel[1]->SetBinContent(iCopy,thisPMTPlacement.y());
          if (pmtDiameter > 0.0) {
            hPanel[2]->SetBinContent(iCopy,thisPMTPlacement.z()-pmtMinorAxis/2.0);
          } else {
            hPanel[2]->SetBinContent(iCopy,thisPMTPlacement.z()-pmtMinorAxis/2.0 - pmtLegLength);
          }
          ++iCopy;
        }
      }
    }

    //x-y plane, -z
    if(iUseBottomPMTs) {
      G4RotationMatrix *rotBottom = new G4RotationMatrix();
      if (!(numInY%2) && !(numInX%2)) {
        G4double yshift = 0.0;
        if (!(numInY%2)) {
          yshift = (pmtMajorAxis / 2.0) / TMath::Sqrt(3.0);
        }
        for (int iY=0; iY<numInY; ++iY) {
          G4int ycount = iY - (numInY/2);
          G4double ycoord = static_cast<G4double>(ycount);
          ycoord *= TMath::Sqrt(3.0) * ( (pmtMajorAxis/2.0) + minGapBetweenPMTs);
          ycoord += yshift;
          G4int nx = numInX - TMath::Abs(ycount);
          G4double xshift = 0.0;
          if (!(nx%2)) {
            xshift = pmtMajorAxis/2.0;
          }
          for (int iX=0; iX<nx; ++iX) {
            G4int xcount = iX - (nx/2);
            G4double xcoord = static_cast<G4double>(xcount);
            xcoord *= pmtMajorAxis + minGapBetweenPMTs;
            xcoord += xshift;
            G4double radius = TMath::Sqrt(TMath::Power(ycoord, 2.0) + TMath::Power(xcoord, 2.0) );
            radius += (pmtMajorAxis / 2.0);
            if (radius >= h2oLength/2.0)
              continue;
            G4ThreeVector thisPMTPlacement;
            if (pmtDiameter > 0.0) {
              thisPMTPlacement = G4ThreeVector(xcoord,ycoord,-(h2oHeight/2.0-tyvekThickness));
            } else {
              thisPMTPlacement = G4ThreeVector(xcoord,ycoord,-(h2oHeight/2.0-tyvekThickness - pmtMinorAxis/2.0 - pmtLegLength));
            }
            hPanel[0]->SetBinContent(iCopy,thisPMTPlacement.x());
            hPanel[1]->SetBinContent(iCopy,thisPMTPlacement.y());
            if (pmtDiameter > 0.0) {
              hPanel[2]->SetBinContent(iCopy,thisPMTPlacement.z()-pmtMinorAxis/2.0);
            } else {
              hPanel[2]->SetBinContent(iCopy,thisPMTPlacement.z()-pmtMinorAxis/2.0-pmtLegLength);
            }
            ++iCopy;
          }
        }
      } else {
        G4double yshift = 0.0;
        if (!(numInY%2)) {
          yshift = pmtMajorAxis / 2.0;
        }
        for (int iY=0; iY<numInY; ++iY) {
          G4int ycount = iY - (numInY/2);
          G4double ycoord = static_cast<G4double>(ycount);
          ycoord *= TMath::Sqrt(3.0) * ( (pmtMajorAxis/2.0) + minGapBetweenPMTs);
          ycoord += yshift;
          G4int nx = numInX - TMath::Abs(ycount);
          G4double xshift = 0.0;
          if (!(nx%2)) {
            xshift = pmtMajorAxis/2.0;
          }
          for (int iX=0; iX<nx; ++iX) {
            G4int xcount = iX - (nx/2);
            G4double xcoord = static_cast<G4double>(xcount);
            xcoord *= pmtMajorAxis + minGapBetweenPMTs;
            xcoord += xshift;
            G4double radius = TMath::Sqrt(TMath::Power(ycoord, 2.0) + TMath::Power(xcoord, 2.0) );
            radius += (pmtMajorAxis / 2.0);
            if (radius >= h2oLength/2.0)
              continue;
            G4ThreeVector thisPMTPlacement;
            if (pmtDiameter > 0.0) {
              thisPMTPlacement = G4ThreeVector(xcoord,ycoord,-(h2oHeight/2.0-tyvekThickness));
            } else {
              thisPMTPlacement = G4ThreeVector(xcoord,ycoord,-(h2oHeight/2.0-tyvekThickness - pmtMinorAxis/2.0 - pmtLegLength));
            }
            hPanel[0]->SetBinContent(iCopy,thisPMTPlacement.x());
            hPanel[1]->SetBinContent(iCopy,thisPMTPlacement.y());
            if (pmtDiameter > 0.0) {
              hPanel[2]->SetBinContent(iCopy,thisPMTPlacement.z()-pmtMinorAxis/2.0);
            } else {
              hPanel[2]->SetBinContent(iCopy,thisPMTPlacement.z()-pmtMinorAxis/2.0-pmtLegLength);
            }
            ++iCopy;
          }
        }
      }
    }

    std::cout << "\n contOuterHeight = \n\n";
    std::cout <<(contOuterHeight);
    std::cout << "\n h2oHeight = \n\n";
    std::cout <<(h2oHeight);
    std::cout << "\n shieldHeight = \n\n";
    std::cout <<(shieldHeight);
    std::cout << "\n vetoOuterLength = \n\n";
    std::cout <<(vetoOuterLength);
    std::cout << "\n offset = \n\n";
}

