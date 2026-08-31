#ifndef G4d2oCylindricalDetector_H
#define G4d2oCylindricalDetector_H 1

#include "G4d2oDetector.hh"

class G4d2oCylindricalDetector : public G4d2oDetector
{
public:
    
    G4d2oCylindricalDetector();
    ~G4d2oCylindricalDetector();
    
    virtual G4LogicalVolume * GetDetector();
    
protected:
    virtual void Initialize();
    virtual void DetermineSpacing();
    
    G4LogicalVolume *GetTotalDetectorLogV();
    G4LogicalVolume *GetShieldingLogV();
    G4LogicalVolume *GetOuterVesselLogV();
    G4LogicalVolume *GetOuterAirselLogV();
    G4LogicalVolume *GetTyvekLiningLogV();
    G4LogicalVolume *GetTyvekLiningCapLogV(const char * name, bool boundPMTs);
    G4LogicalVolume *GetH2OLogV();
    G4LogicalVolume *GetAcrylicLogV();
    G4LogicalVolume *GetD2OLogV();
    G4LogicalVolume *GetbubbleLogV();
    G4LogicalVolume *GetTopVetoLogV();
    G4LogicalVolume *GetInnerVetoLogV();
    G4LogicalVolume *GetNarrowVetoLogV();
    G4LogicalVolume *GetWideVetoLogV();
    G4LogicalVolume *GetTopAluminumLogV();
    G4LogicalVolume *GetBasePlateLogV();
    G4LogicalVolume *GetRibsLogV();
    G4LogicalVolume *GetAirStrawLogV();
    G4LogicalVolume *GetLeadCapLogV();
    G4LogicalVolume *GetStrawLogV();  
    G4LogicalVolume *GetFeetLogV();
    G4LogicalVolume *GetTopAirLogV();
    G4LogicalVolume *GetNarrowAluminumLogV();
    G4LogicalVolume *GetWideAluminumLogV();
    G4LogicalVolume *GetWideAirLogV();
    G4LogicalVolume *GetNarrowAirLogV();

    // ============================================================
    // NEW: Module 2 target volume (both H2O)
    // ============================================================
    G4LogicalVolume *GetH2OModuleLogV();
  
    //-----------------------------------------------------------
    
    G4LogicalVolume *GetTyvekLiningLLogV();
    G4LogicalVolume *GetTyvekLiningFLogV();
    G4LogicalVolume *GetTyvekLiningTLogV();
    
    //-----------------------------------------------------------

    void PlacePMTs(G4LogicalVolume *thePMTLogV, G4LogicalVolume *theMotherLogV);

    G4int iUseBottomPMTs;
    
    G4double acrylicEndCapThickness;
    
    //-----------------------------------------------------------
    
    G4double h2oInnerRadius;
    G4double TCRadius; 
    
    G4double tyvekHeight; 
    G4double tyvekLength; 
    G4double tyvekOuterR; 
    G4double tyvekInnerR; 
    G4int numInCircle;
    G4double spacingInCircle;
    G4double totalAngle;
    G4double bottomThickness;
  
    //-----------------------------------------------------------

    
};//END of class G4d2oCylindricalDetector

#endif
