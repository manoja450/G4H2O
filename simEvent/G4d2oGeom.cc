#include "G4d2oGeom.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TDirectoryFile.h"
#include "TH1.h"
#include <iostream>
#include "TFile.h"
#include "TH2.h"
#include "G4SystemOfUnits.hh"

ClassImp(G4d2oGeom)

G4d2oGeom* G4d2oGeom::_geom = 0;

G4d2oGeom::G4d2oGeom()
{
  TString xsnDFile(gSystem->Getenv("G4D2OXSND"));
  TString xsnOFile(gSystem->Getenv("G4D2OXSNO"));
  
  if(gDirectory&&gDirectory->Get("pmtPositions"))
  {
    auto pmtDir = (TDirectoryFile*)gDirectory->Get("pmtPositions");
    TH1 *pmtPosX = 0, *pmtPosY = 0, *pmtPosZ = 0;
    if(pmtDir&&pmtDir->Get("pmtPosX")){
      pmtPosX = (TH1*)pmtDir->Get("pmtPosX");
    }
    if(pmtDir&&pmtDir->Get("pmtPosY")){
      pmtPosY = (TH1*)pmtDir->Get("pmtPosY");
    }
    if(pmtDir&&pmtDir->Get("pmtPosZ")){
      pmtPosZ = (TH1*)pmtDir->Get("pmtPosZ");
    }
    
    if(pmtPosX&&pmtPosY&&pmtPosZ){
      _detpos.resize(pmtPosX->GetNbinsX());
      for(int ibin = 1; ibin <= pmtPosX->GetNbinsX(); ibin++){
        _detpos[ibin-1].SetXYZ(pmtPosX->GetBinContent(ibin),
                               pmtPosY->GetBinContent(ibin),
                               pmtPosZ->GetBinContent(ibin));
      }
    }else{
      std::cerr<<"Warning unable to initialize pmt positions"<<std::endl;
    }
    
  }
  
  if(xsnDFile!=""){
    TDirectory* prevDir = 0;
    if(gDirectory) prevDir = gDirectory;
    auto fin = TFile::Open(xsnDFile.Data());
    if(fin&&fin->Get("fluxW")){
      _hxscnD = (TH2*)fin->Get("fluxW");
      _hxscnD->SetDirectory(0);
      _hxscnD->Scale(1.0/_hxscnD->Integral("width"));
    }else{
      _hxscnD = 0;
    }
    delete fin;
    prevDir->cd();
  }
  if(xsnOFile!=""){
    TDirectory* prevDir = 0;
    if(gDirectory) prevDir = gDirectory;
    auto fin = TFile::Open(xsnOFile.Data());
    if(fin&&fin->Get("haxtondd2")){
      _hxscnO = (TH2*)fin->Get("haxtondd2");
      _hxscnO->SetDirectory(0);
      _hxscnO->Scale(1.0/_hxscnO->Integral("width"));
    }else{
      _hxscnO = 0;
    }
    delete fin;
    prevDir->cd();
  }
  
//  Bottom left: (x, y, z), (17.845, -8.076, 4.320)
//  Bottom right: (17.845, -8.066, 3.291)
//  Top right: (17.847, -6.402, 3.282)
//  Bottom right against wall*: (16.943, -8.074, 4.321)
// * not sure about this one
// In this simulation, the x and y coordinates are swapped.
  _neutrinoSourcePos.SetXYZ(-2.0/m,-(17.845 - 0.5)/m , 8.066/m+(140/2.0+20.0+2.54*(1+10))/cm);
}


double G4d2oGeom::ProbO(const TVector3 epos, const TVector3 edir, const double eEnergy) const{
  if(!_hxscnO) return 0.0;
  if(eEnergy<_hxscnO->GetXaxis()->GetXmin() || eEnergy>_hxscnO->GetXaxis()->GetXmax())
    return 0.0;
  TVector3 nuDir = (epos-_neutrinoSourcePos);
  double costhe = nuDir.Dot(edir)/nuDir.Mag()/edir.Mag();
  return _hxscnO->Interpolate(eEnergy, costhe);
}

double G4d2oGeom::ProbD(const TVector3 epos, const TVector3 edir, const double eEnergy) const{
  if(!_hxscnD) return 0.0;
  if(eEnergy<_hxscnD->GetXaxis()->GetXmin() || eEnergy>_hxscnD->GetXaxis()->GetXmax())
    return 0.0;
  TVector3 nuDir = (epos-_neutrinoSourcePos);
  double costhe = nuDir.Dot(edir)/nuDir.Mag()/edir.Mag();
  return _hxscnD->Interpolate(eEnergy, costhe);
}

double G4d2oGeom::ParticleCosThe( const TVector3 epos, const TVector3 edir) const{
  TVector3 nuDir = (epos-_neutrinoSourcePos);
  double costhe = nuDir.Dot(edir)/nuDir.Mag()/edir.Mag();
  return costhe;
}

G4d2oGeom::~G4d2oGeom()
{
  //TODO: Recall how to clean up a singleton
}
