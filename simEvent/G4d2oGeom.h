#pragma once

#include "TVector3.h"
#include "TObject.h"

class TH2;

class G4d2oGeom : public TObject
{
public:
  virtual ~G4d2oGeom();
  static G4d2oGeom* Instance(){
    if(!_geom) _geom = new G4d2oGeom();
    return _geom;
  }
  
  TVector3& DetPos(int ipmt) { return _detpos[ipmt]; }
  double ProbO( const TVector3 epos, const TVector3 edir, const double eEnergy) const;
  double ProbD( const TVector3 epos, const TVector3 edir, const double eEnergy) const;
  double ParticleCosThe( const TVector3 epos, const TVector3 edir) const;
  
private:
  G4d2oGeom();
  static G4d2oGeom* _geom;
  TVector3 _neutrinoSourcePos;
  std::vector<TVector3> _detpos;
  TH2* _hxscnD = 0;
  TH2* _hxscnO = 0;
  
  ClassDef(G4d2oGeom,0)

  
};

