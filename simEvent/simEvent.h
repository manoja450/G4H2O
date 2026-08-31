#ifndef __simEvent_H__
#define __simEvent_H__

#include "TVector3.h"
#include "TClonesArray.h"

class simHit : public TObject
{
public:
    simHit() { ; }
    simHit(Int_t pNum, Double_t eTime, Double_t phEn) { Set(pNum, eTime, phEn); }

    void Set(Int_t pNum, Double_t eTime, Double_t phEn)
    {
        pmtNum = pNum;
        eventTime = eTime;
        photonEnergy = phEn;
    }

    Int_t pmtNum;
    Double_t eventTime;
    Double_t photonEnergy;

    ClassDef(simHit, 1)
};

class simAreaHit : public simHit
{
public:
    simAreaHit() { ; }
    simAreaHit(Int_t pNum, Double_t eTime, Double_t phEn, TVector3 hitPos) : simHit(pNum, eTime, phEn) { Set(hitPos); }

    void Set(TVector3 hitPos)
    {
        photonPosition = hitPos;
    }

    void Set(Int_t pNum, Double_t eTime, Double_t phEn, TVector3 hitPos)
    {
        pmtNum = pNum;
        eventTime = eTime;
        photonEnergy = phEn;
        photonPosition = hitPos;
    }

    TVector3 photonPosition;

    ClassDef(simAreaHit, 1)
};

class simEvent : public TObject
{
public:
    simEvent(Int_t maxHits = 10000);
    virtual ~simEvent();

    static const int kNarrowPanels = 4;
    double narrow_veto_edep[kNarrowPanels] = {0, 0, 0, 0};
    static const int kWidePanels = 4;
    double wide_veto_edep[kWidePanels] = {0, 0, 0, 0};
    Int_t eventNumber;
    TVector3 direction0;
    TVector3 position0;
    int vol0;
    bool veto_tag{false};
    double veto_edep{0.0};

    Double_t sourceParticleEnergy;

    Int_t numHits;
    Int_t numHitsArea;

    Double_t muVetoEnergy[12];

    TClonesArray *pmtHits;
    static TClonesArray *sPMTHits;

    TClonesArray *areaPMTHits;
    static TClonesArray *sAreaPMTHits;

    // ============================================================
    // NEW: Photon instrumentation fields
    // ============================================================
    Int_t nReflections;          // Total reflections in event
    Double_t totalPathLength;    // Total path length of all photons (mm)
    Int_t nPhotons;              // Number of optical photons tracked
    Double_t meanReflections;    // Average reflections per photon
    Double_t meanPathLength;     // Average path length per photon (mm)

    void AddPMTHit(Int_t pNum, Double_t eTime, Double_t phEn);
    void AddAreaPMTHit(Int_t pNum, Double_t eTime, Double_t phEn, TVector3 hitPos);
    double MeanX() const;
    double MeanY() const;
    double MeanZ() const;
    double TimeRMS() const;
    double WeightO() const;
    double WeightD() const;
    double SourceCosThe() const;
    int NumVetoPairs(double threshold = 5.0 /*MeV*/) const;

    const simHit *GetHit(int ihit) const { return dynamic_cast<simHit *>(pmtHits->At(ihit)); }
    const simAreaHit *GetAHit(int ihit) const { return dynamic_cast<simAreaHit *>(areaPMTHits->At(ihit)); }

    virtual void ClearData();
    virtual void CopyData(simEvent *dataToCopy);
    ClassDef(simEvent, 7)
};

#endif
