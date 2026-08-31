#include "simEvent.h"
#include "TClonesArray.h"
#include "G4d2oGeom.h"
#include <vector>
#include <numeric>
#include <algorithm>
#include <functional>

using namespace std;

ClassImp(simEvent)

TClonesArray *simEvent::sPMTHits = 0;
TClonesArray *simEvent::sAreaPMTHits = 0;

simEvent::simEvent(Int_t maxPMT)
{
  if (!sPMTHits)
    sPMTHits = new TClonesArray("simHit", maxPMT);
  sPMTHits->SetOwner(true);
  pmtHits = sPMTHits;

  if (!sAreaPMTHits)
    sAreaPMTHits = new TClonesArray("simAreaHit", maxPMT);
  sAreaPMTHits->SetOwner(true);
  areaPMTHits = sAreaPMTHits;

  ClearData();
}

simEvent::~simEvent()
{
  pmtHits->Clear("C");
  areaPMTHits->Clear("C");
}

void simEvent::AddPMTHit(Int_t pNum, Double_t eTime, Double_t phEn)
{
  simHit *theHit = (simHit *)pmtHits->ConstructedAt(numHits++);
  theHit->Set(pNum, eTime, phEn);
}

void simEvent::AddAreaPMTHit(Int_t pNum, Double_t eTime, Double_t phEn, TVector3 hitPos)
{
  simAreaHit *theHit = (simAreaHit *)areaPMTHits->ConstructedAt(numHitsArea++);
  theHit->Set(pNum, eTime, phEn, hitPos);
}

void simEvent::ClearData()
{
  eventNumber = -1;
  direction0.SetXYZ(0.0, 0.0, 0.0);
  position0.SetXYZ(0.0, 0.0, 0.0);
  sourceParticleEnergy = 0;
  vol0 = 0;
  veto_tag = false;
  veto_edep = 0.0;

  for (int i=0;i<kNarrowPanels;++i) narrow_veto_edep[i] = 0.0;
  for (int i=0;i<kWidePanels;++i) wide_veto_edep[i] = 0.0;
  numHits = 0;
  pmtHits->Clear();

  // NEW: Clear instrumentation fields
  nReflections = 0;
  totalPathLength = 0.0;
  nPhotons = 0;
  meanReflections = 0.0;
  meanPathLength = 0.0;

  numHitsArea = 0;
  areaPMTHits->Clear();
  for (int i = 0; i < 12; i++)
    muVetoEnergy[i] = 0.0;
}

void simEvent::CopyData(simEvent *dataToCopy)
{
  eventNumber = dataToCopy->eventNumber;
  direction0 = TVector3(dataToCopy->direction0);
  position0 = TVector3(dataToCopy->position0);
  sourceParticleEnergy = dataToCopy->sourceParticleEnergy;
  vol0 = dataToCopy->vol0;

  numHits = dataToCopy->numHits;
  pmtHits = (TClonesArray *)dataToCopy->pmtHits->Clone();

  numHitsArea = dataToCopy->numHitsArea;
  areaPMTHits = (TClonesArray *)dataToCopy->areaPMTHits->Clone();
  for (int i = 0; i < 12; i++)
    muVetoEnergy[i] = dataToCopy->muVetoEnergy[i];

  // NEW: Copy instrumentation fields
  nReflections = dataToCopy->nReflections;
  totalPathLength = dataToCopy->totalPathLength;
  nPhotons = dataToCopy->nPhotons;
  meanReflections = dataToCopy->meanReflections;
  meanPathLength = dataToCopy->meanPathLength;
}

double simEvent::MeanX() const
{
  double x = 0.0;
  for (int ihit = 0; ihit < numHits; ihit++)
  {
    x += G4d2oGeom::Instance()->DetPos(GetHit(ihit)->pmtNum).X();
  }
  x /= numHits;
  return x;
}

double simEvent::MeanY() const
{
  double x = 0.0;
  for (int ihit = 0; ihit < numHits; ihit++)
  {
    x += G4d2oGeom::Instance()->DetPos(GetHit(ihit)->pmtNum).Y();
  }
  x /= numHits;
  return x;
}

double simEvent::TimeRMS() const
{
  std::vector<double> vtimes;
  vtimes.resize(numHits);
  for (int ihit = 0; ihit < numHits; ihit++)
  {
    vtimes[ihit] = GetHit(ihit)->eventTime;
  }
  double sum = std::accumulate(vtimes.begin(), vtimes.end(), 0.0);
  double mean = sum / vtimes.size();

  std::vector<double> diff(vtimes.size());
  std::transform(vtimes.begin(), vtimes.end(), diff.begin(), [mean](double x)
                 { return x - mean; });
  double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  double stdev = std::sqrt(sq_sum / vtimes.size());
  return stdev;
}

double simEvent::MeanZ() const
{
  double x = 0.0;
  for (int ihit = 0; ihit < numHits; ihit++)
  {
    x += G4d2oGeom::Instance()->DetPos(GetHit(ihit)->pmtNum).Z();
  }
  x /= numHits;
  return x;
}

double simEvent::WeightO() const
{
  return G4d2oGeom::Instance()->ProbO(position0, direction0, sourceParticleEnergy);
}

double simEvent::WeightD() const
{
  return G4d2oGeom::Instance()->ProbD(position0, direction0, sourceParticleEnergy);
}

double simEvent::SourceCosThe() const
{
  return G4d2oGeom::Instance()->ParticleCosThe(position0, direction0);
}

int simEvent::NumVetoPairs(double threshold /*MeV*/) const
{
  int numpairs = 0;
  for (int i = 0; i < 12; i += 2)
  {
    if (muVetoEnergy[i] > threshold && muVetoEnergy[i + 1] > threshold)
      numpairs++;
  }
  return numpairs;
}
