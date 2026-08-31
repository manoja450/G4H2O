#ifndef G4d2oStackingAction_hh
#define G4d2oStackingAction_hh 1

#include "G4UserStackingAction.hh"

#include <vector>

class G4d2oDetector;

class G4d2oStackingAction : public G4UserStackingAction
{
public:
  G4d2oStackingAction(G4d2oDetector * det);
  ~G4d2oStackingAction();

  G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track * the_track);

  void NewStage();
  void PrepareNewEvent();

private:
  G4d2oDetector * fDetector;
  std::vector<int> pmt_instance_ids;
};

#endif
