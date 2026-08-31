#ifndef G4d2oRunAction_hh
#define G4d2oRunAction_hh 1

#include <fstream>
#include "globals.hh"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TStopwatch.h"
#include "TMacro.h"
#include "TVector3.h"
#include "G4UserRunAction.hh"
#include "G4Run.hh"
#include "G4d2oMaterialsDefinition.hh"
#include "inputVariables.hh"
#include "ReplayTools.h"

class G4d2oRunAction : public G4UserRunAction
{
private:
    static G4d2oMaterialsDefinition *materialsPtr;
    inputVariables *input;
    G4int runno;
    Int_t numPrimaryEvents;
    struct Date_t { UInt_t month, day, year; };
    Date_t runDate;
    Char_t path[1000];
    Char_t name[100];
    TString outFileName;
    Int_t pmtRows;
    Int_t totPMTs;
    TVector3 tankSize;
    TTree *setupTree;
    TMacro mBeamOn;
    ReplayTools *theRT;

protected:

public:
    G4d2oRunAction();
    ~G4d2oRunAction();
    void BeginOfRunAction( const G4Run *aRun);
    void EndOfRunAction( const G4Run*);
    Char_t * GetFilePath() {return path;}
    Char_t * GetFileName() {return name;}
    TString GetOutFileName() {return outFileName;}
    static G4d2oMaterialsDefinition* GetMaterialsPointer( void );
};

#endif
