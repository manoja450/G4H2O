#ifndef inputVariables_H
#define inputVariables_H 1

#include "globals.hh"
#include <fstream>

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TStopwatch.h"

#define NUMINT 13
#define NUMDOUBLE 6
#define NUMCHAR 1
#define NUMLONG 1

#define MAX_SEN_DET 20

class GenericInputHolder
{
public:
    GenericInputHolder() {;}

    Int_t iPar[NUMINT];
    Char_t iDesc[NUMINT][100];

    Long64_t lPar[NUMLONG];
    Char_t lDesc[NUMLONG][100];

    Double_t dPar[NUMDOUBLE];
    Char_t dDesc[NUMDOUBLE][100];

    TString sPar[NUMCHAR];
    Char_t sDesc[NUMCHAR][100];
};

using namespace std;

class inputVariables : public TObject
{
public:
    static inputVariables * ptrIV;
    static ofstream * ptrOF;

    inputVariables(int argc=0, char** argv=0);
    ~inputVariables();

    static inputVariables* GetIVPointer( int argc=0, char** argv=0 );
    ofstream* GetOutputFile( void );

    inline void SetNumberOfEvents(G4int numEntries) {nevents = numEntries;}

    inline G4int GetRunNumber() {return runno;}
    inline G4int GetNumberOfEvents() {return nevents;}
    inline G4int GetVisualization() {return ivis;}
    inline G4int GetPrintStatus() {return iprint;}
    inline G4int GetPGAType() {return iPGA;}
    inline G4int GetRandomStatus() {return irand;}
    inline G4int GetPhysicsType() {return iPhysics;}
    inline G4int GetNeutronHP() {return iNeutronHP;}
    inline G4int GetSideLining() {return iSideLining;}
    inline G4int GetBottomPMTs() {return iUseBottomPMTs;}
    inline G4int GetBottomVeto() {return iUseBottomVeto;}
    inline G4int GetBottomShielding() {return iUseBottomShielding;}
    inline G4int GetPMTQE() {return iPMTQE;}
    inline G4double GetPMTDiameter() {return pmtDiameter;}
    inline G4long GetUserSEED() {return userSEED;}
    inline TString GetOutputDirectory() {return outputDir;}
    inline G4double GetTailCatcherThick() {return tailThick;}
    inline G4double GetShieldThickness() {return shieldThick;}
    inline G4double GetH2oRefl() {return h2oRefl;}
    inline G4double GetReflectivity() {return reflectivity;}
    inline G4double GetSpecialEnergy() {return specialEnergy;}

    Bool_t ReadCurrentFile();
    Bool_t UpdateValues(int argc, char** argv);
    void WriteValues(GenericInputHolder *theVals, Bool_t bFile=true, TString sBOFileName="beamOn.dat");
    void SetValues(GenericInputHolder *theVals);
    void PrintHelp(Char_t *exeName);
    void SetupDirectory();
    void PrintCommandLine(char** argv, GenericInputHolder *theVals);

    G4int numHitsColl, numDetHC, numSrcHC, numTargHC;
    Int_t detCollID[MAX_SEN_DET];
    Int_t srcCollID[MAX_SEN_DET];
    Int_t targCollID[MAX_SEN_DET];
    G4int GetNewCollectionID(G4int iType);
    G4int GetNumHitsCollections() {return numHitsColl;}

private:
    ofstream outfile;
    TObjArray *thePTracFiles;

    G4int runno, nevents, ivis, iprint, irand, iPGA;
    G4int iPhysics, iNeutronHP, iSideLining, iUseBottomPMTs, iPMTQE, iUseBottomVeto, iUseBottomShielding;
    G4int isave;
    G4double pmtDiameter, tailThick, shieldThick, h2oRefl, reflectivity, specialEnergy;
    G4long userSEED;
    TString outputDir;

    TDatime *dateTime;
    GenericInputHolder *currentVals;
    GenericInputHolder *newVals;

    const Int_t indexRunNum = 0;
    const Int_t indexNumEvents = 1;
    const Int_t indexVis = 2;
    const Int_t indexDebug = 3;
    const Int_t indexSeed = 4;
    const Int_t indexPhysics = 5;
    const Int_t indexNeutronHP = 6;
    const Int_t indexPGA = 7;
    const Int_t indexSideLining = 8;
    const Int_t indexUseBottomPMTs = 9;
    const Int_t indexPMTQE = 10;
    const Int_t indexUseBottomVeto = 11;
    const Int_t indexUseBottomShielding = 12;

    const Int_t indexPMTDiam = 0;
    const Int_t indexTailThick = 1;
    const Int_t indexShieldThick = 2;
    const Int_t indexH2oRefl = 3;
    const Int_t indexReflectivity = 4;
    const Int_t indexSpecialEnergy = 5;

    const Int_t indexOutDir = 0;
    const Int_t indexUserSEED = 0;
};

#endif

