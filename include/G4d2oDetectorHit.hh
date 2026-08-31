#ifndef G4d2oDetectorHit_h
#define G4d2oDetectorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4Allocator.hh"

class G4d2oDetectorHit : public G4VHit
{
private:
    
    G4int hitPixelNumber;    // pixel number
    G4int hitDetectorNumber; // detector number
    G4int hitTrackID;        // ID Number of the track that hit the detector
    G4int hitTrackParentID;  // ID Number of the parent track for the particle
    G4int hitScatA;          // A value of G4d2oing nucleus
    G4int hitScatZ;          // Z value of G4d2oing nucleus
    
    G4double hitTotalEnergyDeposit; // Energy deposited in det for a single hit
    G4double hitGlobalTime;         // Kinetic energy of particle that hit the detector
    G4double hitKineticEnergy;      // Kinetic energy of particle that hit the detector
    G4double hitPDGCharge;          // Charge of particle
    G4double hitBaryonNumber;       // Baryon Number of particle
    G4double hitLeptonNumber;       // Lepton Number of particle
    G4double hitDeltaTime;          // Time between pre- and post-step-points for particle that hit the detector
    G4double hitCalcA;              // Time between pre- and post-step-points for particle that hit the detector
    
    char hitParticleName[100];     // Name of particle that hit the detector
    char hitProcessName[100];      // Name of process in the pixel
    char hitProcessTypeName[100];  // Name of process type in the pixel
    char hitMaterialName[100];     // Name of material in which the hit occurred
    char hitCreatorProcessName[100];  // Name of process which created particle
    char hitCreatorProcessType[100];  // Name of process type which created particle
    
    G4ThreeVector hitPosition;
    G4ThreeVector hitPrePosition;
    G4ThreeVector hitVertex;
    G4ThreeVector hitMomentumDirection;
    G4ThreeVector hitPostMomentumDirection;
    
protected:
    
    
public:
    
    G4d2oDetectorHit();
    ~G4d2oDetectorHit();
    
    inline void *operator new(size_t);
    inline void operator delete(void *aHit);
    
    void Draw(void);
    void Print(void);
    
    void SetPixelNumber( G4int thePixelNumber );
    G4int GetPixelNumber( void );
    void SetDetectorNumber( G4int theDetectorNumber );
    G4int GetDetectorNumber( void );
    void SetTrackID( G4int theTrackID );
    G4int GetTrackID( void );
    void SetTrackParentID( G4int theTrackParentID );
    G4int GetTrackParentID( void );
    void SetScatA( G4int theScatA );
    G4int GetScatA( void );
    void SetScatZ( G4int theScatZ );
    G4int GetScatZ( void );
    
    void SetTotalEnergyDeposit( G4double theEnergyDeposit );
    G4double GetTotalEnergyDeposit( void );
    void SetGlobalTime( G4double theGlobalTime );
    G4double GetGlobalTime( void );
    void SetKineticEnergy( G4double theKineticEnergy );
    G4double GetKineticEnergy( void );
    void SetPDGCharge( G4double thePDGCharge );
    G4double GetPDGCharge( void );
    void SetBaryonNumber( G4double theBaryonNumber );
    G4double GetBaryonNumber( void );
    void SetLeptonNumber( G4double theLeptonNumber );
    G4double GetLeptonNumber( void );
    void SetDeltaTime( G4double theDeltaTime );
    G4double GetDeltaTime( void );
    void SetCalcA( G4double theCalcA );
    G4double GetCalcA( void );
    
    void SetParticleName( char *theParticleName );
    char* GetParticleName( void );
    void SetProcessName( char *theProcessName );
    char* GetProcessName( void );
    void SetProcessTypeName( char *theProcessTypeName );
    char* GetProcessTypeName( void );
    void SetMaterialName( char *theMaterialName );
    char* GetMaterialName( void );
    void SetCreatorProcessName( char *theCreatorProcessName );
    char* GetCreatorProcessName( void );
    void SetCreatorProcessType( char *theCreatorProcessType );
    char* GetCreatorProcessType( void );
    
    void SetMomentumDirection( G4ThreeVector theMomentumDirection );
    G4ThreeVector GetMomentumDirection( void );
    void SetPostMomentumDirection( G4ThreeVector thePostMomentumDirection );
    G4ThreeVector GetPostMomentumDirection( void );
    void SetPosition( G4ThreeVector thePosition );
    G4ThreeVector GetPosition( void );
    void SetPrePosition( G4ThreeVector thePrePosition );
    G4ThreeVector GetPrePosition( void );
    void SetVertex( G4ThreeVector theVertex );
    G4ThreeVector GetVertex( void );
    
}; //END of class G4d2oDetectorHit

// Templated hits collection
typedef G4THitsCollection<G4d2oDetectorHit> G4d2oDetectorHitsCollection;

extern G4Allocator<G4d2oDetectorHit> G4d2oDetectorHitAllocator;

inline void* G4d2oDetectorHit::operator new(size_t)
{
    void* aHit;
    aHit = (void*)G4d2oDetectorHitAllocator.MallocSingle();
    return aHit;
}

inline void G4d2oDetectorHit::operator delete(void* aHit)
{
    G4d2oDetectorHitAllocator.FreeSingle((G4d2oDetectorHit*) aHit);
}

#endif
