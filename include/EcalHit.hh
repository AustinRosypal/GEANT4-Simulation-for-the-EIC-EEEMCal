#ifndef EcalHit_h
#define EcalHit_h 1

#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4VHit.hh"
#include "globals.hh"

class EcalHit : public G4VHit
{
  public:
    EcalHit() = default;
    ~EcalHit() override = default;

    EcalHit(const EcalHit&) = default;
    EcalHit& operator=(const EcalHit&) = default;

    G4bool operator==(const EcalHit& right) const { return (this == &right); }

    inline void* operator new(size_t);
    inline void operator delete(void*);

    void SetCopyID(G4int value)      { fCopyID = value; }
    void AddEdep(G4double value)     { fEdep += value; }
    void SetRegionID(G4int value)    { fRegionID = value; }

    G4int GetCopyID() const          { return fCopyID; }
    G4double GetEdep() const         { return fEdep; }
    G4int GetRegionID() const        { return fRegionID; }

    void Print() override;

  private:
    G4int    fCopyID   = -1;
    G4double fEdep     = 0.;
    G4int    fRegionID = -1; // 0 = PbWO4, 1 = SciGlass
};

using EcalHitsCollection = G4THitsCollection<EcalHit>;

extern G4ThreadLocal G4Allocator<EcalHit>* EcalHitAllocator;

inline void* EcalHit::operator new(size_t)
{
    if (!EcalHitAllocator) {
        EcalHitAllocator = new G4Allocator<EcalHit>;
    }
    return EcalHitAllocator->MallocSingle();
}

inline void EcalHit::operator delete(void* hit)
{
    EcalHitAllocator->FreeSingle(static_cast<EcalHit*>(hit));
}

#endif


/*#ifndef EcalHit_h
#define EcalHit_h 1

#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4VHit.hh"
#include "globals.hh"

class EcalHit : public G4VHit
{
  public:
    EcalHit() = default;
    ~EcalHit() override = default;

    EcalHit(const EcalHit&) = default;
    EcalHit& operator=(const EcalHit&) = default;

    G4bool operator==(const EcalHit& right) const { return (this == &right); }

    inline void* operator new(size_t);
    inline void operator delete(void*);

    void SetLayerID(G4int value) { fLayerID = value; }
    void AddEdep(G4double value) { fEdep += value; }

    G4int GetLayerID() const { return fLayerID; }
    G4double GetEdep() const { return fEdep; }

    void Print() override;

  private:
    G4int    fLayerID = -1;
    G4double fEdep    = 0.;
};

using EcalHitsCollection = G4THitsCollection<EcalHit>;

extern G4ThreadLocal G4Allocator<EcalHit>* EcalHitAllocator;

inline void* EcalHit::operator new(size_t)
{
    if (!EcalHitAllocator) {
        EcalHitAllocator = new G4Allocator<EcalHit>;
    }
    return EcalHitAllocator->MallocSingle();
}

inline void EcalHit::operator delete(void* hit)
{
    EcalHitAllocator->FreeSingle(static_cast<EcalHit*>(hit));
}

#endif */
