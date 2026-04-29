#ifndef EcalSD_h
#define EcalSD_h 1

#include "EcalHit.hh"
#include "G4VSensitiveDetector.hh"

class G4HCofThisEvent;
class G4Step;

class EcalSD : public G4VSensitiveDetector
{
  public:
    EcalSD(const G4String& name, const G4String& hitsCollectionName);
    ~EcalSD() override = default;

    void Initialize(G4HCofThisEvent* hce) override;
    G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;

  private:
    EcalHitsCollection* fHitsCollection = nullptr;
    G4int fHCID = -1;
};

#endif
