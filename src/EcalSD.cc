#include "EcalSD.hh"

#include "G4HCofThisEvent.hh"
#include "G4LogicalVolume.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "G4VPhysicalVolume.hh"

EcalSD::EcalSD(const G4String& name, const G4String& hitsCollectionName)
  : G4VSensitiveDetector(name)
{
    collectionName.insert(hitsCollectionName);
}

void EcalSD::Initialize(G4HCofThisEvent* hce)
{
    fHitsCollection =
        new EcalHitsCollection(SensitiveDetectorName, collectionName[0]);

    if (fHCID < 0) {
        fHCID = GetCollectionID(0);
    }

    hce->AddHitsCollection(fHCID, fHitsCollection);
}

G4bool EcalSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
    const G4double edep = step->GetTotalEnergyDeposit();
    if (edep <= 0.) {
        return false;
    }

    auto* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
    G4String volName = volume->GetName();
    G4int copyID = volume->GetCopyNo();

    G4int regionID = -1;
    if (volName == "PWOBlock") {
        regionID = 0;
    } else if (volName == "SciGlassBlock") {
        regionID = 1;
    } else {
        return false;
    }

    auto* hit = new EcalHit();
    hit->SetCopyID(copyID);
    hit->AddEdep(edep);
    hit->SetRegionID(regionID);

    fHitsCollection->insert(hit);

    return true;
}

