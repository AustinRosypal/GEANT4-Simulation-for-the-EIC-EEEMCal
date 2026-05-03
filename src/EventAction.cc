#include "EventAction.hh"
#include "EcalHit.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"

#include <cmath>

void EventAction::BeginOfEventAction(const G4Event*)
{
    fTotalEdep    = 0.;
    fPWOEdep      = 0.;
    fSciGlassEdep = 0.;
    fNHits        = 0;
}

void EventAction::EndOfEventAction(const G4Event* event)
{
    auto* hce = event->GetHCofThisEvent();
    if (!hce) {
        return;
    }

    static G4int hcID = -1;
    if (hcID < 0) {
        hcID = G4SDManager::GetSDMpointer()->GetCollectionID("EcalHitsCollection");
    }

    auto* hitsCollection =
        static_cast<EcalHitsCollection*>(hce->GetHC(hcID));

    if (!hitsCollection) {
        return;
    }

    const G4int nHits = hitsCollection->entries();
    fNHits = nHits;

    for (G4int i = 0; i < nHits; ++i) {
        G4double edep = (*hitsCollection)[i]->GetEdep();
        G4int region  = (*hitsCollection)[i]->GetRegionID();

        fTotalEdep += edep;

        if (region == 0) {
            fPWOEdep += edep;
        } else if (region == 1) {
            fSciGlassEdep += edep;
        }
    }

    G4int initialRegion = -1;
    G4double x0_mm = 0.0;
    G4double y0_mm = 0.0;
    G4double r0_mm = 0.0;
    G4double beamEnergy_GeV = 0.0;

    auto* pv = event->GetPrimaryVertex(0);
    if (pv) {
        const G4double x0 = pv->GetX0();
        const G4double y0 = pv->GetY0();
        const G4double r0 = std::sqrt(x0*x0 + y0*y0);

        x0_mm = x0 / mm;
        y0_mm = y0 / mm;
        r0_mm = r0 / mm;

        const G4double pwoOuterR = 260.0 * mm; // must match DetectorConstruction
        if (r0 < pwoOuterR) {
            initialRegion = 0;
        } else {
            initialRegion = 1;
        }

        auto* pp = pv->GetPrimary();
        if (pp) {
            beamEnergy_GeV = pp->GetKineticEnergy() / GeV;
        }
    }

    auto* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleDColumn(0, fTotalEdep / MeV);
    analysisManager->FillNtupleDColumn(1, fPWOEdep / MeV);
    analysisManager->FillNtupleDColumn(2, fSciGlassEdep / MeV);
    analysisManager->FillNtupleIColumn(3, fNHits);
    analysisManager->FillNtupleIColumn(4, initialRegion);
    analysisManager->FillNtupleDColumn(5, x0_mm);
    analysisManager->FillNtupleDColumn(6, y0_mm);
    analysisManager->FillNtupleDColumn(7, r0_mm);
    analysisManager->FillNtupleDColumn(8, beamEnergy_GeV);
    analysisManager->AddNtupleRow();
}

