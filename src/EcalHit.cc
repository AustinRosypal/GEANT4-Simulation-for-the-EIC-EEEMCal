#include "EcalHit.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"

G4ThreadLocal G4Allocator<EcalHit>* EcalHitAllocator = nullptr;

void EcalHit::Print()
{
    G4cout << "Copy " << fCopyID
           << "  Region " << fRegionID
           << "  Edep = " << G4BestUnit(fEdep, "Energy")
           << G4endl;
}

