#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction() = default;
    ~DetectorConstruction() override = default;

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

  private:
    G4LogicalVolume* fPWOLogical      = nullptr;
    G4LogicalVolume* fSciGlassLogical = nullptr;

    // ---- Publicly motivated dimensions ----
    // PbWO4 inner region: 20 mm transverse, 200 mm long
    const G4double fPWOCellXY = 20.0;   // mm
    const G4double fPWOLength = 200.0;  // mm

    // SciGlass outer region: 40 mm transverse
    // Use 455 mm as a practical 17-ish X0 scale surrogate from ECCE-era studies
    const G4double fSciCellXY = 40.0;   // mm
    const G4double fSciLength = 455.0;  // mm

    // Small packing gaps
    const G4double fGapXY = 0.5;        // mm

    // Annular envelope
    const G4double fInnerRadius    = 95.0;   // mm  beam pipe hole
    const G4double fPWOOuterRadius = 260.0;  // mm  tunable boundary
    const G4double fOuterRadius    = 650.0;  // mm  detector outer edge

    // Center position of endcap
    const G4double fEEEMCalZ = -2000.0;      // mm

    // Simple front/back support
    const G4double fCarbonFrameThickness = 0.5; // mm
};

#endif
