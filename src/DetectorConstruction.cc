#include "DetectorConstruction.hh"
#include "EcalSD.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ThreeVector.hh"
#include "G4ios.hh"

#include <cmath>
#include <algorithm>

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    auto* nist = G4NistManager::Instance();

    //
    // Materials
    //
    auto* worldMat  = nist->FindOrBuildMaterial("G4_AIR");
    auto* carbonMat = nist->FindOrBuildMaterial("G4_C");

    // PbWO4
    auto* elPb = nist->FindOrBuildElement("Pb");
    auto* elW  = nist->FindOrBuildElement("W");
    auto* elO  = nist->FindOrBuildElement("O");
    auto* elBa = nist->FindOrBuildElement("Ba");
    auto* elSi = nist->FindOrBuildElement("Si");
    auto* elCe = nist->FindOrBuildElement("Ce");

    auto* pbwo4 = new G4Material("PbWO4", 8.28 * g/cm3, 3);
    pbwo4->AddElement(elPb, 1);
    pbwo4->AddElement(elW,  1);
    pbwo4->AddElement(elO,  4);

    // SciGlass surrogate:
    // Public sources say SciGlass is a Ce-based, barium-silicate-based scintillating glass,
    // but do not provide the exact full production stoichiometry publicly.
    // This is therefore an effective showering material for geometry studies.
    auto* sciGlass = new G4Material("SciGlassEffective", 4.2 * g/cm3, 4);
    sciGlass->AddElement(elBa, 0.45);
    sciGlass->AddElement(elSi, 0.18);
    sciGlass->AddElement(elO,  0.36);
    sciGlass->AddElement(elCe, 0.01);

    //
    // Dimensions
    //
    const G4double innerR    = fInnerRadius    * mm;
    const G4double pwoOuterR = fPWOOuterRadius * mm;
    const G4double outerR    = fOuterRadius    * mm;
    const G4double ecalZ     = fEEEMCalZ       * mm;

    const G4double pwoXY     = fPWOCellXY * mm;
    const G4double pwoZ      = fPWOLength * mm;

    const G4double sciXY     = fSciCellXY * mm;
    const G4double sciZ      = fSciLength * mm;

    const G4double gapXY     = fGapXY * mm;
    const G4double frameT    = fCarbonFrameThickness * mm;

    // Common front face. Put both radiators behind the same front plane.
    const G4double maxDepth = std::max(pwoZ, sciZ);
    const G4double envelopeThickness = maxDepth + 2.0 * frameT;

    const G4double worldXY = 2.2 * outerR;
    const G4double worldZ  = 5.0 * m;

    //
    // World
    //
    auto* solidWorld = new G4Box("World", worldXY/2., worldXY/2., worldZ/2.);
    auto* logicWorld = new G4LogicalVolume(solidWorld, worldMat, "World");
    auto* physWorld =
        new G4PVPlacement(nullptr, {}, logicWorld, "World", nullptr, false, 0, true);

    //
    // Mother volume
    //
    auto* solidEEEMCal =
        new G4Tubs("EEEMCalMother",
                   innerR, outerR,
                   envelopeThickness/2.,
                   0.0*deg, 360.0*deg);

    auto* logicEEEMCal =
        new G4LogicalVolume(solidEEEMCal, worldMat, "EEEMCalMother");

    new G4PVPlacement(nullptr,
                      G4ThreeVector(0., 0., ecalZ),
                      logicEEEMCal,
                      "EEEMCalMother",
                      logicWorld,
                      false,
                      0,
                      true);

    //
    // Front / back support plates
    //
    auto* solidFrame =
        new G4Tubs("CarbonFrame",
                   innerR, outerR,
                   frameT/2.,
                   0.0*deg, 360.0*deg);

    auto* logicFrame =
        new G4LogicalVolume(solidFrame, carbonMat, "CarbonFrame");

    new G4PVPlacement(nullptr,
                      G4ThreeVector(0., 0., +envelopeThickness/2. - frameT/2.),
                      logicFrame,
                      "CarbonFrameFront",
                      logicEEEMCal,
                      false,
                      0,
                      true);

    new G4PVPlacement(nullptr,
                      G4ThreeVector(0., 0., -envelopeThickness/2. + frameT/2.),
                      logicFrame,
                      "CarbonFrameBack",
                      logicEEEMCal,
                      false,
                      1,
                      true);

    //
    // Logical volumes
    //
    auto* solidPWO = new G4Box("PWOBlock", pwoXY/2., pwoXY/2., pwoZ/2.);
    auto* solidSci = new G4Box("SciGlassBlock", sciXY/2., sciXY/2., sciZ/2.);

    fPWOLogical      = new G4LogicalVolume(solidPWO, pbwo4,   "PWOLogical");
    fSciGlassLogical = new G4LogicalVolume(solidSci, sciGlass, "SciGlassLogical");

    //
    // Front-face alignment:
    // local +z is toward the IP-facing side of the mother volume
    //
    const G4double frontFaceZ = +envelopeThickness/2. - frameT;
    const G4double pwoCenterZ = frontFaceZ - pwoZ/2.;
    const G4double sciCenterZ = frontFaceZ - sciZ/2.;

    //
    // Inner PbWO4 annulus
    //
    const G4double pwoPitch = pwoXY + gapXY;
    const G4int nP = static_cast<G4int>((2.0 * pwoOuterR) / pwoPitch) + 4;

    G4int copyNo = 0;
    G4int nPWO   = 0;
    for (G4int ix = -nP/2; ix <= nP/2; ++ix) {
        for (G4int iy = -nP/2; iy <= nP/2; ++iy) {
            const G4double x = ix * pwoPitch;
            const G4double y = iy * pwoPitch;
            const G4double r = std::sqrt(x*x + y*y);

            if (r < innerR)    continue;
            if (r > pwoOuterR) continue;

            new G4PVPlacement(nullptr,
                              G4ThreeVector(x, y, pwoCenterZ),
                              fPWOLogical,
                              "PWOBlock",
                              logicEEEMCal,
                              false,
                              copyNo,
                              true);
            ++copyNo;
            ++nPWO;
        }
    }

    //
    // Outer SciGlass annulus
    //
    const G4double sciPitch = sciXY + gapXY;
    const G4int nS = static_cast<G4int>((2.0 * outerR) / sciPitch) + 4;

    G4int nSci = 0;
    for (G4int ix = -nS/2; ix <= nS/2; ++ix) {
        for (G4int iy = -nS/2; iy <= nS/2; ++iy) {
            const G4double x = ix * sciPitch;
            const G4double y = iy * sciPitch;
            const G4double r = std::sqrt(x*x + y*y);

            if (r < pwoOuterR) continue;
            if (r > outerR)    continue;

            new G4PVPlacement(nullptr,
                              G4ThreeVector(x, y, sciCenterZ),
                              fSciGlassLogical,
                              "SciGlassBlock",
                              logicEEEMCal,
                              false,
                              copyNo,
                              true);
            ++copyNo;
            ++nSci;
        }
    }

    G4cout << "Placed " << nPWO << " PbWO4 blocks and "
           << nSci << " SciGlass blocks." << G4endl;

    //
    // Visualization
    //
    logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());

    auto* motherVis = new G4VisAttributes(G4Colour(0.8, 0.8, 0.8, 0.05));  // Fourth param is alpha
    motherVis->SetVisibility(false);
    logicEEEMCal->SetVisAttributes(motherVis);

    //auto* pwoVis = new G4VisAttributes(G4Colour(0.15, 0.55, 0.95));
    auto* pwoVis = new G4VisAttributes(G4Colour(1.0,0.5,0.0));
    pwoVis->SetForceSolid(true);
    fPWOLogical->SetVisAttributes(pwoVis);

    //auto* sciVis = new G4VisAttributes(G4Colour(0.15, 0.85, 0.35));
    auto* sciVis = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
    sciVis->SetForceSolid(true);
    fSciGlassLogical->SetVisAttributes(sciVis);

    auto* frameVis = new G4VisAttributes(G4Colour(0.3, 0.3, 0.3));
    //frameVis->SetForceSolid(true);
    frameVis->SetLineWidth(2.0);
    frameVis->SetForceWireframe(true);
    logicFrame->SetVisAttributes(frameVis);

    return physWorld;
}

void DetectorConstruction::ConstructSDandField()
{
    auto* sdManager = G4SDManager::GetSDMpointer();
    auto* ecalSD = new EcalSD("EcalSD", "EcalHitsCollection");
    sdManager->AddNewDetector(ecalSD);

    // Preserve the same energy-deposition method in both regions
    SetSensitiveDetector(fPWOLogical, ecalSD);
    SetSensitiveDetector(fSciGlassLogical, ecalSD);
}
