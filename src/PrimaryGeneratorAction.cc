#include "PrimaryGeneratorAction.hh"

#include "G4Electron.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

#include <cmath>

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
    fParticleGun = new G4ParticleGun(1);

    fParticleGun->SetParticleDefinition(G4Electron::Definition());
    fParticleGun->SetParticleEnergy(10.0 * GeV);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
    const G4double innerR    = 95.0  * mm;
    const G4double outerR    = 650.0 * mm;

    const G4double zStart = 50.0 * cm;
    const G4double zFace  = -200.0 * cm;

    const G4double u   = G4UniformRand();
    const G4double phi = CLHEP::twopi * G4UniformRand();

    const G4double r = std::sqrt(innerR*innerR +
                         u * (outerR*outerR - innerR*innerR));

    const G4double x = r * std::cos(phi);
    const G4double y = r * std::sin(phi);

    const G4ThreeVector startPos(x, y, zStart);
    const G4ThreeVector target(x, y, zFace);
    G4ThreeVector dir = (target - startPos).unit();

    fParticleGun->SetParticlePosition(startPos);
    fParticleGun->SetParticleMomentumDirection(dir);

    fParticleGun->GeneratePrimaryVertex(event);
}

/*#include "PrimaryGeneratorAction.hh"

#include "G4Electron.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

#include <cmath>

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
    fParticleGun = new G4ParticleGun(1);

    fParticleGun->SetParticleDefinition(G4Electron::Definition());
    fParticleGun->SetParticleEnergy(10.0 * GeV);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
    // Match the approximate EEEMCal annular geometry
    const G4double innerR = 95.0  * mm;
    const G4double outerR = 650.0 * mm;

    // Start upstream of the backward endcap and shoot toward -z
    const G4double zStart = 50.0 * cm;
    const G4double zFace  = -200.0 * cm;  // approx EEEMCal center plane

    // Sample uniformly in area over the annulus
    const G4double u   = G4UniformRand();
    const G4double phi = CLHEP::twopi * G4UniformRand();

    const G4double r = std::sqrt(innerR*innerR +
                         u * (outerR*outerR - innerR*innerR));

    const G4double x = r * std::cos(phi);
    const G4double y = r * std::sin(phi);

    // Start on the same x,y line so the particle heads straight into the chosen point
    const G4ThreeVector startPos(x, y, zStart);
    const G4ThreeVector target(x, y, zFace);

    G4ThreeVector dir = (target - startPos).unit();

    fParticleGun->SetParticlePosition(startPos);
    fParticleGun->SetParticleMomentumDirection(dir);

    fParticleGun->GeneratePrimaryVertex(event);
} */
