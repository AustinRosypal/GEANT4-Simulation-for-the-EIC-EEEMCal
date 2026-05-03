#include "RunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"

RunAction::RunAction()
{
    auto* analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetVerboseLevel(1);

    analysisManager->SetDefaultFileType("root");
    analysisManager->SetFileName("ecal");

    analysisManager->CreateNtuple("events", "EEEMCal event response");
    analysisManager->CreateNtupleDColumn("totalEdep_MeV");
    analysisManager->CreateNtupleDColumn("pwoEdep_MeV");
    analysisManager->CreateNtupleDColumn("sciGlassEdep_MeV");
    analysisManager->CreateNtupleIColumn("nHits");
    analysisManager->CreateNtupleIColumn("initialRegion");
    analysisManager->CreateNtupleDColumn("impactX_mm");
    analysisManager->CreateNtupleDColumn("impactY_mm");
    analysisManager->CreateNtupleDColumn("impactR_mm");
    analysisManager->CreateNtupleDColumn("beamEnergy_GeV");
    analysisManager->FinishNtuple();
}

void RunAction::BeginOfRunAction(const G4Run*)
{
    G4RunManager::GetRunManager()->SetRandomNumberStore(false);
    auto* analysisManager = G4AnalysisManager::Instance();
    analysisManager->OpenFile();
}

void RunAction::EndOfRunAction(const G4Run*)
{
    auto* analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();
}

