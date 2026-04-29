#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"

#include "G4PhysListFactory.hh"
#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"

int main(int argc, char** argv)
{
    auto* runManager =
        G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

    runManager->SetUserInitialization(new DetectorConstruction());

    G4PhysListFactory factory;
    runManager->SetUserInitialization(factory.GetReferencePhysList("FTFP_BERT_EMZ"));

    runManager->SetUserInitialization(new ActionInitialization());

    G4VisManager* visManager = nullptr;

    if (argc == 1) {
        auto* ui = new G4UIExecutive(argc, argv);
        visManager = new G4VisExecutive();
        visManager->Initialize();

        auto* UImanager = G4UImanager::GetUIpointer();
        UImanager->ApplyCommand("/control/execute macros/init_vis.mac");

        ui->SessionStart();
        delete ui;
    } else {
        auto* UImanager = G4UImanager::GetUIpointer();
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command + fileName);
    }

    delete visManager;
    delete runManager;
    return 0;
}
