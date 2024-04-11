//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "Run.hh"
#include "HistoManager.hh"

#include "G4SteppingManager.hh"
#include "G4RunManager.hh"
#include "Randomize.hh"

#include <G4HadronicInteraction.hh>
#include "G4HadronicProcess.hh"
#include "G4SystemOfUnits.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4LogicalVolume.hh"
#include "G4IonTable.hh"
#include "G4Alpha.hh"
#include "G4ParticleTable.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, EventAction* event)
:G4UserSteppingAction(),fDetector(det), fEventAction(event)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  // count processes
  // 
  Run* run = static_cast<Run*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  const G4StepPoint* endPoint = step->GetPostStepPoint();
  //const G4VProcess* process   = endPoint->GetProcessDefinedStep();
  G4VProcess* process   = 
                   const_cast<G4VProcess*>(endPoint->GetProcessDefinedStep());
  
  run->CountProcesses(process);
  //add
  // check that an real interaction occured (eg. not a transportation)
  G4StepStatus stepStatus = endPoint->GetStepStatus();
  G4bool transmit = (stepStatus==fGeomBoundary || stepStatus==fWorldBoundary);
  if (transmit) return;
  //end

 
 G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();     

 //longitudinal profile of deposited energy
 //randomize point of energy deposition
 //
 G4StepPoint* prePoint  = step->GetPreStepPoint();
 G4StepPoint* postPoint = step->GetPostStepPoint(); 
 G4ThreeVector P1 = prePoint ->GetPosition();
 G4ThreeVector P2 = postPoint->GetPosition();
 G4ThreeVector point = P1 + G4UniformRand()*(P2 - P1);


//add
 //initialisation of the nuclear channel identification
  //
  G4double Q             = - prePoint->GetKineticEnergy();
  G4ParticleDefinition* particle = step->GetTrack()->GetDynamicParticle()->GetDefinition();
  G4String partName = particle->GetParticleName();
  G4String nuclearChannel = partName;
  G4HadronicProcess* hproc = dynamic_cast<G4HadronicProcess*>(process);
  const G4Isotope* target = NULL;
  if (hproc) target = hproc->GetTargetIsotope();
  G4String targetName = "XXXX";  
  if (target) targetName = target->GetName();
  nuclearChannel += " + " + targetName + " --> ";
  if (targetName == "XXXX") run->SetTargetXXX(true);
  
   //scattered primary particle (if any)
  if (step->GetTrack()->GetTrackStatus() == fAlive) {
    G4double energy = endPoint->GetKineticEnergy();      
    //
    Q        += energy;
    //
    nuclearChannel += partName + " + ";
  }  

//secondaries
  //
  const std::vector<const G4Track*>* secondary 
                                    = step->GetSecondaryInCurrentStep();  
  G4int nbOfAbsor   = fDetector->GetNbOfAbsor();
  for (size_t lp=0; lp<(*secondary).size(); lp++) {
      for (G4int k=1; k<= nbOfAbsor; k++) {
    particle = (*secondary)[lp]->GetDefinition(); 
    G4String name   = particle->GetParticleName();
    G4String type   = particle->GetParticleType();      
    G4double energy = (*secondary)[lp]->GetKineticEnergy();
    //run->ParticleCount(k,name,energy,0);
    Q        += energy;
    //count e- from internal conversion together with gamma

    if (particle == G4Electron::Electron()) particle = G4Gamma::Gamma();
    //particle flag
    fParticleFlag[particle]++;
  }
  }

  //end


 // add
 // nuclear channel
  const G4int kMax = 16;  
  const G4String conver[] = {"0","","2 ","3 ","4 ","5 ","6 ","7 ","8 ","9 ",
                             "10 ","11 ","12 ","13 ","14 ","15 ","16 "};
  std::map<G4ParticleDefinition*,G4int>::iterator ip;               
  for (ip = fParticleFlag.begin(); ip != fParticleFlag.end(); ip++) {
    particle = ip->first;
    G4String name = particle->GetParticleName();      
    G4int nb = ip->second;
    if (nb > kMax) nb = kMax;   
    G4String Nb = conver[nb];    
    if (particle == G4Gamma::Gamma()) {
     run->CountGamma(nb);
     Nb = "N ";
     name = "gamma or e-";
    } 
    if (ip != fParticleFlag.begin()) nuclearChannel += " + ";
    //nuclearChannel += Nb + name;
    nuclearChannel += name;
  }
 
  //G4cout << "\n nuclear channel: " << nuclearChannel << G4endl;
  run->CountNuclearChannel(nuclearChannel, Q);

  G4double edep = step->GetTotalEnergyDeposit();
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* alpha = particleTable->FindParticle("alpha");
  G4int iabs = prePoint->GetTouchableHandle()->GetCopyNumber(0); 
//   for (int i = 0; i < 8; i++) 
//  {
//             if (particle == G4IonTable::GetIonTable()->GetIon(3, 7) && iabs == i) {
//                 G4double kine = step->GetPreStepPoint()->GetKineticEnergy();
//                 //fEventAction->AddKine(i, kine);
// // G4cout<<kine/MeV<<G4endl;
//                 G4double xLi = point.x();
//                 G4double xshiftedLi = xLi + 0.5*fDetector->GetAbsorSizeX();  
//                 analysisManager->FillH1(14, xshiftedLi, kine/MeV);
//             }
//             else if 
//             (particle == alpha && iabs == i) {
//                 G4double xAl = point.x();
//                 G4double xshiftedAl = xAl + 0.5*fDetector->GetAbsorSizeX();  
//                 G4double kine = step->GetPreStepPoint()->GetKineticEnergy();
//                 //fEventAction->AddKine(i, kine);
// // G4cout<<kine/MeV<<G4endl;
//                 //G4cout << " kinetic energy " << kine/MeV << " Xalpha " << xshiftedAl << G4endl;
//                 analysisManager->FillH1(15, xshiftedAl, kine/MeV);
//             }


  // if (edep <= 0.) return;
  if (edep > 0.){
  if (step->GetTrack()->GetDefinition()->GetPDGCharge() == 0.) point = P2;
  G4double x = point.x();
  G4double y = point.y();
  G4double z = point.z();
  G4double xshifted = x + 0.5*fDetector->GetAbsorSizeX();  
  G4double yshifted = y + 0.5*fDetector->GetAbsorSizeYZ();  
  G4double zshifted = z + 0.5*fDetector->GetAbsorSizeYZ();  
  
  analysisManager->FillH1(22, xshifted, edep);
  analysisManager->FillH1(23, yshifted, edep);
    
  //   if (xshifted > 0.){
  //     G4cout
  //     <<edep*2250* keV/um
  // //G4cout << edep*225* keV/um
  //    // << x/um
  //    << "\n xshifted" << x/um
  //     << G4endl;
  //   }

  analysisManager->FillH1(24, zshifted, edep);

  // G4double stepLength = 0;
  // stepLength = step->GetStepLength();
  
  // analysisManager->FillNtupleDColumn(0,xshifted/um);
  // analysisManager->FillNtupleDColumn(1,yshifted/um);
  // analysisManager->FillNtupleDColumn(2,zshifted/um);
   // analysisManager->FillNtupleDColumn(3,edep/stepLength/keV);
  // analysisManager->AddNtupleRow();
  
  //add
  // get volume of the current step
    G4LogicalVolume *volume
            = step->GetPreStepPoint()->GetTouchableHandle()
                    ->GetVolume()->GetLogicalVolume();

    G4String volname = volume->GetName();
  
    // check if we are in scoring volume

    // if (volname == "Target1") {
      
    //     G4StepPoint *point2 = step->GetPreStepPoint();
    //     G4TouchableHandle touch2 = point2->GetTouchableHandle();
    //     G4int copyno = touch2->GetCopyNumber();
 //end

 //total energy deposit in absorbers
 //
 if (iabs > 0) fEventAction->AddEdep(iabs, edep);
  fParticleFlag.clear();
  //G4RunManager::GetRunManager()->AbortEvent();  
}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
