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
/// \file HistoManager.cc
/// \brief Implementation of the HistoManager class
//
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "DetectorConstruction.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
  : fFileName("Hadr07")
{
  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");
  analysisManager->SetFileName(fFileName);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);

  // Define histograms start values
  const G4int kMaxHisto = 30;
  const G4String id[] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
                         "10","11","12","13","14","15","16","17","18","19",
                         "20","21","23","24","25","26","27","28","29","30"};

  const G4String title[] = 
                { "dummy",                                        //0
                  "total Energy deposited in absorber 1 per event",         //1
                  "total Energy deposited in absorber 2 per event",         //2
                  "total Energy deposited in absorber 3 per event",         //3
                  "total Energy deposited in absorber 4 per event",         //4
                  "total Energy deposited in absorber 5 per event",         //5
                  "total Energy deposited in absorber 6 per event",         //6
                  "total Energy deposited in absorber 7 per event",         //7
                  "total Energy deposited in absorber 8 per event",         //8
                  "total Energy deposited in absorber 9 per event",         //9
                  "total Energy deposited in absorber 10 per event",         //10
                  "total Energy deposited in absorber 11 per event",         //11
                  "total Energy deposited in absorber 12 per event",         //12
                  "total Energy deposited in absorber 13 per event",         //13
                  "total Energy deposited in absorber 14 per event",         //14
                  "total Energy deposited in absorber 15 per event",         //15
                  "total Energy deposited in absorber 16 per event",         //16
                  "total Energy deposited in absorber 17 per event",         //17
                  "total Energy deposited in absorber 18 per event",         //18
                  "total Energy deposited in absorber 19 per event",         //19
                  "total Energy deposited in absorber 20 per event",         //20
                  "total Energy deposited in absorbers per event",         //21
                  "Edep (keV/um) along x axix in absorbers" ,                //22
                  "Edep (keV/um) along y axix in absorbers" ,                //23
                  "Edep (keV/um) along z axix in absorbers" ,                //24
                  // "x y position" ,                //25
                  // "kinetic energy of lithium" ,                //26
                  // "kinetic energy of alpha" ,                //27
                 };

  // Default values (to be reset via /analysis/h1/set command)
  G4int nbins = 100;
  G4double vmin = 0.;
  G4double vmax = 100.;

  // Create all histograms as inactivated 
  // as we have not yet set nbins, vmin, vmax
  for (G4int k=0; k<kMaxHisto; k++) {
    G4int ih = analysisManager->CreateH1(id[k], title[k], nbins, vmin, vmax);
    analysisManager->SetH1Activation(ih, false);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
