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
/// \file electromagnetic/TestEm11/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "B3aHistoManager.hh"
#include "B3DetectorConstruction.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager(B3DetectorConstruction* patient)
  : fFileName("B3A"),fpatient(patient)
{
  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFileName(fFileName);
  analysisManager->SetVerboseLevel(0);
  analysisManager->SetActivation(true);

  // Define histograms start values
  const G4int kMaxHisto = 2;
  const G4String id[] = { "0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
                         "10","11","12","13","14","15","16","17","18","19",
                         "20","21","22"};

  /*const G4String title[] = 
                { "dummy",                                        //0
                  "Edep (MeV/mm) along absorbers",                //1
                  "total Energy deposited in absorbers",          //2
                  "true track length of the primary particle",    //3
                  "true step size of the primary particle",       //4
                  "projected range of the primary particle",      //5
                  "true track length of charged secondaries",     //6
                  "true step size of charged secondaries",        //7
                  "Edep (MeV.cm2/g) along x/r0",                  //8
                  "dummy",                                        //9
                  "dummy"                                         //10
                 };
*/
   const G4String title[] = 
                { "E-Deposition (MeV/mm)",
		  "E-Deposition 2D (MeV/mm)",
		  "E-Deposition 3D (MeV/mm)",
		  "Secondary List"};

   const std::string second[] = { "positron","electron","photons","gammas","proton","alpha","Li6","Be7","C11","C12","N15","O15","O16"};
   const int secondSize = sizeof(second)/sizeof(second[0]);
  // Default values (to be reset via /analysis/h1/set command)
  G4double* worldsizeP;
  G4double worldsize = fpatient->GetWorldSizeXY();
  worldsizeP = new G4double(worldsize);

  G4int nbinsx = 100;
  G4double xmin = 0.;
  G4double xmax = 237.958;
  G4int nbinsy = 100;
  G4double ymin = 0.;
  G4double ymax = 237.958;
  G4int nbinsz = 100;
  G4double zmin = 0.;
  G4double zmax = 324;

  G4double sbins = secondSize;

    G4int ih = analysisManager->CreateH1(id[0], title[0], nbinsx, xmin, xmax);
    analysisManager->SetH1Activation(ih, true);
G4cout << "### ####### " << worldsize << " worldsize." << G4endl;
    G4int ih2 = analysisManager->CreateH2(id[1], title[1], nbinsx, xmin, xmax, nbinsy, ymin, ymax);
    analysisManager->SetH2Activation(ih2, true);
    G4int ih3 = analysisManager->CreateH3(id[2], title[2], nbinsx, xmin, xmax, nbinsy, ymin, ymax, nbinsz, zmin, zmax);
    analysisManager->SetH2Activation(ih3, true);

    /*analysisManager->CreateNtuple(id[3], title[3]);
	for(int i = 0; i < secondSize; i++)
	{
		analysisManager->CreateNtupleDColumn(second[i]);
	}
	analysisManager->FinishNtuple();*/
    

    G4int ih4 = analysisManager->CreateH1(id[3], title[3], sbins, 0, sbins);
    analysisManager->SetH1Activation(ih4, true);
/*
  // Create all histograms as inactivated 
  // as we have not yet set nbins, vmin, vmax
  for (G4int k=0; k<=kMaxHisto; k++) {

  }
  
 G4String title2;
  for (G4int k=1; k<kMaxAbsor; k++) {
    title2 = "Edep in absorber " + id[k];
    G4int ih 
      = analysisManager->CreateH1(id[kMaxHisto+k], title2, nbins, vmin, vmax);
    analysisManager->SetH1Activation(ih, false);
  }*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
