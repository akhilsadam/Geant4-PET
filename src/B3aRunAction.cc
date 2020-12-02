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
//
/// \file B3aRunAction.cc
/// \brief Implementation of the B3aRunAction class

#include "B3aRunAction.hh"
#include "B3PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "B3aRun.hh"//"G4Run.hh"
#include "G4AccumulableManager.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "B3aHistoManager.hh"
#include "B3DetectorConstruction.hh"
#include "B3SteppingAction.hh"
#include <mutex>
//#include "B3aHistoManager.cc" ///  NEED TO FIX -----
//#include "B3SteppingAction.cc" ///  NEED TO FIX -----
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::mutex foo;
std::mutex barL;
G4int B3aRunAction::id = 0;

B3aRunAction::B3aRunAction(B3DetectorConstruction* patient)
 : G4UserRunAction(),
   fGoodEvents(0),
   fSumDose(0.),
   fpatient(patient),
   fHistoManager(0)
{  
  //add new units for dose
  // 
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;  
  const G4double picogray  = 1.e-12*gray;
   
  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray);

  // Register accumulable to the accumulable manager
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->RegisterAccumulable(fGoodEvents);
  accumulableManager->RegisterAccumulable(fSumDose); 

  fHistoManager = new HistoManager(fpatient);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3aRunAction::~B3aRunAction()
{
   delete fHistoManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3aRunAction::BeginOfRunAction(const G4Run* run)
{ 
  G4cout << "### Run " << run->GetRunID() << " start." << G4endl;
  
  // reset accumulables to their initial values
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();

  //histograms
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->IsActive() ) {
    analysisManager->OpenFile();
  }  

  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3aRunAction::EndOfRunAction(const G4Run* run)
{  
  std::lock(foo,barL);

  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) 
	{
		foo.unlock();
		barL.unlock();	
		return;	
	}

  // Merge accumulables 
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
//integrate histograms sanitycheck

  if (analysisManager->IsActive() && (!IsMaster())) {
	
	//integrate the hist for Econs. sanity check
	//Primary
	const std::vector<double> prmE = (const std::vector<double>) analysisManager->GetH1(9)->bins_sum_w();
	//Secondary
	const std::vector<double> secE = (const std::vector<double>) analysisManager->GetH1(10)->bins_sum_w();

	//G4cout << "stepmax = " << stepMax << G4endl;
	

	G4double totalsent = 0;
	double lastEnergy = 100;
	G4double secondE = 0;
	G4int stepMax = fHistoManager->stepMaxV;
	bool proceed = true;	

	for(int i = 0; i < (stepMax); i++)
	{
		/*if((i>0) && (prmE[i] > prmE[(i-1)]))
		{
			//secondE = 0; //reset
			G4cout << "RESET -- PRMI: " << prmE[i] << " PRMI-1: " << prmE[(i-1)] << G4endl;
		}*/
		secondE += secE[i];
		/*if(prmE[i] == 0)
		{	
			if((i < (stepMax+1))&&( prmE[i+1] > lastEnergy ))
			{
				secondE = 0;
				G4cout << "REset" << G4endl;
			}
			else
			{
				///G4cout << "HIT zero -- BREAK" << (secE[i]) << G4endl;
				G4cout << "HIT prmE = " << (prmE[i]) << G4endl;
				G4cout << "HIT secE = " << (secondE) << G4endl;///
			}
		}*/

		double diff = (prmE[i] - lastEnergy);
		if((i>0) && (diff > 0.0000001))
		{
			secondE = secE[i];
			G4cout << "||\\-----RESET RESET -- secondE = " << secondE << " diff = "<< diff << G4endl;
			if(prmE[i] < 19)
			{
				proceed = false;
			}
			else
			{
				proceed = true;
			}
			
		}
		//if(prmE[i]+secondE > 40 ){secondE = 0;}

		lastEnergy = prmE[i];	

		if (proceed)//(prmE[i] > 0)
		{	
			totalsent = (prmE[i] + secondE);
			/*if(abs(totalsent - 40.272013) > 1.)
			{
				G4cout << "totalsent = " << (totalsent) << G4endl;
				G4cout << "prmE = " << (prmE[i]) << G4endl;
				G4cout << "secE = " << (secondE) << G4endl;
			}*/
		

			analysisManager->FillH1(11,id,totalsent);
			id++;
			G4cout << "totalsent = " << (totalsent) << " prmE = " << (prmE[i]) << " secondE = " << (secondE) << " id = " << (id-1) << G4endl;

		}
		

		if((diff == 0) && (secE[i] ==0))
		{
			secondE =0;
		}
	

		/*if((secE[i] != 0) && (prmE[i] == 0))
		{
			secondE = 0;
		}*/
		//G4cout << "-- PRMI: " << prmE[i] << " PRMI-1: " << prmE[(i-2)] << G4endl;
	}
	G4cout << "EVENT ACTION BREAK" << G4endl;


  } 


  // save histograms
  if ( analysisManager->IsActive() ) {  
    analysisManager->Write();
    analysisManager->CloseFile();
  }    


  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const B3PrimaryGeneratorAction* generatorAction
    = static_cast<const B3PrimaryGeneratorAction*>(
        G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String partName;
  if (generatorAction) 
  {
    G4ParticleDefinition* particle 
      = generatorAction->GetParticleGun()->GetParticleDefinition();
    partName = particle->GetParticleName();
  }  
 
          
  // Print results
  //
  if (IsMaster())
  {
    G4cout
     << G4endl
     << "--------------------End of Global Run-----------------------"
     << G4endl
     << "  The run was " << nofEvents << " events ";
     B3SteppingAction::id = 0;
  }
  else
  {
    G4cout
     << G4endl
     << "--------------------End of Local Run------------------------"
     << G4endl
     << "  The run was " << nofEvents << " "<< partName;
  }      
  G4cout
     << "; Nb of 'good' e+ annihilations: " << fGoodEvents.GetValue()  << G4endl
     << " Total dose in patient : " << G4BestUnit(fSumDose.GetValue(),"Dose") 
     << G4endl 
     << "------------------------------------------------------------" << G4endl 
     << G4endl;

foo.unlock();
barL.unlock();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
