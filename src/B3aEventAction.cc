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
/// \file B3aEventAction.cc
/// \brief Implementation of the B3aEventAction class

#include "B3aEventAction.hh"
#include "B3aRunAction.hh"
#include "B3Analysis.hh"
#include "B3aHistoManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"

#include "G4EventManager.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


G4int fCollID_cryst_p;
G4int fCollID_cryst_ep;
G4int fCollID_cryst_en;
G4int fCollID_cryst_y;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3aEventAction::B3aEventAction(B3aRunAction* runAction)
 : G4UserEventAction(), 
   fRunAction(runAction),
   fCollID_cryst(-1),
/* do we need this ??
   fCollID_cryst_p(-1),
   fCollID_cryst_ep(-1),
   fCollID_cryst_en(-1),
*/
   fCollID_patient(-1)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3aEventAction::~B3aEventAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3aEventAction::BeginOfEventAction(const G4Event* /*evt*/)
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3aEventAction::EndOfEventAction(const G4Event* evt )
{
   //Hits collections
  //  
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  if(!HCE) return;
               
   // Get hits collections IDs
  if (fCollID_cryst < 0) {
    G4SDManager* SDMan = G4SDManager::GetSDMpointer();  
    fCollID_cryst   = SDMan->GetCollectionID("crystal/edep");

    fCollID_cryst_p    = SDMan->GetCollectionID("crystal/edep_p");
    fCollID_cryst_ep   = SDMan->GetCollectionID("crystal/edep_ep");
    fCollID_cryst_en   = SDMan->GetCollectionID("crystal/edep_en");
    fCollID_cryst_y   = SDMan->GetCollectionID("crystal/edep_y");

    fCollID_patient = SDMan->GetCollectionID("patient/dose");    
  }

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  //Energy in crystals : identify 'good events'
  //
  const G4double eThreshold = 500*keV;
  G4int nbOfFired = 0;
   
  G4THitsMap<G4double>* evtMap = 
                     (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_cryst));
               
  std::map<G4int,G4double*>::iterator itr;
  for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++) {
    //G4int copyNb  = (itr->first);
    G4double edep = *(itr->second);
    if (edep > eThreshold) nbOfFired++;
    //G4cout << "\n  cryst" << copyNb << ": " << edep/keV << " keV ";
    *(itr->second) = edep/MeV; // MODIFIED TO CONVERT UNITS
	analysisManager->FillH1(2, *(itr->second));
  }  
  if (nbOfFired == 2) fRunAction->CountEvent();
 // MODIFIED ----------------------------------------------------------------------------------------
  evtMap = (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_cryst_p));
  for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++) {
    G4double edep_p = *(itr->second);
    *(itr->second) = edep_p/MeV; // MODIFIED TO CONVERT UNITS
	analysisManager->FillH1(3, *(itr->second));
  }
  // MODIFIED ----------------------------------------------------------------------------------------
  evtMap = (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_cryst_ep));
  for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++) {
    G4double edep_ep = *(itr->second);
    *(itr->second) = edep_ep/MeV; // MODIFIED TO CONVERT UNITS
	analysisManager->FillH1(4, *(itr->second));
  }  
  // MODIFIED ----------------------------------------------------------------------------------------
  evtMap = (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_cryst_en));
  for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++) {
    G4double edep_en = *(itr->second);
    *(itr->second) = edep_en/MeV; // MODIFIED TO CONVERT UNITS
	analysisManager->FillH1(5, *(itr->second));
  }
  // MODIFIED ----------------------------------------------------------------------------------------
  evtMap = (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_cryst_y));
  for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++) {
    G4double edep_y = *(itr->second);
    *(itr->second) = edep_y/MeV; // MODIFIED TO CONVERT UNITS
	analysisManager->FillH1(6, *(itr->second));
  }
  // ------------------------------------------------------------------------------------------------

  //Dose deposit in patient
  //
  G4double dose = 0.;
     
  evtMap = (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_patient));
               
  for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++) {
    ///G4int copyNb  = (itr->first);
    dose = *(itr->second);
    *(itr->second) = dose/gray; //MODIFIED TO CONVERT UNITS
	analysisManager->FillH1(7, *(itr->second));
  }
  if (dose > 0.) fRunAction->SumDose(dose);
// ------------------------------------------------------------------------------------------------
  //
  //analysisManager->FillH1(1, fTotalEnergyDeposit/MeV);


//integrate histograms sanitycheck
  if (analysisManager->IsActive()) {

	//integrate the hist for Econs. sanity check
	//Primary
	const std::vector<double> prmE = (const std::vector<double>) analysisManager->GetH1(9)->bins_sum_w();
	//Secondary
	const std::vector<double> secE = (const std::vector<double>) analysisManager->GetH1(10)->bins_sum_w();

	//G4cout << "stepmax = " << stepMax << G4endl;
	G4int stepMax = fRunAction->fstepMax;
	

	G4double totalsent = 0;
	G4double lastEnergy =0;
	for(int i = 0; i < (stepMax); i++)
	{
		/*if((i>0) && (prmE[i] > prmE[(i-1)]))
		{
			//secondE = 0; //reset
			G4cout << "RESET -- PRMI: " << prmE[i] << " PRMI-1: " << prmE[(i-1)] << G4endl;
		}*/
		secondE += secE[i];
		if(prmE[i] == 0)
		{	
			if((i < (stepMax+1))&&( prmE[i+1] > lastEnergy ))
			{
				secondE = 0;
				G4cout << "REset" << G4endl;
			}
			else
			{
				/*G4cout << "HIT zero -- BREAK" << (secE[i]) << G4endl;
				G4cout << "HIT prmE = " << (prmE[i]) << G4endl;
				G4cout << "HIT secE = " << (secondE) << G4endl;*/
			}
		}



		if (prmE[i] > 0)
		{	
			totalsent = (prmE[i] + secondE);
			if(abs(totalsent - 40.272013) > 100.)
			{
				G4cout << "totalsent = " << (totalsent) << G4endl;
				G4cout << "prmE = " << (prmE[i]) << G4endl;
				G4cout << "secE = " << (secondE) << G4endl;
			}
			else
			{

				analysisManager->FillH1(11,i,totalsent);
			}
			lastEnergy = prmE[i];
			
		}
	

		/*if((secE[i] != 0) && (prmE[i] == 0))
		{
			secondE = 0;
		}*/
		//G4cout << "-- PRMI: " << prmE[i] << " PRMI-1: " << prmE[(i-2)] << G4endl;
	}
	G4cout << "EVENT ACTION BREAK" << G4endl;


  } 
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
