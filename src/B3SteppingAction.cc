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
/// \file electromagnetic/TestEm4/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "B3SteppingAction.hh"
#include "B3aEventAction.hh"
#include "G4SteppingManager.hh"
#include "G4RunManager.hh"
#include "Randomize.hh"

#include "B3DetectorConstruction.hh"
#include "B3aHistoManager.hh"

#include "G4Step.hh"
#include "G4Event.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


std::mutex foo2;
std::mutex barL2;
G4int B3SteppingAction::id = 0;
G4double lastEnergy=0;
//G4double seconds = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3SteppingAction::B3SteppingAction(B3aEventAction* EvAct, B3DetectorConstruction* patient)
:G4UserSteppingAction(),fEventAction(EvAct),fpatient(patient)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3SteppingAction::~B3SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3SteppingAction::UserSteppingAction(const G4Step* step)
{
std::lock(foo2,barL2);

G4double edep = step->GetTotalEnergyDeposit();

//if (edep <= 0.) return;

 G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();     

 //longitudinal profile of deposited energy
 //randomize point of energy deposotion
 //

 G4StepPoint* prePoint  = step->GetPreStepPoint();
 G4StepPoint* postPoint = step->GetPostStepPoint(); 
 G4ThreeVector P1 = prePoint ->GetPosition();
 G4ThreeVector P2 = postPoint->GetPosition();
 G4ThreeVector point = P1 + G4UniformRand()*(P2 - P1);
 if (step->GetTrack()->GetDefinition()->GetPDGCharge() == 0.) point = P2;
 G4double x = point.x();
 G4double y = point.y();
 G4double z = point.y();
 G4double xshifted = x + 0.5*(fpatient->GetWorldSizeXY());
 G4double yshifted = y + 0.5*(fpatient->GetWorldSizeXY());
 G4double zshifted = z + 0.5*(fpatient->GetWorldSizeZ());

/*G4cout
     << G4endl
     << "--------------------"<< fpatient->GetWorldSizeXY() <<"-----------------------"
     << G4endl;*/
 //G4double zshifted = x + 0.5*patient->GetAbsorSizeX();

 analysisManager->FillH1(0, xshifted, edep);//x,edep
 analysisManager->FillH2(0, xshifted, yshifted, edep);
 analysisManager->FillH3(0, xshifted, yshifted, zshifted, edep);

//implemented later: log secondaries only when they have been created
const std::vector< const G4Track* >* secondaries = step->GetSecondaryInCurrentStep();  
const int size = secondaries->size();
//std::vector< G4double> times = std::vector<G4double>();
//std::vector< const G4DynamicParticle*>* secondparticles = new std::vector< const G4DynamicParticle* >(size);

 //int nOfe = 0;
const std::string title[] = { "e+","e-","gamma","gamma","proton","alpha","Li6","Be7","C11","C12","N15","O15","O16"};
const int titleSize = sizeof(title)/sizeof(title[0]);

G4double timeL = prePoint->GetGlobalTime();

//implemented here: log secondaries only when they have been created
//log the primary only once!
//std::vector<G4double>::iterator it = find(times.begin(),times.end(),(timeL/ns));
// note make sure primary is proton!!
std::string prim = step->GetTrack()->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();


G4double Eprim = (step->GetTrack()->GetKineticEnergy()/MeV);
G4double Esec = 0;

if(prim.compare("proton")==0)
{
//G4cout << "//--Time (step): " << (id) << " step"<< G4endl;
	if(Eprim <=0)
		Eprim = 0;
	//ADD ENERGY TO TOTAL LIST
	G4cout << "-- Primary: " << step->GetTrack()->GetDynamicParticle()->GetParticleDefinition()->GetParticleName()<<" || Ep: "<<(Eprim) << "ID " << id << G4endl;
	G4cout << " -- Edep: " <<(edep/MeV)<<G4endl;
	Esec += (edep/MeV);	
	//times.push_back((timeL/ns));	
	//log the secondaries:
}
 for(int n = 0; n < size; n++)
{
	G4Track* ptrG4Track = (G4Track*)(*secondaries)[n];
	const G4DynamicParticle* pp = (const G4DynamicParticle*) ptrG4Track->GetDynamicParticle();
 	const G4ParticleDefinition* pd = (const G4ParticleDefinition*) pp->GetParticleDefinition();
	const std::string pdVar = pd->GetParticleName();

	G4cout << "-- Secondary: " << pdVar << " || Es: " << (pp->GetKineticEnergy()/MeV) << G4endl;
	const G4VProcess* prePointStep = prePoint->GetProcessDefinedStep();
	if(prePointStep)
	{	
		G4cout << "-- - Secondary process name: " << (prePoint->GetProcessDefinedStep()->GetProcessName()) << G4endl;	
	}
	Esec += (pp->GetKineticEnergy()/MeV);

//--------------------------------------///-------------------------//

	//analysisManager->FillH1(8, (timeL/ns), deltaE);
	if(pdVar.compare("e+")==0){
	G4cout << pdVar <<" " << (pp->GetKineticEnergy()/MeV) << G4endl;}
	bool filled = false;
	for(int i = 0; i < titleSize; i++)
	{
           //G4cout << pdVar << " " << title[i] << " " << titleSize << G4endl;
		if(pdVar.compare(title[i])==0)
		{
			if(i==2)
			{
				if(pp->GetTotalEnergy() > 500*keV)
				{analysisManager->FillH1(1, (0.5+3), 1);
				}
				else
				{analysisManager->FillH1(1, (0.5+2), 1);
				}
				//G4cout << "true " << pdVar << " " << title[i] << G4endl;
				filled = true;
				//add to histogram!!
				
				//analysisManager->FillNtupleDColumn(i,1);
				//analysisManager->AddNtupleRow();
				break;

			}
			filled = true;
			//add to histogram!!
			analysisManager->FillH1(1, (0.5+i), 1);
			//analysisManager->FillNtupleDColumn(i,1);
			//analysisManager->AddNtupleRow();
			break;

		}
	}
		
	if(!filled)
	{
	/*G4cout
		<< G4endl
		<< "------//\\---------- broke if statement "<< pdVar <<" -///////////\\\\\\\\\\"
		<< G4endl;*/
	}

//RECURSIVE SECONDARY ENERGIES::::



	//(*secondparticles)[n] = pp;

/*G4cout
     << G4endl
     << "-------------------- "<< nOfe <<" ----"<< n <<"-------------------"
     << "-------------------- "<< size <<" -----------------------"
     << G4endl;*/

}
if(prim.compare("proton")==0)
{
	//G4cout << "/|\\--- Total (should equal the previous primary): " << (Esec + Eprim) <<G4endl;
	if((Esec + Eprim)>80)
	{
		G4cout << "/|\\***********\\|/--- Error: E=" << ((Esec + Eprim)) << " MeV T=" << (timeL/ns) << " ns XYZ=" << xshifted << " " << yshifted  << " " << zshifted <<G4endl;
		G4cout << "/|\\***********\\|/---" <<G4endl;
		G4cout << "/|\\\\|/---" <<G4endl;
		G4cout << "/|\\\\|/---" <<G4endl;
		G4cout << "/|\\\\|/---" <<G4endl;
		G4cout << "/|\\***********\\|/---" <<G4endl;
		G4cout << "/|\\***********\\|/---" <<G4endl;
		G4cout << "/|\\***********\\|/---" <<G4endl;
		G4cout << "/|\\***********\\|/---" <<G4endl;
		G4cout << "/|\\***********\\|/---" <<G4endl;
	}
	else
	{	
		analysisManager->FillH1(8, (timeL/ns), (Eprim));
		analysisManager->FillH1(9, (id), (Eprim));
		analysisManager->FillH1(10, (id), (Esec));
		analysisManager->FillH1(12, (id), (Eprim-lastEnergy));
		analysisManager->FillH1(13, (id), (Eprim-lastEnergy+Esec));
		/*seconds+=Esec;
		if(Eprim > 30.272)//lastEnergy)
		{
			seconds = 0;//Esec;
			G4cout << "/|\\\\|/--- esec: " << Esec << "LastE: " << lastEnergy << " Eprim: " << Eprim <<G4endl;
		}
		analysisManager->FillH1(11, (id), (Eprim + seconds));
		lastEnergy = Eprim;*/

		id += 1;
		lastEnergy = Eprim;
	}
}
  
 //example of saving random number seed of this event, under condition
 //// if (condition) G4RunManager::GetRunManager()->rndmSaveThisEvent();  
foo2.unlock();
barL2.unlock();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

