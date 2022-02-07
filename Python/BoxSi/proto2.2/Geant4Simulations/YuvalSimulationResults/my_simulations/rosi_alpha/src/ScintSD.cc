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
/// \file ScintSD.cc
/// \brief Implementation of the ScintSD class

#include "ScintHit.hh"
#include "ScintSD.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4Gamma.hh"
#include "G4Proton.hh"
#include "G4Gamma.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ScintSD::ScintSD(G4String name) :
		G4VSensitiveDetector(name), fHitsCollection(nullptr), fHCID(-1)
{
	collectionName.insert("scintColl");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ScintSD::~ScintSD()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ScintSD::Initialize(G4HCofThisEvent *hce)
{
	fHitsCollection = new ScintHitsCollection(SensitiveDetectorName, collectionName[0]);
	if (fHCID < 0)
	{
		fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
	}
	hce->AddHitsCollection(fHCID, fHitsCollection);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool ScintSD::ProcessHits(G4Step *step, G4TouchableHistory*)
{
	/*
	 bool save = false;
	 auto edep = step->GetTotalEnergyDeposit();
	 if (edep == 0.)
	 return true;

	 auto touchable = step->GetPreStepPoint()->GetTouchable();
	 auto physical = touchable->GetVolume();
	 auto copyNo = physical->GetCopyNo();



	 hit->SetEdep(step->GetTotalEnergyDeposit());
	 auto preStepPoint = step->GetPreStepPoint();
	 if (step->GetTrack()->GetParentID() == 1 && step->GetDeltaEnergy() != 0)
	 hit->SetTime(preStepPoint->GetLocalTime());
	 else
	 hit->SetTime(-1.0 * s);
	 hit->SetID(copyNo);

	 */
	if (step->GetTrack()->GetParentID() != 1)
		return false;
	auto hit = new ScintHit();
	hit->SetEdep(
			step->GetPreStepPoint()->GetKineticEnergy()
					- step->GetPostStepPoint()->GetKineticEnergy());
	bool keep_pulse = false;
	auto particle_def = step->GetTrack()->GetParticleDefinition();

	if (particle_def == G4Gamma::Gamma())
		keep_pulse = true;
	auto seconderies = step->GetSecondaryInCurrentStep();
	//G4cout<<particle_def->GetParticleName()<<G4endl;
	for (auto it = seconderies->begin(); it != seconderies->end(); ++it)
	{
		//G4cout<<(*it)->GetParticleDefinition()->GetParticleName()<<G4endl;
		//if(particle_def==G4Neutron::Neutron())
		//G4cout<<particle_def->GetParticleName()<<G4endl;

		if ((*it)->GetParticleDefinition() == G4Proton::Proton()
				&& particle_def == G4Neutron::Neutron())
		{
			keep_pulse = true;
			break;
		}

	}

	G4int PDGEcode = particle_def->GetPDGEncoding();
	hit->SetPDGCode(PDGEcode);

	if (hit->GetEdep() / keV > 0.1 && keep_pulse)
	{
		hit->SetTrackID(step->GetTrack()->GetTrackID());
		hit->SetTime(step->GetPreStepPoint()->GetLocalTime());
		fHitsCollection->insert(hit);
	}

	return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
