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
/// \file B1EventAction.cc
/// \brief Implementation of the B1EventAction class

#include <RunAction.hh>
#include "B1EventAction.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
using std::array;
using std::vector;

#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "B1Analysis.hh"
#include <numeric>
#include "ScintHit.hh"

namespace
{

// Utility function which finds a hit collection with the given Id
// and print warnings if not found
G4VHitsCollection* GetHC(const G4Event *event, G4int collId)
{
	auto hce = event->GetHCofThisEvent();
	if (!hce)
	{
		G4ExceptionDescription msg;
		msg << "No hits collection of this event found." << G4endl;
		G4Exception("B2EventAction::EndOfEventAction()", "B2Code001", JustWarning, msg);
		return nullptr;
	}

	auto hc = hce->GetHC(collId);
	if (!hc)
	{
		G4ExceptionDescription msg;
		msg << "Hits collection " << collId << " of this event not found." << G4endl;
		G4Exception("B1EventAction::EndOfEventAction()", "B1Code001", JustWarning, msg);
	}
	return hc;
}

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::B1EventAction() :
		G4UserEventAction(), fEdep(0.), fNeutron
		{
		{ std::vector<G4double>(0, 0.0), std::vector<G4double>(0, 0.0), std::vector<G4double>(0,
				0.0) } }, fGamma
		{
		{ std::vector<G4double>(0, 0.0), std::vector<G4double>(0, 0.0), std::vector<G4double>(0,
				0.0) } }, fScintHCID
		{
		{ -1, -1, -1, -1 } }, fNneutrons(0), fNgammas(0),fIsFission(true)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::~B1EventAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::BeginOfEventAction(const G4Event*)
{
	fEdep = 0.;
	fNneutrons = fNgammas = 0;
	fIsFission=true;
	if (fScintHCID[0] == -1)
	{
		auto sdManager = G4SDManager::GetSDMpointer();

		// hits collections names
		array<G4String, kDim> scintHCName =
		{
		{ "Scint1/scintColl", "Scint2/scintColl", "Scint3/scintColl", "Scint4/scintColl" } };

		for (G4int iDet = 0; iDet < kDim; ++iDet)
		{
			// hit collections IDs
			fScintHCID[iDet] = sdManager->GetCollectionID(scintHCName[iDet]);
		}
	}
	fNeutron =
	{ std::vector<G4double>(0, 0.0), std::vector<G4double>(0, 0.0), std::vector<G4double>(0, 0.0) };
	fGamma =
	{ std::vector<G4double>(0, 0.0), std::vector<G4double>(0, 0.0), std::vector<G4double>(0, 0.0) };
}
struct PTrack
{
	std::vector<G4double> times;
	std::vector<G4double> energies;
	G4int PDGcode = -1;
	G4int detector_id = -1;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void B1EventAction::EndOfEventAction(const G4Event *event)
{
	G4int PDGEcode;
	int trackID;
	G4bool add_row = false;

	// Old Analysis
	for (G4int iDet = 0; iDet < kDim; ++iDet)
	{
		auto hc = GetHC(event, fScintHCID[iDet]);
		if (!hc)
			return; //CHECK THIS

		for (unsigned int i = 0; i < hc->GetSize(); ++i)
		{
			auto hit = static_cast<ScintHit*>(hc->GetHit(i));
			PDGEcode = hit->GetPDGCode();
			if (hit->GetPDGCode() == 2112)
			{
				fNeutron[0].push_back(hit->GetTime() / s);
				fNeutron[1].push_back(hit->GetEdep() / MeV);
				fNeutron[2].push_back(iDet);
				add_row = true;
			}
			if (hit->GetPDGCode() == 22)
			{
				fGamma[0].push_back(hit->GetTime() / s);
				fGamma[1].push_back(hit->GetEdep() / MeV);
				fGamma[2].push_back(iDet);
				add_row = true;
			}
		}
	}
	auto analysisManager = G4AnalysisManager::Instance();
	if (add_row)
	{

		analysisManager->AddNtupleRow();
	}
	if (fIsFission)
	{
		analysisManager->FillNtupleDColumn(1, 0, fNneutrons);
		analysisManager->AddNtupleRow(1);
		analysisManager->FillNtupleDColumn(2, 0, fNgammas);
		analysisManager->AddNtupleRow(2);
	}
	/*
	 //Track analysis
	 std::map<std::pair<int,int>, PTrack> tracks;
	 tracks.clear();

	 for (int iDet = 0; iDet < kDim; ++iDet)
	 {
	 auto hc = GetHC(event, fScintHCID[iDet]);

	 for (unsigned int i = 0; i < hc->GetSize(); ++i)
	 {
	 auto hit = static_cast<ScintHit*>(hc->GetHit(i));
	 trackID = hit->GetTrackID();
	 if (!tracks.count(std::make_pair(iDet,trackID)))
	 {
	 tracks[std::make_pair(iDet,trackID)] = PTrack();

	 }
	 tracks[std::make_pair(iDet,trackID)].energies.push_back(hit->GetEdep() / MeV);
	 tracks[std::make_pair(iDet,trackID)].times.push_back(hit->GetTime() / s);
	 tracks[std::make_pair(iDet,trackID)].PDGcode = hit->GetPDGCode();;
	 tracks[std::make_pair(iDet,trackID)].detector_id = iDet+5;
	 }
	 }
	 std::map<std::pair<int,int>, PTrack>::iterator it;
	 G4double energy,time;
	 int time_counter = 0;
	 for (it = tracks.begin(); it != tracks.end(); it++)
	 {
	 auto temp_det_track = it->second;
	 time = energy = 0;
	 time_counter=0;

	 for (auto &n : temp_det_track.energies)
	 energy += n;

	 for (unsigned int i = 0; i < temp_det_track.times.size(); i++) //FIX average only for significant deposition to avoid thermalization.
	 {
	 if (temp_det_track.energies[i]>0.001)
	 {
	 time += temp_det_track.times[i];
	 time_counter++;
	 }
	 }
	 time = time / time_counter;
	 if (temp_det_track.PDGcode == 2112)
	 {
	 fNeutron[0].push_back(time);
	 fNeutron[1].push_back(energy);
	 fNeutron[2].push_back(temp_det_track.detector_id);
	 add_row = true;
	 }
	 else if (temp_det_track.PDGcode == 22)
	 {
	 fGamma[0].push_back(time);
	 fGamma[1].push_back(energy);
	 fGamma[2].push_back(temp_det_track.detector_id);
	 add_row = true;
	 }

	 }
	 if (add_row)
	 {
	 auto analysisManager = G4AnalysisManager::Instance();
	 analysisManager->AddNtupleRow();
	 }
	 */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
