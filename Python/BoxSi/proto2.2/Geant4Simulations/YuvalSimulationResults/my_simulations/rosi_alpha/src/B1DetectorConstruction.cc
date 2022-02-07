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
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"
#include "ScintSD.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4GenericMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction() :
		G4VUserDetectorConstruction(), fScoringVolume(0), fLogicalVolumes(0), fSolids(0), fScint1(
				0), fScint2(0)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{
	int no_of_logical = 4;
	fLogicalVolumes.resize(no_of_logical);
	fSolids.resize(no_of_logical);
	// Get nist material manager
	G4NistManager *nist = G4NistManager::Instance();

	G4bool checkOverlaps = true;

	//
	// World
	//
	G4double world_sizeXY = 5 * m;
	G4double world_sizeZ = 5 * m;
	G4Material *world_mat = nist->FindOrBuildMaterial("G4_AIR");

	G4Box *solidWorld = new G4Box("World", 0.5 * world_sizeXY, 0.5 * world_sizeXY,
			0.5 * world_sizeZ);     //its size

	G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld, world_mat, "World");    //its name
	G4VPhysicalVolume *physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0,
			false, 0, checkOverlaps); //overlaps checking

	fScoringVolume = logicWorld;

	G4Material *scint_mat = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
	auto scnit1Solid1 = new G4Box("scintLarge", 51 * mm / 2, 51 * mm / 2, 100.5 * mm / 2);
	fLogicalVolumes[0] = new G4LogicalVolume(scnit1Solid1, scint_mat, "scint1Logical");
	fLogicalVolumes[1] = new G4LogicalVolume(scnit1Solid1, scint_mat, "scint2Logical");
	auto position1 = G4ThreeVector((95.23 - 52 / 2) * mm, 0, 0);
	new G4PVPlacement(0, position1, fLogicalVolumes[0], "scint1Phys", logicWorld, false, 0,checkOverlaps);
	new G4PVPlacement(0, -position1, fLogicalVolumes[1], "scint2Phys", logicWorld, false, 1,checkOverlaps);

	auto scnit1Solid2 = new G4Box("scintLarge", 29.2 * mm / 2, 28.8 * mm / 2, 100.68 * mm / 2);
	fLogicalVolumes[2] = new G4LogicalVolume(scnit1Solid2, scint_mat, "scint3Logical");
	fLogicalVolumes[3] = new G4LogicalVolume(scnit1Solid2, scint_mat, "scint4Logical");
	auto position2 = G4ThreeVector(0, (95.23 - 52 / 2) * mm, 0);
	new G4PVPlacement(0, position2, fLogicalVolumes[2], "scint3Phys", logicWorld, false, 0,checkOverlaps);
	new G4PVPlacement(0, -position2, fLogicalVolumes[3], "scint4Phys", logicWorld, false, 1,checkOverlaps);

	return physWorld;
}
void B1DetectorConstruction::ConstructSDandField()
{

	auto sdManager = G4SDManager::GetSDMpointer();
	G4String SDname;
	ScintSD *scintSD;
	for (std::size_t i = 0; i < fLogicalVolumes.size(); i++)
	{
		scintSD=new ScintSD(SDname = "/Scint"+std::to_string(i+1));
		sdManager->AddNewDetector(scintSD);
		fLogicalVolumes[i]->SetSensitiveDetector(scintSD);
	}

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
