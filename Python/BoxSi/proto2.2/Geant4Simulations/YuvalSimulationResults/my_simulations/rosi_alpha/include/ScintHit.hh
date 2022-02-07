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
/// \file B5EmCalorimeterHit.hh
/// \brief Definition of the B5EmCalorimeterHit class

#ifndef ScintHit_h
#define ScintHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4AttDef;
class G4AttValue;

/// EM Calorimeter hit
///
/// It records:
/// - the cell ID
/// - the energy deposit 
/// - the cell logical volume, its position and rotation

class ScintHit : public G4VHit
{
  public:
    ScintHit();
    ScintHit(G4int cellID);
    ScintHit(const ScintHit &right);
    virtual ~ScintHit();

    const ScintHit& operator=(const ScintHit &right);
    G4bool operator==(const ScintHit &right) const;
    
    inline void *operator new(size_t);
    inline void operator delete(void *aHit);

    void SetEdep(G4double de) { fEdep = de; }
    void AddEdep(G4double de) { fEdep += de; }
    G4double GetEdep() const { return fEdep; }

    void SetTime(G4double t) { fTime = t; }
    G4double GetTime() const { return fTime; }

    void SetPDGCode(G4int c) { fPDGCode = c; }
    G4double GetPDGCode() const { return fPDGCode; }

    void SetTrackID(G4int id) { fTrackID = id; }
    G4double GetTrackID() const { return fTrackID; }

    
  private:
    G4int fID;
    G4double fEdep;
    G4double fTime;
    G4int fPDGCode;
    G4int fTrackID;
};

using ScintHitsCollection = G4THitsCollection<ScintHit>;

extern G4ThreadLocal G4Allocator<ScintHit>* ScintHitAllocator;

inline void* ScintHit::operator new(size_t)
{
  if (!ScintHitAllocator) {
       ScintHitAllocator = new G4Allocator<ScintHit>;
  }
  return (void*)ScintHitAllocator->MallocSingle();
}

inline void ScintHit::operator delete(void* aHit)
{
  ScintHitAllocator->FreeSingle((ScintHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
