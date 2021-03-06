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
/// \file B5CsISD.cc
/// \brief Implementation of the B5CsISD class

#include "B5CsISD.hh"
#include "B5EmCalorimeterHit.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
extern bool gSaveStepLevel;
extern int gNumberOfScintillators;

B5CsISD::B5CsISD(G4String name, G4int xid, G4int yid)
  : G4VSensitiveDetector(name), fNameSD(name), fHitsCollection(nullptr), fHCID(-1)
{
  fEdep = 0;
  fxid = xid;
  fyid = yid;
  collectionName.insert("CsIHitCollection"); 
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5CsISD::~B5CsISD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5CsISD::Initialize(G4HCofThisEvent* hce)
{
  fHitsCollection = new B5CsIHitsCollection(fNameSD,collectionName[0]);
  if (fHCID<0) {
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
  }
  hce->AddHitsCollection(fHCID,fHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool B5CsISD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  auto edep = step->GetTotalEnergyDeposit();
  if (edep != 0){
    fEdep += edep;
  }
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5CsISD::EndOfEvent(G4HCofThisEvent* hce){
  if (fEdep == 0) return;
  fHitsCollection->insert(new B5CsIHit(fHCID));
  auto hit = (B5CsIHit*) ((hce -> GetHC(fHCID)) -> GetHit(0));
  
  hit -> SetEdep(fEdep);
  hit -> SetXID(fxid);
  hit -> SetYID(fyid);
  
  fEdep = 0;
}
