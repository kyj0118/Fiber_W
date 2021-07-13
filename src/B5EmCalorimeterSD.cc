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
/// \file B5EmCalorimeterSD.cc
/// \brief Implementation of the B5EmCalorimeterSD class

#include "B5EmCalorimeterSD.hh"
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

B5EmCalorimeterSD::B5EmCalorimeterSD(G4String name, G4int layerNumber, G4int scintNumber)
  : G4VSensitiveDetector(name), fNameSD(name), fLayerId(layerNumber), fSegmentId(scintNumber), fHitsCollection(nullptr), fHCID(-1)
{
  fEdep.clear();
  fEweightedx.clear();
  fEweightedy.clear();
  fEweightedz.clear();
  fEweightedt.clear();

  fEdep.resize(gNumberOfScintillators,0);
  fEweightedx.resize(gNumberOfScintillators,0);
  fEweightedy.resize(gNumberOfScintillators,0);
  fEweightedz.resize(gNumberOfScintillators,0);
  fEweightedt.resize(gNumberOfScintillators,0);
  
  collectionName.insert("EMCalHitCollection"); 
  
  fStepEdep.clear();
  
  fPreStepx.clear();
  fPreStepy.clear();
  fPreStepz.clear();
  fPreStept.clear();
  
  fPostStepx.clear();
  fPostStepy.clear();
  fPostStepz.clear();
  fPostStept.clear();

  fParticlePx.clear();
  fParticlePy.clear();
  fParticlePz.clear();
  fParticleTrackID.clear();
  fParticleParentID.clear();
  fParticleCharge.clear();
  fParticleMass.clear();
  fParticlePDGID.clear();

  fStepEdep.resize(gNumberOfScintillators);
  
  fPreStepx.resize(gNumberOfScintillators);
  fPreStepy.resize(gNumberOfScintillators);
  fPreStepz.resize(gNumberOfScintillators);
  fPreStept.resize(gNumberOfScintillators);
  
  fPostStepx.resize(gNumberOfScintillators);
  fPostStepy.resize(gNumberOfScintillators);
  fPostStepz.resize(gNumberOfScintillators);
  fPostStept.resize(gNumberOfScintillators);

  fParticlePx.resize(gNumberOfScintillators);
  fParticlePy.resize(gNumberOfScintillators);
  fParticlePz.resize(gNumberOfScintillators);
  fParticleTrackID.resize(gNumberOfScintillators);
  fParticleParentID.resize(gNumberOfScintillators);
  fParticleCharge.resize(gNumberOfScintillators);
  fParticleMass.resize(gNumberOfScintillators);
  fParticlePDGID.resize(gNumberOfScintillators);

  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5EmCalorimeterSD::~B5EmCalorimeterSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5EmCalorimeterSD::Initialize(G4HCofThisEvent* hce)
{
  fHitsCollection = new B5EmCalorimeterHitsCollection(fNameSD,collectionName[0]);
  if (fHCID<0) {
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
  }
  hce->AddHitsCollection(fHCID,fHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool B5EmCalorimeterSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  auto edep = step->GetTotalEnergyDeposit();

  if (edep != 0){
    auto prepoint = step -> GetPreStepPoint();
    auto postpoint = step -> GetPostStepPoint();
    G4int sid = prepoint -> GetTouchable() -> GetReplicaNumber();
    
    G4double prex = (prepoint -> GetPosition()).x();
    G4double prey = (prepoint -> GetPosition()).y();
    G4double prez = (prepoint -> GetPosition()).z();
    G4double pret = prepoint -> GetGlobalTime();
    
    G4double postx = (postpoint -> GetPosition()).x();
    G4double posty = (postpoint -> GetPosition()).y();
    G4double postz = (postpoint -> GetPosition()).z();
    G4double postt = postpoint -> GetGlobalTime();
    
    G4double x = (prex + postx)/2.0;
    G4double y = (prey + posty)/2.0;
    G4double z = (prez + postz)/2.0;
    G4double t = (pret + postt)/2.0;
    
    fEdep[sid] += edep;
    fEweightedx[sid] += x * edep;
    fEweightedy[sid] += y * edep;
    fEweightedz[sid] += z * edep;
    fEweightedt[sid] += t * edep; 
    
    if(gSaveStepLevel == true){
      fStepEdep[sid].push_back(edep);
    
      fPreStepx[sid].push_back(prex);
      fPreStepy[sid].push_back(prey);
      fPreStepz[sid].push_back(prez);
      fPreStept[sid].push_back(pret);

      fPostStepx[sid].push_back(postx);
      fPostStepy[sid].push_back(posty);
      fPostStepz[sid].push_back(postz);
      fPostStept[sid].push_back(postt);

      // particle info
      G4Track *tr = step -> GetTrack();
      const G4ParticleDefinition *pdef = tr -> GetParticleDefinition();
      G4ThreeVector pp = tr -> GetMomentum();
    
      G4double ppx = pp.x();
      G4double ppy = pp.y();
      G4double ppz = pp.z();
      G4int trackid = tr -> GetTrackID();
      G4int parentid = tr -> GetParentID();
      G4double pcharge = pdef -> GetPDGCharge();
      G4double pmass = pdef -> GetPDGMass();
      G4int pid = pdef -> GetPDGEncoding();
    
      fParticlePx[sid].push_back(ppx);
      fParticlePy[sid].push_back(ppy);
      fParticlePz[sid].push_back(ppz);
      fParticleTrackID[sid].push_back(trackid);
      fParticleParentID[sid].push_back(parentid);
      fParticleCharge[sid].push_back(pcharge);
      fParticleMass[sid].push_back(pmass);
      fParticlePDGID[sid].push_back(pid);
    }
  }
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5EmCalorimeterSD::EndOfEvent(G4HCofThisEvent* hce){
  for (int i = 0; i < gNumberOfScintillators; i++){
    if (fEdep[i] == 0) continue;
    
    fEweightedx[i]/=fEdep[i];
    fEweightedy[i]/=fEdep[i];
    fEweightedz[i]/=fEdep[i];
    fEweightedt[i]/=fEdep[i];

    fHitsCollection->insert(new B5EmCalorimeterHit(fHCID));
    G4int CurrentHitID = fHitsCollection->GetSize()-1;
    auto hit = (B5EmCalorimeterHit*) ((hce -> GetHC(fHCID)) -> GetHit(CurrentHitID));
    
    hit -> SetXYZTE(fEweightedx[i], fEweightedy[i], fEweightedz[i], fEweightedt[i], fEdep[i]);
    hit -> SetLayerID(fLayerId);
    hit -> SetSegmentID(i);
    
    if(gSaveStepLevel == true){
      hit -> SetPreStepPos(fPreStepx[i],fPreStepy[i],fPreStepz[i],fPreStept[i]);
      hit -> SetPostStepPos(fPostStepx[i],fPostStepy[i],fPostStepz[i],fPostStept[i]);
      hit -> SetStepEdep(fStepEdep[i]);
      hit -> SetParticleTrackInfo(fParticlePx[i],fParticlePy[i],fParticlePz[i],fParticleTrackID[i],fParticleParentID[i],
				  fParticleCharge[i],fParticleMass[i],fParticlePDGID[i]);

      fStepEdep[i].clear();
      
      fPreStepx[i].clear();
      fPreStepy[i].clear();
      fPreStepz[i].clear();
      fPreStept[i].clear();
      
      fPostStepx[i].clear();
      fPostStepy[i].clear();
      fPostStepz[i].clear();
      fPostStept[i].clear();
      
      fParticlePx[i].clear();
      fParticlePy[i].clear();
      fParticlePz[i].clear();
      fParticleTrackID[i].clear();
      fParticleParentID[i].clear();
      fParticleCharge[i].clear();
      fParticleMass[i].clear();
      fParticlePDGID[i].clear();
    }
    fEdep[i] = 0.0; fEweightedx[i] = 0.0; fEweightedy[i] = 0.0; fEweightedz[i] = 0.0; fEweightedt[i] = 0.0;
  }
}
