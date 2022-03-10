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
/// \file B4PrimaryGeneratorAction.cc
/// \brief Implementation of the B4PrimaryGeneratorAction class
#include "TF2.h"
#include "TH2.h"
#include "TMath.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

#include <iostream>
#include <iomanip>
#include <math.h>

#include "B4PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4UnitsTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleMomentum.hh"
#include "TRandom.h"



Double_t muI(Double_t *x, Double_t *par) {

  Double_t p=x[0]; // p

  Double_t ksi=p*x[1]; // p*cost

  const Double_t c1=0.00253;
  const Double_t c2=0.2455;
  const Double_t c3=1.288;
  const Double_t c4=-0.2555;
  const Double_t c5=0.0209;

  Double_t lksi=TMath::Log10(ksi);

  Double_t Iv = par[0]*c1*pow(ksi,-1.*(c2+c3*lksi+c4*lksi*lksi+c5*lksi*lksi*lksi));

  return Iv*x[1]*x[1]*x[1];


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4PrimaryGeneratorAction::B4PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(nullptr)
{
  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);

  gRandom->SetSeed(0);

const Int_t npar = 1;


f2 = new TF2("f2",muI,0,2000,0,1, npar);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4PrimaryGeneratorAction::~B4PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  
  using namespace std; 

////////////////////////////////////////////////////////////////////////////////////
////////// For angle and momentum muons ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

const Int_t npar = 1;

 Double_t f2params[npar] = {0.8915};

 Double_t cost, p;

f2->SetParameters(f2params);

 f2->SetNpx(5000);

   f2->GetRandom2(p, cost ); 
 
p=p*GeV;

G4cout<<"                                                  cos()   ===  "<<cost<<G4endl;
G4cout<<"                                                  impuls  ===  "<<p<<G4endl;


std::ofstream output_file;

double sint=1.-cost*cost;

   double phi=2.*TMath::Pi()*gRandom->Rndm();


    double xDir=sint*std::cos(phi);
    double yDir=sint*std::sin(phi);
    double zDir=cost;

    double pz = (p*zDir);
    double px = (p*xDir);
    double py = (p*yDir);


 G4double CMEn = sqrt(pow(105.6583745,2)+pow(p,2)); 

  auto particleDefinition 
    = G4ParticleTable::GetParticleTable()->FindParticle("mu-");
  fParticleGun->SetParticleDefinition(particleDefinition);
  fParticleGun->SetParticleEnergy(CMEn);
  

  fParticleGun->SetParticleMomentumDirection(G4ThreeVector( px, py, pz));

  // This function is called at the begining of event

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get world volume 
  // from G4LogicalVolumeStore
  //
  G4double worldZHalfLength = 0.;
  auto worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");

  // Check that the world volume has box shape
  G4Box* worldBox = nullptr;
  if (  worldLV ) {
    worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
  }

  if ( worldBox ) {
    worldZHalfLength = worldBox->GetZHalfLength();  
  }
  else  {
    G4ExceptionDescription msg;
    msg << "World volume of box shape not found." << G4endl;
    msg << "Perhaps you have changed geometry." << G4endl;
    msg << "The gun will be place in the center.";
    G4Exception("B4PrimaryGeneratorAction::GeneratePrimaries()",
      "MyCode0002", JustWarning, msg);
  } 


        G4double positionX = G4UniformRand();
        G4double positionY = G4UniformRand();
        G4double XX = G4UniformRand()*1.25*m;
        G4double YY = G4UniformRand()*1.4*m;
   if (positionX >= 0.5)
   {
    XX = XX*(-1);
   }
    if (positionY >= 0.5)
   {
    YY = YY*(-1);
   }
    
  // Set gun position
  fParticleGun
    ->SetParticlePosition(G4ThreeVector(XX, YY, -worldZHalfLength));

  fParticleGun->GeneratePrimaryVertex(anEvent);

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......