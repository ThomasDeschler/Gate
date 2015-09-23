/*----------------------
   Copyright (C): OpenGATE Collaboration

This software is distributed under the terms
of the GNU Lesser General  Public Licence (LGPL)
See GATE/LICENSE.txt for further details
----------------------*/


/*!
  \class  GateDoseActor
  \author thibault.frisson@creatis.insa-lyon.fr
          laurent.guigues@creatis.insa-lyon.fr
	  david.sarrut@creatis.insa-lyon.fr

	  DoseToWater option added by Loïc Grevillot
  \date	March 2011
 */

  // WORK IN PROGRESS (september 2015) : Dose calculation in inhomogeneous volume added by Thomas DESCHLER (thomas.deschler@iphc.cnrs.fr)

#ifndef GATEDOSEACTOR_HH
#define GATEDOSEACTOR_HH

#include <G4NistManager.hh>

#include "GateVImageActor.hh"
#include "GateActorManager.hh"
#include "G4UnitsTable.hh"
#include "GateDoseActorMessenger.hh"
#include "GateImageWithStatistic.hh"


class G4EmCalculator;

class GateDoseActor : public GateVImageActor
{
 public:

  //-----------------------------------------------------------------------------
  // Actor name
  virtual ~GateDoseActor();

  FCT_FOR_AUTO_CREATOR_ACTOR(GateDoseActor)

  //-----------------------------------------------------------------------------
  // Constructs the sensor
  virtual void Construct();

  void EnableEdepImage(bool b) { mIsEdepImageEnabled = b; }
  void EnableEdepSquaredImage(bool b) { mIsEdepSquaredImageEnabled = b; }
  void EnableEdepUncertaintyImage(bool b) { mIsEdepUncertaintyImageEnabled = b; }
  void EnableDoseImage(bool b) { mIsDoseImageEnabled = b; }
  void EnableDoseSquaredImage(bool b) { mIsDoseSquaredImageEnabled = b; }
  void EnableDoseUncertaintyImage(bool b) { mIsDoseUncertaintyImageEnabled = b; }
  void EnableDoseToWaterImage(bool b) { mIsDoseToWaterImageEnabled = b; }
  void EnableDoseToWaterSquaredImage(bool b) { mIsDoseToWaterSquaredImageEnabled = b; }
  void EnableDoseToWaterUncertaintyImage(bool b) { mIsDoseToWaterUncertaintyImageEnabled = b; }
  void EnableNumberOfHitsImage(bool b) { mIsNumberOfHitsImageEnabled = b; }
  void EnableDoseNormalisationToMax(bool b);
  void EnableDoseNormalisationToIntegral(bool b);
  void EnableDoseToWaterNormalisation(bool b) { mIsDoseToWaterNormalisationEnabled = b; mDoseToWaterImage.SetScaleFactor(1.0); }

  virtual void BeginOfRunAction(const G4Run*r);
  virtual void BeginOfEventAction(const G4Event * event);

  virtual void UserSteppingActionInVoxel(const int index, const G4Step* step);
  virtual void UserPreTrackActionInVoxel(const int /*index*/, const G4Track* track);
  virtual void UserPostTrackActionInVoxel(const int /*index*/, const G4Track* /*t*/) {}

  //  Saves the data collected to the file
  virtual void SaveData();
  virtual void ResetData();

  // Scorer related
  virtual void Initialize(G4HCofThisEvent*){}
  virtual void EndOfEvent(G4HCofThisEvent*){}

  //Added for Voxel algorithm
  bool HasDaughter(const G4LogicalVolume* doseactorLV);
  std::pair<double,G4VSolid*> VolumeIteration(const G4VPhysicalVolume* motherPV,int generation,G4RotationMatrix motherRotation,G4ThreeVector motherTranslation);
  double VoxelIteration(G4VPhysicalVolume* motherPV,const int Generation,G4RotationMatrix MotherRotation,G4ThreeVector MotherTranslation,const int index);

protected:
  GateDoseActor(G4String name, G4int depth=0);
  GateDoseActorMessenger * pMessenger;

  int mCurrentEvent;
  StepHitType mUserStepHitType;

  bool mIsLastHitEventImageEnabled;
  bool mIsEdepImageEnabled;
  bool mIsEdepSquaredImageEnabled;
  bool mIsEdepUncertaintyImageEnabled;
  bool mIsDoseImageEnabled;
  bool mIsDoseSquaredImageEnabled;
  bool mIsDoseUncertaintyImageEnabled;
  bool mIsDoseToWaterImageEnabled;
  bool mIsDoseToWaterSquaredImageEnabled;
  bool mIsDoseToWaterUncertaintyImageEnabled;
  bool mIsNumberOfHitsImageEnabled;
  bool mIsDoseNormalisationEnabled;
  bool mIsDoseToWaterNormalisationEnabled;

  GateImageWithStatistic mEdepImage;
  GateImageWithStatistic mDoseImage;
  GateImageWithStatistic mDoseToWaterImage;
  GateImageInt mNumberOfHitsImage;
  GateImageInt mLastHitEventImage;

  G4String mEdepFilename;
  G4String mDoseFilename;
  G4String mDoseToWaterFilename;
  G4String mNbOfHitsFilename;

  G4EmCalculator * emcalc;

  //Added for Voxel algorithm
  std::vector<double> doselMass;
  std::vector<long int> voxelMatrix;
  std::vector<G4ThreeVector> voxelAbsoluteTranslation;
  std::vector<double> voxelCubicVolume;
  std::vector<double> voxelMass;
  double VoxelDensity;
  double VoxelVolume;
  double daughterVolume;
  bool   FirstTime;
  bool   matrixFirstTime;
  bool   mHasDaughter;
  G4VPhysicalVolume* DAPV;
  G4LogicalVolume* DALV;
  G4String space;

};

MAKE_AUTO_CREATOR_ACTOR(DoseActor,GateDoseActor)

#endif /* end #define GATESIMULATIONSTATISTICACTOR_HH */
