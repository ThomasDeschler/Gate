/*----------------------
   Copyright (C): OpenGATE Collaboration

This software is distributed under the terms
of the GNU Lesser General  Public Licence (LGPL)
See GATE/LICENSE.txt for further details
----------------------*/


/*!
  \class  GateVoxelizedMassActor
  \author Thomas DESCHLER (thomas.deschler@iphc.cnrs.fr)
  \date	October 2015
 */

#ifndef GATEVOXELIZEDMASS_HH
#define GATEVOXELIZEDMASS_HH

#include "GateVImageActor.hh"
#include "GateActorManager.hh"
#include "G4UnitsTable.hh"
#include "GateVoxelizedMassActorMessenger.hh"
#include "GateImageWithStatistic.hh"

#include <G4Box.hh>

class GateVoxelizedMass
{
 public:

  GateVoxelizedMass();

  virtual ~GateVoxelizedMass() {}

  void Initialize(const G4String mExtVolumeName, const GateImageDouble mExtImage);
  double GetVoxelMass(const int index);
  std::vector<double> GetVoxelMassVector();
  virtual void GenerateVectors();
  virtual void GenerateVoxels();
  virtual void ParameterizedVolume(const int index);
  virtual std::pair<double,double> VoxelIteration(G4VPhysicalVolume* motherPV,const int Generation,G4RotationMatrix MotherRotation,G4ThreeVector MotherTranslation,const int index);

 protected:

  G4VPhysicalVolume* DAPV;
  G4LogicalVolume* DALV;
  G4Box* doselSV;
  G4Box* DABox;
  G4Box* voxelBox;

  double doselReconstructedTotalCubicVolume;
  double doselReconstructedTotalMass;
  std::pair<double,double> doselReconstructedData;
  std::vector<double> doselReconstructedCubicVolume;
  std::vector<double> doselReconstructedMass;

  std::vector<std::vector<std::vector<double> > > voxelCubicVolume;
  std::vector<std::vector<std::vector<double> > > voxelMass;

  GateImageDouble mImage;
  G4String mVolumeName;
  bool mIsParameterised;
  bool mIsVecGenerated;

  G4String space;
  int seconds;

};

#endif
