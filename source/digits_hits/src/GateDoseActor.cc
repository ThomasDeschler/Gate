/*----------------------
  Copyright (C): OpenGATE Collaboration

  This software is distributed under the terms
  of the GNU Lesser General  Public Licence (LGPL)
  See GATE/LICENSE.txt for further details
  ----------------------*/

  // WORK IN PROGRESS (september 2015) : Dose calculation in inhomogeneous volume added by Thomas DESCHLER (thomas.deschler@iphc.cnrs.fr)

/*
  \brief Class GateDoseActor :
  \brief
*/

// gate
#include "GateDoseActor.hh"
#include "GateMiscFunctions.hh"

// g4
#include <G4EmCalculator.hh>
#include <G4VoxelLimits.hh>
#include <G4NistManager.hh>
#include <G4PhysicalConstants.hh>

// Added for algorithm
#include <G4IntersectionSolid.hh>
#include <G4UnionSolid.hh>
#include <G4SubtractionSolid.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4VSolid.hh>
#include <G4Transform3D.hh>

//-----------------------------------------------------------------------------
GateDoseActor::GateDoseActor(G4String name, G4int depth):
  GateVImageActor(name,depth) {
  GateDebugMessageInc("Actor",4,"GateDoseActor() -- begin\n");

  mCurrentEvent=-1;
  mIsEdepImageEnabled = false;
  mIsLastHitEventImageEnabled = false;
  mIsEdepSquaredImageEnabled = false;
  mIsEdepUncertaintyImageEnabled = false;
  mIsDoseImageEnabled = true;
  mIsDoseSquaredImageEnabled = false;
  mIsDoseUncertaintyImageEnabled = false;
  mIsDoseToWaterImageEnabled = false;
  mIsDoseToWaterSquaredImageEnabled = false;
  mIsDoseToWaterUncertaintyImageEnabled = false;
  mIsNumberOfHitsImageEnabled = false;
  mIsDoseNormalisationEnabled = false;
  mIsDoseToWaterNormalisationEnabled = false;

  pMessenger = new GateDoseActorMessenger(this);
  GateDebugMessageDec("Actor",4,"GateDoseActor() -- end\n");
  emcalc = new G4EmCalculator;

  //Added for algorithm
  FirstTime=true;
  space="";
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
/// Destructor
GateDoseActor::~GateDoseActor()  {
  delete pMessenger;
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void GateDoseActor::EnableDoseNormalisationToMax(bool b) {
  mIsDoseNormalisationEnabled = b;
  mDoseImage.SetNormalizeToMax(b);
  mDoseImage.SetScaleFactor(1.0);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void GateDoseActor::EnableDoseNormalisationToIntegral(bool b) {
  mIsDoseNormalisationEnabled = b;
  mDoseImage.SetNormalizeToIntegral(b);
  mDoseImage.SetScaleFactor(1.0);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
/// Construct
void GateDoseActor::Construct() {
  GateDebugMessageInc("Actor", 4, "GateDoseActor -- Construct - begin\n");
  GateVImageActor::Construct();

  // Find G4_WATER. This it needed here because we will used this
  // material for dedx computation for DoseToWater.
  G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");

  // Record the stepHitType
  mUserStepHitType = mStepHitType;

  // Enable callbacks
  EnableBeginOfRunAction(true);
  EnableBeginOfEventAction(true);
  EnablePreUserTrackingAction(true);
  EnableUserSteppingAction(true);

  // Check if at least one image is enabled
  if (!mIsEdepImageEnabled &&
      !mIsDoseImageEnabled &&
      !mIsDoseToWaterImageEnabled &&
      !mIsNumberOfHitsImageEnabled)  {
    GateError("The DoseActor " << GetObjectName()
              << " does not have any image enabled ...\n Please select at least one ('enableEdep true' for example)");
  }

  // Output Filename
  mEdepFilename = G4String(removeExtension(mSaveFilename))+"-Edep."+G4String(getExtension(mSaveFilename));
  mDoseFilename = G4String(removeExtension(mSaveFilename))+"-Dose."+G4String(getExtension(mSaveFilename));
  mDoseToWaterFilename = G4String(removeExtension(mSaveFilename))+"-DoseToWater."+G4String(getExtension(mSaveFilename));
  mNbOfHitsFilename = G4String(removeExtension(mSaveFilename))+"-NbOfHits."+G4String(getExtension(mSaveFilename));

  // Set origin, transform, flag
  SetOriginTransformAndFlagToImage(mEdepImage);
  SetOriginTransformAndFlagToImage(mDoseImage);
  SetOriginTransformAndFlagToImage(mNumberOfHitsImage);
  SetOriginTransformAndFlagToImage(mLastHitEventImage);
  SetOriginTransformAndFlagToImage(mDoseToWaterImage);

  // Resize and allocate images
  if (mIsEdepSquaredImageEnabled || mIsEdepUncertaintyImageEnabled ||
      mIsDoseSquaredImageEnabled || mIsDoseUncertaintyImageEnabled ||
      mIsDoseToWaterSquaredImageEnabled || mIsDoseToWaterUncertaintyImageEnabled) {
    mLastHitEventImage.SetResolutionAndHalfSize(mResolution, mHalfSize, mPosition);
    mLastHitEventImage.Allocate();
    mIsLastHitEventImageEnabled = true;
  }
  if (mIsEdepImageEnabled) {
    //  mEdepImage.SetLastHitEventImage(&mLastHitEventImage);
    mEdepImage.EnableSquaredImage(mIsEdepSquaredImageEnabled);
    mEdepImage.EnableUncertaintyImage(mIsEdepUncertaintyImageEnabled);
    // Force the computation of squared image if uncertainty is enabled
    if (mIsEdepUncertaintyImageEnabled) mEdepImage.EnableSquaredImage(true);
    mEdepImage.SetResolutionAndHalfSize(mResolution, mHalfSize, mPosition);
    mEdepImage.Allocate();
    mEdepImage.SetFilename(mEdepFilename);
  }
  if (mIsDoseImageEnabled) {
    // mDoseImage.SetLastHitEventImage(&mLastHitEventImage);
    mDoseImage.EnableSquaredImage(mIsDoseSquaredImageEnabled);
    mDoseImage.EnableUncertaintyImage(mIsDoseUncertaintyImageEnabled);
    mDoseImage.SetResolutionAndHalfSize(mResolution, mHalfSize, mPosition);
    // Force the computation of squared image if uncertainty is enabled
    if (mIsDoseUncertaintyImageEnabled) mDoseImage.EnableSquaredImage(true);

    // DD(mDoseImage.GetVoxelVolume());
    //mDoseImage.SetScaleFactor(1e12/mDoseImage.GetVoxelVolume());
    mDoseImage.Allocate();
    mDoseImage.SetFilename(mDoseFilename);
  }
  if (mIsDoseToWaterImageEnabled) {
    mDoseToWaterImage.EnableSquaredImage(mIsDoseToWaterSquaredImageEnabled);
    mDoseToWaterImage.EnableUncertaintyImage(mIsDoseToWaterUncertaintyImageEnabled);
    // Force the computation of squared image if uncertainty is enabled
    if (mIsDoseToWaterUncertaintyImageEnabled) mDoseToWaterImage.EnableSquaredImage(true);
    mDoseToWaterImage.SetResolutionAndHalfSize(mResolution, mHalfSize, mPosition);
    mDoseToWaterImage.Allocate();
    mDoseToWaterImage.SetFilename(mDoseToWaterFilename);
  }
  if (mIsNumberOfHitsImageEnabled) {
    mNumberOfHitsImage.SetResolutionAndHalfSize(mResolution, mHalfSize, mPosition);
    mNumberOfHitsImage.Allocate();
  }

  // Print information
  GateMessage("Actor", 1,
              "\tDose DoseActor    = '" << GetObjectName() << "'\n" <<
              "\tDose image        = " << mIsDoseImageEnabled << Gateendl <<
              "\tDose squared      = " << mIsDoseSquaredImageEnabled << Gateendl <<
              "\tDose uncertainty  = " << mIsDoseUncertaintyImageEnabled << Gateendl <<
              "\tDose to water image        = " << mIsDoseToWaterImageEnabled << Gateendl <<
              "\tDose to water squared      = " << mIsDoseToWaterSquaredImageEnabled << Gateendl <<
              "\tDose to wateruncertainty  = " << mIsDoseToWaterUncertaintyImageEnabled << Gateendl <<
              "\tEdep image        = " << mIsEdepImageEnabled << Gateendl <<
              "\tEdep squared      = " << mIsEdepSquaredImageEnabled << Gateendl <<
              "\tEdep uncertainty  = " << mIsEdepUncertaintyImageEnabled << Gateendl <<
              "\tNumber of hit     = " << mIsNumberOfHitsImageEnabled << Gateendl <<
              "\t     (last hit)   = " << mIsLastHitEventImageEnabled << Gateendl <<
              "\tedepFilename      = " << mEdepFilename << Gateendl <<
              "\tdoseFilename      = " << mDoseFilename << Gateendl <<
              "\tNb Hits filename  = " << mNbOfHitsFilename << Gateendl);

  ResetData();
  GateMessageDec("Actor", 4, "GateDoseActor -- Construct - end\n");
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
/// Save data
void GateDoseActor::SaveData() {
  GateVActor::SaveData(); // (not needed because done into GateImageWithStatistic)

  if (mIsEdepImageEnabled) mEdepImage.SaveData(mCurrentEvent+1);
  if (mIsDoseImageEnabled) {
    if (mIsDoseNormalisationEnabled)
      mDoseImage.SaveData(mCurrentEvent+1, true);
    else
      mDoseImage.SaveData(mCurrentEvent+1, false);
  }

  if (mIsDoseToWaterImageEnabled) {
    if (mIsDoseToWaterNormalisationEnabled)
      mDoseToWaterImage.SaveData(mCurrentEvent+1, true);
    else
      mDoseToWaterImage.SaveData(mCurrentEvent+1, false);
  }

  if (mIsLastHitEventImageEnabled) {
    mLastHitEventImage.Fill(-1); // reset
  }

  if (mIsNumberOfHitsImageEnabled) {
    mNumberOfHitsImage.Write(mNbOfHitsFilename);
  }
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void GateDoseActor::ResetData() {
  if (mIsLastHitEventImageEnabled) mLastHitEventImage.Fill(-1);
  if (mIsEdepImageEnabled) mEdepImage.Reset();
  if (mIsDoseImageEnabled) mDoseImage.Reset();
  if (mIsDoseToWaterImageEnabled) mDoseToWaterImage.Reset();
  if (mIsNumberOfHitsImageEnabled) mNumberOfHitsImage.Fill(0);
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void GateDoseActor::BeginOfRunAction(const G4Run * r) {
  GateVActor::BeginOfRunAction(r);
  GateDebugMessage("Actor", 3, "GateDoseActor -- Begin of Run\n");
  // ResetData(); // Do no reset here !! (when multiple run);
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Callback at each event
void GateDoseActor::BeginOfEventAction(const G4Event * e) {
  GateVActor::BeginOfEventAction(e);
  mCurrentEvent++;
  GateDebugMessage("Actor", 3, "GateDoseActor -- Begin of Event: "<<mCurrentEvent << Gateendl);
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void GateDoseActor::UserPreTrackActionInVoxel(const int /*index*/, const G4Track* track)
{
  if(track->GetDefinition()->GetParticleName() == "gamma") { mStepHitType = PostStepHitType; }
  else { mStepHitType = mUserStepHitType; }
}
//-----------------------------------------------------------------------------

std::pair<double,G4VSolid*> GateDoseActor::VolumeIteration(G4LogicalVolume* MotherLV,int Generation,G4RotationMatrix MotherRotation,G4ThreeVector MotherTranslation)//,int Generation)
{
  double MotherMass(0.);
  double MotherDensity(MotherLV->GetMaterial()->GetDensity());
  double MotherProgenyMass(0.);
  double DaughtersMass(0.);
  double DaughtersCubicVolume(0.);
  double MotherDaughtersMass(0.);
  G4VSolid* MotherSV(MotherLV->GetSolid());
  G4VSolid* MotherDaughterSV(MotherLV->GetSolid());

  if(Generation>0)
    MotherSV=new G4IntersectionSolid(MotherSV->GetName()+"∩"+DALV->GetSolid()->GetName(),
                                 DALV->GetSolid(),
                                 MotherSV,
                                 &MotherRotation,
                                 MotherTranslation);
  MotherDaughterSV=MotherSV;

  space+="    ";

  G4cout<<space<<"* Mother : "<<MotherLV->GetName()<<" (generation n°"<<Generation<<") :"<< G4endl;
  G4cout<<space<<"    Volume             : "<<G4BestUnit(MotherLV->GetSolid()->GetCubicVolume(),"Volume")<<G4endl;
  G4cout<<space<<"    Volume overlap DA  : "<<G4BestUnit(MotherSV->GetCubicVolume(),"Volume")<<G4endl;
  G4cout<<space<<"    Density            : "<<G4BestUnit(MotherDensity,"Volumic Mass")<<G4endl;
  G4cout<<space<<"    Has "<<MotherLV->GetNoDaughters()<<" daughter(s)"<<G4endl;

  for(int i=0;i<MotherLV->GetNoDaughters();i++)
  {
    G4VPhysicalVolume*  DaughterPV(MotherLV->GetDaughter(i));
    G4LogicalVolume*    DaughterLV(DaughterPV->GetLogicalVolume());
    G4VSolid*           DaughterSV(DaughterLV->GetSolid());
    G4VSolid*           DaughterUnionSV(DaughterLV->GetSolid());

    G4cout<<space<<"   * Daughter n°"<<i<<" : "<<DaughterLV->GetName()<<" :"<<G4endl;
    G4cout<<space<<"       Has "<<DaughterLV->GetNoDaughters()<<" daughter(s)"<<G4endl;

    G4RotationMatrix DaughterRotation(DaughterPV->GetObjectRotationValue());
    G4ThreeVector    DaughterTranslation(DaughterPV->GetObjectTranslation());

    if(Generation>0)
    {
      DaughterTranslation=MotherTranslation+DaughterPV->GetObjectTranslation().transform(MotherRotation);
      DaughterRotation=DaughterRotation.transform(MotherRotation);
    }


    // Overlap Daughter-DoseActor
    DaughterSV=new G4IntersectionSolid(DaughterSV->GetName()+"∩"+DALV->GetSolid()->GetName(),
                                       DALV->GetSolid(),
                                       DaughterSV, 
                                       &DaughterRotation, // Global rotation
                                       DaughterTranslation); //Global translation
    double DaughterProgenyMass(0.);

    double DaughterCubicVolume(DaughterSV->GetCubicVolume());
    double DaughterDensity(DaughterLV->GetMaterial()->GetDensity());

    if(DaughterLV->GetNoDaughters()>0)
    {
      std::pair<double,G4VSolid*> DaughterIteration(VolumeIteration(DaughterLV,Generation+1,DaughterRotation,DaughterTranslation));
      DaughterProgenyMass=DaughterIteration.first;
      DaughterUnionSV=DaughterIteration.second;
    }
    else
    {
      DaughterUnionSV=DaughterSV;
      DaughterProgenyMass=DaughterCubicVolume*DaughterDensity;
    }

    //Union de la mère et de l'union de sa lignée de filles.
    MotherDaughterSV=new G4UnionSolid(MotherSV->GetName()+"∪"+DaughterUnionSV->GetName(),
                                      MotherSV, // Already overlapped with DA volume
                                      DaughterUnionSV, // From DaughterIteration
                                      DaughterPV->GetObjectRotation(), // Local rotation
                                      DaughterPV->GetObjectTranslation()); // Local translation

    // Substraction Mother-Daughter
    //G4cout<<space<<"  MotherSV ("<<MotherSV->GetName()<<") volume      : "<<G4BestUnit(MotherSV->GetCubicVolume(),"Volume")<<G4endl;
    //G4cout<<space<<"  DaughterSV ("<<DaughterSV->GetName()<<") volume      : "<<G4BestUnit(DaughterSV->GetCubicVolume(),"Volume")<<G4endl;
    MotherSV=new G4SubtractionSolid(MotherSV->GetName()+"-"+DaughterSV->GetName(),
                                     MotherSV, // Already overlapped with DA volume
                                     DaughterSV);/*, // Already overlapped with DA volume
                                     DaughterPV->GetFrameRotation(), // Local rotation
                                     DaughterPV->GetFrameTranslation()); // Local translation*/
    //G4cout<<space<<"  "<<MotherSV->GetName()<<" volume      : "<<G4BestUnit(MotherSV->GetCubicVolume(),"Volume")<<G4endl;

 

    G4cout<<space<<"       Original volume   : "<<G4BestUnit(DaughterLV->GetSolid()->GetCubicVolume(),"Volume")<<G4endl
          <<space<<"       DA overlap volume : "<<G4BestUnit(DaughterCubicVolume,"Volume")<<G4endl
          <<space<<"       Density           : "<<G4BestUnit(DaughterDensity,"Volumic Mass")<<G4endl
          <<space<<"       Mass              : "<<G4BestUnit(DaughterCubicVolume*DaughterDensity,"Mass")<<G4endl;
    G4cout<<space<<"       Progeny Mass      : "<<G4BestUnit(DaughterProgenyMass,"Mass")<<G4endl;
    G4cout<<space<<"       Relative translation: X="<<G4BestUnit(DaughterPV->GetTranslation().getX(),"Length")<<",Y="<<G4BestUnit(DaughterPV->GetTranslation().getY(),"Length")<<",Z="<<G4BestUnit(DaughterPV->GetTranslation().getZ(),"Length")<<G4endl;
    G4cout<<space<<"       Absolute translation: X="<<G4BestUnit(DaughterTranslation.getX(),"Length")<<",Y="<<G4BestUnit(DaughterTranslation.getY(),"Length")<<",Z="<<G4BestUnit(DaughterTranslation.getZ(),"Length")<<G4endl;
    G4cout<<space<<"       Relative rotation   : Phi= "<<DaughterPV->GetObjectRotationValue().getPhi()<<", Theta= "<<DaughterPV->GetObjectRotationValue().getTheta()<<", Psi= "<<DaughterPV->GetObjectRotationValue().getPsi()<<G4endl;
    G4cout<<space<<"       Absolute rotation   : Phi= "<<DaughterRotation.getPhi()<<", Theta= "<<DaughterRotation.getTheta()<<", Psi= "<<DaughterRotation.getPsi()<<G4endl;

    MotherProgenyMass+=DaughterProgenyMass;
    DaughtersCubicVolume=DaughtersCubicVolume+DaughterCubicVolume;
  }

  G4cout<<space<<"  * Total daughters :"<<G4endl;
  G4cout<<space<<"    Daughters Volume : "<<G4BestUnit(DaughtersCubicVolume,"Volume")<<G4endl;

  //FIXME Faire que le volume mère soit soustrait du volume en overlap avec sa fille.
  double MotherCubicVolume(MotherSV->GetCubicVolume());
  MotherMass=MotherCubicVolume*MotherDensity;

  MotherProgenyMass+=MotherMass;

  G4cout<<space<<"  Mother original volume     : "<<G4BestUnit(MotherLV->GetSolid()->GetCubicVolume(),"Volume")<<G4endl;
  G4cout<<space<<"  Mother w/o daughters volume: "<<G4BestUnit(MotherCubicVolume,"Volume")<<G4endl;
  G4cout<<space<<"  Mother  +  daughters volume: "<<G4BestUnit(MotherCubicVolume+DaughtersCubicVolume,"Volume")<<G4endl;
  G4cout<<space<<"  Mother Mass                : "<<G4BestUnit(MotherMass,"Mass")<<G4endl;
  G4cout<<space<<"  Mother  +  daughters mass  : "<<G4BestUnit(MotherProgenyMass,"Mass")<<" (MotherProgenyMass)"<<G4endl;
  //FIXME :
  G4cout<<space<<"  MotherDaughtersUnion volume: "<<G4BestUnit(MotherDaughterSV->GetCubicVolume(),"Volume")<<G4endl;

  space.resize(space.size()-3);

  return std::make_pair(MotherProgenyMass,MotherDaughterSV);
}

//-----------------------------------------------------------------------------
void GateDoseActor::UserSteppingActionInVoxel(const int index, const G4Step* step) {
  GateDebugMessageInc("Actor", 4, "GateDoseActor -- UserSteppingActionInVoxel - begin\n");
  GateDebugMessageInc("Actor", 4, "enedepo = " << step->GetTotalEnergyDeposit() << Gateendl);
  GateDebugMessageInc("Actor", 4, "weight = " <<  step->GetTrack()->GetWeight() << Gateendl);
  const double weight = step->GetTrack()->GetWeight();
  const double edep = step->GetTotalEnergyDeposit()*weight;//*step->GetTrack()->GetWeight();

  // if no energy is deposited or energy is deposited outside image => do nothing
  if (edep == 0) {
    GateDebugMessage("Actor", 5, "edep == 0 : do nothing\n");
    GateDebugMessageDec("Actor", 4, "GateDoseActor -- UserSteppingActionInVoxel -- end\n");
    return;
  }
  if (index < 0) {
    GateDebugMessage("Actor", 5, "index < 0 : do nothing\n");
    GateDebugMessageDec("Actor", 4, "GateDoseActor -- UserSteppingActionInVoxel -- end\n");
    return;
  }

  // compute sameEvent
  // sameEvent is false the first time some energy is deposited for each primary particle
  bool sameEvent=true;
  if (mIsLastHitEventImageEnabled) {
    GateDebugMessage("Actor", 2,  "GateDoseActor -- UserSteppingActionInVoxel: Last event in index = " << mLastHitEventImage.GetValue(index) << Gateendl);
    if (mCurrentEvent != mLastHitEventImage.GetValue(index)) {
      sameEvent = false;
      mLastHitEventImage.SetValue(index, mCurrentEvent);
    }
  }

  double dose=0.;
  double density = step->GetPreStepPoint()->GetMaterial()->GetDensity();
  //double volume = step->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetSolid()->GetCubicVolume();

  if(FirstTime)
  {
    FirstTime=false;
    VoxelMass=0.;
    VoxelDensity=0.;
    VoxelVolume=0.;

    DAPV=G4PhysicalVolumeStore::GetInstance()->GetVolume(mVolumeName+"_phys");
    DALV=G4LogicalVolumeStore::GetInstance()->GetVolume(mVolumeName+"_log");

    G4cout<<G4endl<<" DEBUG of DoseActor "<< GetObjectName()<<" :"<<G4endl
                  <<" DA volume name : "<<mVolumeName<<G4endl;
    //G4cout        <<"   Number of voxels  : "<<DALV->GetVoxelHeader()->GetNoSlices()<<G4endl;
    //G4cout        <<"   Physical volume : "<<step->GetPreStepPoint()->GetPhysicalVolume()->GetName()<<G4endl;


    //if(DALV->GetNoDaughters()>0)
      VoxelMass=VolumeIteration(DALV,0,DAPV->GetObjectRotationValue(),DAPV->GetObjectTranslation()).first;
    //else VoxelMass=DALV->GetMass();
      //VolumeIteration(DALV,DAPV->GetRotation(),DAPV->GetTranslation());

    double DALVCubicVolume(DALV->GetSolid()->GetCubicVolume());
    VoxelDensity=VoxelMass/DALVCubicVolume;

    G4cout<<" DA volume           : "<<G4BestUnit(DALVCubicVolume,"Volume")<<G4endl;
    G4cout<<" DA original mass    : "<<G4BestUnit(DALV->GetMass(true,true,0),"Mass")<<" (propagated)"<<G4endl;
    G4cout<<" DA original mass    : "<<G4BestUnit(DALV->GetMass(false,false,0),"Mass")<<" (not propagated)"<<G4endl;
    G4cout<<" DA original density : "<<G4BestUnit(DALV->GetMass()/DALVCubicVolume,"Volumic Mass")<<G4endl;
    G4cout<<" DA corrected mass   : "<<G4BestUnit(VoxelMass,"Mass")<<G4endl;
    G4cout<<" DA corrected density: "<<G4BestUnit(VoxelDensity,"Volumic Mass")<<G4endl;
    G4cout<<G4endl;
  }


  if (mIsDoseImageEnabled) {

    // ------------------------------------
    // Convert deposited energy into Gray

    // OLD version (correct but not clear)
    // dose = edep/density*1e12/mDoseImage.GetVoxelVolume();

    // NEW version (same results but more clear)
    //dose = edep/density/mDoseImage.GetVoxelVolume()/gray;
    dose = edep/VoxelDensity/mDoseImage.GetVoxelVolume()/gray;
    // ------------------------------------

    GateDebugMessage("Actor", 2,  "GateDoseActor -- UserSteppingActionInVoxel:\tdose = "
		     << G4BestUnit(dose, "Dose")
		     << " rho = "
		     << G4BestUnit(density, "Volumic Mass")<< Gateendl );
  }

  double doseToWater = 0;
  if (mIsDoseToWaterImageEnabled) {

    // to get nuclear inelastic cross-section, see "geant4.9.4.p01/examples/extended/hadronic/Hadr00/"
    // #include "G4HadronicProcessStore.hh"
    // G4HadronicProcessStore* store = G4HadronicProcessStore::Instance();
    // store->GetInelasticCrossSectionPerAtom(particle,e,elm);

    double cut = DBL_MAX;
    cut=1;
    G4String material = step->GetPreStepPoint()->GetMaterial()->GetName();
    double Energy = step->GetPreStepPoint()->GetKineticEnergy();
    G4String PartName = step->GetTrack()->GetDefinition()->GetParticleName();
    //    const G4ParticleDefinition * PartDef = step->GetTrack()->GetParticleDefinition();
    //    G4Material  * MatDef = step->GetTrack()->GetMaterial();
    double DEDX=0, DEDX_Water=0;
    //    G4cout<<PartName<<"\t";//Gateendl;//"  "<<edep<<"  "<<NonIonizingEdep<< Gateendl;


    // Dose to water: it could be possible to make this process more
    // generic by choosing any material in place of water
    double Volume = mDoseToWaterImage.GetVoxelVolume();

    // Other particles should be taken into account (Helium etc), but bug ? FIXME
    if (PartName== "proton" || PartName== "e-" || PartName== "e+" || PartName== "deuteron"){
      //if (PartName != "O16[0.0]" && PartName != "alpha" && PartName != "Be7[0.0]" && PartName != "C12[0.0]"){

      DEDX = emcalc->ComputeTotalDEDX(Energy, PartName, material, cut);
      DEDX_Water = emcalc->ComputeTotalDEDX(Energy, PartName, "G4_WATER", cut);

      doseToWater=edep/density/Volume/gray*(DEDX_Water/1.)/(DEDX/(density*e_SI));

    }
    else {
      DEDX = emcalc->ComputeTotalDEDX(100, "proton", material, cut);
      DEDX_Water = emcalc->ComputeTotalDEDX(100, "proton", "G4_WATER", cut);
      doseToWater=edep/density/Volume/gray*(DEDX_Water/1.)/(DEDX/(density*e_SI));
    }

    GateDebugMessage("Actor", 2,  "GateDoseActor -- UserSteppingActionInVoxel:\tdose to water = "
		     << G4BestUnit(doseToWater, "Dose to water")
		     << " rho = "
		     << G4BestUnit(density, "Volumic Mass")<< Gateendl );
  }


  if (mIsEdepImageEnabled) {
    GateDebugMessage("Actor", 2, "GateDoseActor -- UserSteppingActionInVoxel:\tedep = " << G4BestUnit(edep, "Energy") << Gateendl);
  }



  if (mIsDoseImageEnabled) {

    if (mIsDoseUncertaintyImageEnabled || mIsDoseSquaredImageEnabled) {
      if (sameEvent) mDoseImage.AddTempValue(index, dose);
      else mDoseImage.AddValueAndUpdate(index, dose);
    }
    else mDoseImage.AddValue(index, dose);
  }

  if (mIsDoseToWaterImageEnabled) {

    if (mIsDoseToWaterUncertaintyImageEnabled || mIsDoseToWaterSquaredImageEnabled) {
      if (sameEvent) mDoseToWaterImage.AddTempValue(index, doseToWater);
      else mDoseToWaterImage.AddValueAndUpdate(index, doseToWater);
    }
    else mDoseToWaterImage.AddValue(index, doseToWater);
  }

  if (mIsEdepImageEnabled) {
    if (mIsEdepUncertaintyImageEnabled || mIsEdepSquaredImageEnabled) {
      if (sameEvent) mEdepImage.AddTempValue(index, edep);
      else mEdepImage.AddValueAndUpdate(index, edep);
    }
    else mEdepImage.AddValue(index, edep);
  }

  if (mIsNumberOfHitsImageEnabled) mNumberOfHitsImage.AddValue(index, weight);

  GateDebugMessageDec("Actor", 4, "GateDoseActor -- UserSteppingActionInVoxel -- end\n");
}
//-----------------------------------------------------------------------------
