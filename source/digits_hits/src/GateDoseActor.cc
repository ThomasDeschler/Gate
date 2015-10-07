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
#include <G4Box.hh>
#include <G4Transform3D.hh>
#include <G4SmartVoxelHeader.hh>
#include <G4VPVParameterisation.hh>
#include <time.h>

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
  mIsMassImageEnabled = false;

  pMessenger = new GateDoseActorMessenger(this);
  GateDebugMessageDec("Actor",4,"GateDoseActor() -- end\n");
  emcalc = new G4EmCalculator;

  //Added for algorithm
  FirstTime=true;
  matrixFirstTime=true;
  mMassFirstTime=true;
  space="";
  nbofRecVoxels=0;
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
      !mIsNumberOfHitsImageEnabled &&
      !mIsMassImageEnabled)  {
    GateError("The DoseActor " << GetObjectName()
              << " does not have any image enabled ...\n Please select at least one ('enableEdep true' for example)");
  }

  // Output Filename
  mEdepFilename = G4String(removeExtension(mSaveFilename))+"-Edep."+G4String(getExtension(mSaveFilename));
  mDoseFilename = G4String(removeExtension(mSaveFilename))+"-Dose."+G4String(getExtension(mSaveFilename));
  mDoseToWaterFilename = G4String(removeExtension(mSaveFilename))+"-DoseToWater."+G4String(getExtension(mSaveFilename));
  mNbOfHitsFilename = G4String(removeExtension(mSaveFilename))+"-NbOfHits."+G4String(getExtension(mSaveFilename));
  mMassFilename = G4String(removeExtension(mSaveFilename))+"-Mass."+G4String(getExtension(mSaveFilename));

  // Set origin, transform, flag
  SetOriginTransformAndFlagToImage(mEdepImage);
  SetOriginTransformAndFlagToImage(mDoseImage);
  SetOriginTransformAndFlagToImage(mNumberOfHitsImage);
  SetOriginTransformAndFlagToImage(mLastHitEventImage);
  SetOriginTransformAndFlagToImage(mDoseToWaterImage);
  SetOriginTransformAndFlagToImage(mMassImage);

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
  if (mIsMassImageEnabled) {
    mMassImage.SetResolutionAndHalfSize(mResolution, mHalfSize, mPosition);
    mMassImage.Allocate();
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
              "\tMass image        = " << mIsMassImageEnabled << Gateendl <<
              "\tEdepFilename      = " << mEdepFilename << Gateendl <<
              "\tDoseFilename      = " << mDoseFilename << Gateendl <<
              "\tNb Hits filename  = " << mNbOfHitsFilename << Gateendl <<
              "\tMassFilename  = " << mMassFilename << Gateendl);

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

  if (mIsMassImageEnabled) {
    mMassImage.Write(mMassFilename);
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
  if (mIsMassImageEnabled) mMassImage.Fill(0);
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

double RealZero(double someValue,double anotherValue=0.,double relativeError=0.05)
{
  if(anotherValue!=0.)
  {
    if(someValue<(anotherValue*relativeError)
    && someValue>(-anotherValue*relativeError))
      someValue=0.0;
    else if(someValue-anotherValue< anotherValue*relativeError
         && someValue-anotherValue>-anotherValue*relativeError)
      someValue=anotherValue;
  }
  else if(someValue< std::numeric_limits<double>::epsilon()
        &&someValue>-std::numeric_limits<double>::epsilon())
      someValue=0.0;
  return someValue;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
std::pair<double,G4VSolid*> GateDoseActor::VolumeIteration(const G4VPhysicalVolume* motherPV,int Generation,G4RotationMatrix MotherRotation,G4ThreeVector MotherTranslation)
{
  G4LogicalVolume* MotherLV(motherPV->GetLogicalVolume());
  G4VSolid*        MotherSV(MotherLV->GetSolid());
  G4VSolid*         MotherDaughterSV(MotherLV->GetSolid());
  double MotherMass(0.);
  double MotherDensity(MotherLV->GetMaterial()->GetDensity());
  double MotherProgenyMass(0.);
  //double DaughtersMass(0.);
  double DaughtersCubicVolume(0.);
  //double MotherDaughtersMass(0.);


  if(Generation>0)
    MotherSV=new G4IntersectionSolid(MotherSV->GetName()+"∩"+DALV->GetSolid()->GetName(),
                                 DALV->GetSolid(),
                                 MotherSV,
                                 &MotherRotation,
                                 MotherTranslation);
  MotherDaughterSV=MotherSV;

  space+="   ";

  if(MotherLV->GetSolid()->GetCubicVolume()<MotherSV->GetCubicVolume())
    GateWarning(MotherLV->GetSolid()->GetName()<<" overlaps the DoseActor associated volume ("<<DALV->GetSolid()->GetName()<<") !"<<Gateendl
              <<"               The results may not be accurate.");

  G4cout<<space<<"* Mother : "<<MotherLV->GetName()<<" (generation n°"<<Generation<<") :"<<G4endl
        <<space<<"    Volume             : "<<G4BestUnit(MotherLV->GetSolid()->GetCubicVolume(),"Volume")<<G4endl
        <<space<<"    Volume overlap DA  : "<<G4BestUnit(MotherSV->GetCubicVolume(),"Volume")<<G4endl
        <<space<<"    Density            : "<<G4BestUnit(MotherDensity,"Volumic Mass")<<G4endl
        <<space<<"    Parameterised ?    : "<<motherPV->IsParameterised()<<G4endl
        <<space<<"    Has "<<MotherLV->GetNoDaughters()<<" daughter(s)"<<G4endl;

  for(int i=0;i<MotherLV->GetNoDaughters();i++)
  {
    G4VPhysicalVolume*  daughterPV(MotherLV->GetDaughter(i));
    G4LogicalVolume*    DaughterLV(daughterPV->GetLogicalVolume());
    G4VSolid*           DaughterSV(DaughterLV->GetSolid());
    G4VSolid*           DaughterUnionSV(DaughterLV->GetSolid());

    G4cout<<space<<"   * Daughter n°"<<i<<" : "<<DaughterLV->GetName()<<" :"<<G4endl
          <<space<<"       Has "<<DaughterLV->GetNoDaughters()<<" daughter(s)"<<G4endl;

    G4RotationMatrix DaughterRotation(daughterPV->GetObjectRotationValue());
    G4ThreeVector    DaughterTranslation(daughterPV->GetObjectTranslation());

    if(Generation>0)
    {
      DaughterTranslation=MotherTranslation+daughterPV->GetObjectTranslation().transform(MotherRotation);
      DaughterRotation=DaughterRotation.transform(MotherRotation);
    }


    // Overlap Daughter-DoseActor
    DaughterSV=new G4IntersectionSolid(DaughterSV->GetName()+"∩"+DALV->GetSolid()->GetName(),
                                       DALV->GetSolid(),
                                       DaughterSV, 
                                       &DaughterRotation, // Global rotation
                                       DaughterTranslation); //Global translation

    if(MotherLV->GetDaughter(i)->GetLogicalVolume()->GetSolid()->GetCubicVolume()>DaughterSV->GetCubicVolume())
      GateWarning(MotherLV->GetDaughter(i)->GetLogicalVolume()->GetSolid()->GetName()<<" overlaps the DoseActor associated volume ("<<DALV->GetSolid()->GetName()<<") !"<<Gateendl
                <<"               The results may not be accurate.");

    double DaughterProgenyMass(0.);

    double DaughterCubicVolume(DaughterSV->GetCubicVolume());
    double DaughterDensity(DaughterLV->GetMaterial()->GetDensity());

    if(DaughterLV->GetNoDaughters()>0)
    {
      std::pair<double,G4VSolid*> DaughterIteration(VolumeIteration(daughterPV,Generation+1,DaughterRotation,DaughterTranslation));
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
                                      daughterPV->GetObjectRotation(), // Local rotation
                                      daughterPV->GetObjectTranslation()); // Local translation

    // Substraction Mother-Daughter
    //G4cout<<space<<"  MotherSV ("<<MotherSV->GetName()<<") volume      : "<<G4BestUnit(MotherSV->GetCubicVolume(),"Volume")<<G4endl;
    //G4cout<<space<<"  DaughterSV ("<<DaughterSV->GetName()<<") volume      : "<<G4BestUnit(DaughterSV->GetCubicVolume(),"Volume")<<G4endl;
    MotherSV=new G4SubtractionSolid(MotherSV->GetName()+"-"+DaughterSV->GetName(),
                                     MotherSV, // Already overlapped with DA volume
                                     DaughterSV);/*, // Already overlapped with DA volume
                                     daughterPV->GetFrameRotation(), // Local rotation
                                     daughterPV->GetFrameTranslation()); // Local translation*/
    //G4cout<<space<<"  "<<MotherSV->GetName()<<" volume      : "<<G4BestUnit(MotherSV->GetCubicVolume(),"Volume")<<G4endl;

 

    G4cout<<space<<"       Original volume   : "<<G4BestUnit(DaughterLV->GetSolid()->GetCubicVolume(),"Volume")<<G4endl
          <<space<<"       DA overlap volume   : "<<G4BestUnit(DaughterCubicVolume,"Volume")<<G4endl
          <<space<<"       Density             : "<<G4BestUnit(DaughterDensity,"Volumic Mass")<<G4endl
          <<space<<"       Mass                : "<<G4BestUnit(DaughterCubicVolume*DaughterDensity,"Mass")<<G4endl
          <<space<<"       Progeny Mass        : "<<G4BestUnit(DaughterProgenyMass,"Mass")<<G4endl
          <<space<<"       Relative translation: X="<<G4BestUnit(daughterPV->GetTranslation().getX(),"Length")<<",Y="<<G4BestUnit(daughterPV->GetTranslation().getY(),"Length")<<",Z="<<G4BestUnit(daughterPV->GetTranslation().getZ(),"Length")<<G4endl
          <<space<<"       Absolute translation: X="<<G4BestUnit(DaughterTranslation.getX(),"Length")<<",Y="<<G4BestUnit(DaughterTranslation.getY(),"Length")<<",Z="<<G4BestUnit(DaughterTranslation.getZ(),"Length")<<G4endl
          <<space<<"       Relative rotation   : Phi= "<<daughterPV->GetObjectRotationValue().getPhi()<<", Theta= "<<daughterPV->GetObjectRotationValue().getTheta()<<", Psi= "<<daughterPV->GetObjectRotationValue().getPsi()<<G4endl
          <<space<<"       Absolute rotation   : Phi= "<<DaughterRotation.getPhi()<<", Theta= "<<DaughterRotation.getTheta()<<", Psi= "<<DaughterRotation.getPsi()<<G4endl;


    G4cout<<space<<"       Parameterised ?     : "<<daughterPV->IsParameterised()<<G4endl;
    G4cout<<space<<"       Replicated ?        : "<<daughterPV->IsReplicated()<<G4endl;
    if(daughterPV->IsParameterised())
    {
      G4cout<<space<<"       daughterPV->GetName()         : "<<daughterPV->GetName() <<G4endl;
      G4cout<<space<<"       daughterPV->GetMultiplicity() : "<<daughterPV->GetMultiplicity() <<G4endl;
      G4VPVParameterisation* daughterParameterisation(daughterPV->GetParameterisation());
      double daughterParameterisationMass(0.);
      for(long int i=0;i<daughterPV->GetMultiplicity();i++)
      {
        G4VSolid* daughterSV(daughterParameterisation->ComputeSolid(i,daughterPV));
        daughterParameterisation->ComputeTransformation(i,daughterPV);
        daughterParameterisationMass+=daughterParameterisation->ComputeMaterial(i,daughterPV)->GetDensity()*daughterSV->GetCubicVolume();

        /*G4cout<<space<<"       Voxel n°"<<i<<" :"<<G4endl;
        G4cout<<space<<"         CubicVolume : "<<G4BestUnit(daughterParameterisation->ComputeSolid(i,daughterPV)->GetCubicVolume(),"Volume")<<G4endl;
        G4cout<<space<<"         Material    : "<<daughterParameterisation->ComputeMaterial(i,daughterPV)->GetName()<<G4endl;
        G4cout<<space<<"         Translation : X="<<G4BestUnit(daughterPV->GetTranslation().getX(),"Length")<<",Y="<<G4BestUnit(daughterPV->GetTranslation().getY(),"Length")<<",Z="<<G4BestUnit(daughterPV->GetTranslation().getZ(),"Length")<<G4endl;*/
      }
      G4cout<<space<<"       daughterPV total mass : "<<G4BestUnit(daughterParameterisationMass,"Mass")<<G4endl;
      //G4cout<<space<<"       voxelHeader->AllSlicesEqual () :"<<voxelHeader->AllSlicesEqual () <<G4endl; // Seg fault
    }





    MotherProgenyMass+=DaughterProgenyMass;
    DaughtersCubicVolume=DaughtersCubicVolume+DaughterCubicVolume;
  }

  //FIXME :
  //G4cout<<space<<"    Daughters Volume : "<<G4BestUnit(DaughtersCubicVolume,"Volume")<<G4endl;

  double MotherCubicVolume(MotherSV->GetCubicVolume());
  MotherMass=MotherCubicVolume*MotherDensity;

  MotherProgenyMass+=MotherMass;

  GateDebugMessage("Actor", 5,space<<"  Mother original volume     : "<<G4BestUnit(MotherLV->GetSolid()->GetCubicVolume(),"Volume")<<Gateendl
        <<space<<"  Mother w/o daughters volume: "<<G4BestUnit(MotherCubicVolume,"Volume")<<Gateendl
        <<space<<"  Mother  +  daughters volume: "<<G4BestUnit(MotherCubicVolume+DaughtersCubicVolume,"Volume")<<Gateendl
        <<space<<"  Mother Mass                : "<<G4BestUnit(MotherMass,"Mass")<<Gateendl
        <<space<<"  Mother  +  daughters mass  : "<<G4BestUnit(MotherProgenyMass,"Mass")<<" (MotherProgenyMass)"<<Gateendl
  //FIXME :
        <<space<<"  MotherDaughtersUnion volume: "<<G4BestUnit(MotherDaughterSV->GetCubicVolume(),"Volume")<<Gateendl);

  space.resize(space.size()-3);

  return std::make_pair(MotherProgenyMass,MotherDaughterSV);
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
double GateDoseActor::VoxelIteration(G4VPhysicalVolume* motherPV,const int Generation,G4RotationMatrix motherRotation,G4ThreeVector motherTranslation,const int index)
{
  //FIXME : Doesn't work with daughter overlapping its mother.
  if(Generation>0)
    space+="   ";

  G4LogicalVolume* motherLV(motherPV->GetLogicalVolume());
  double motherMass(0.);
  double motherDensity(motherLV->GetMaterial()->GetDensity());
  double motherProgenyMass(0.);
  double motherProgenyCubicVolume(0.);
  G4VSolid* motherSV(motherLV->GetSolid());

  G4Box doselSV("DoselSV",//FIXME convertir i en string
                mDoseImage.GetValueImage().GetVoxelSize().getX()/2.0,
                mDoseImage.GetValueImage().GetVoxelSize().getY()/2.0,
                mDoseImage.GetValueImage().GetVoxelSize().getZ()/2.0);

  // Calculation of voxel's local rotation and translation
  G4RotationMatrix doselRotation(mDoseImage.GetValueImage().GetTransformMatrix());
  G4ThreeVector    doselTranslation(mDoseImage.GetValueImage().GetVoxelCenterFromIndex(index));

  if(Generation>0)
    doselTranslation-=motherTranslation;

  // Dosel information dump ////////////////////////////////////////////////
  /*space+="   ";
  G4cout<<space<<"* Doxel n°"<<index<<" informations :"<<G4endl
        <<space<<" --> Translation :"<<G4endl
        <<space<<"       Relative : X="<<G4BestUnit(doselTranslation.getX(),"Length")
                                <<",Y="<<G4BestUnit(doselTranslation.getY(),"Length")
                                <<",Z="<<G4BestUnit(doselTranslation.getZ(),"Length")<<G4endl
        <<space<<"       Absolute : X="<<G4BestUnit(mDoseImage.GetValueImage().GetVoxelCenterFromIndex(index).getX(),"Length")
                                <<",Y="<<G4BestUnit(mDoseImage.GetValueImage().GetVoxelCenterFromIndex(index).getY(),"Length")
                                <<",Z="<<G4BestUnit(mDoseImage.GetValueImage().GetVoxelCenterFromIndex(index).getZ(),"Length")<<G4endl
        <<space<<" --> Rotation :"<<G4endl
        <<space<<"       Relative : Phi="  <<doselRotation.getPhi()
                                <<",Theta="<<doselRotation.getTheta()
                                <<",Psi="  <<doselRotation.getPsi()<<G4endl
        <<space<<"       Absolute : Phi="<<mDoseImage.GetValueImage().GetTransformMatrix().getPhi()
                                <<",Theta="<<mDoseImage.GetValueImage().GetTransformMatrix().getTheta()
                                <<",Psi="<<mDoseImage.GetValueImage().GetTransformMatrix().getPsi()<<G4endl;
  space.resize(space.size()-3);*/
  //////////////////////////////////////////////////////////////////////////

  // Parameterisation //////////////////////////////////////////////////////
  /*G4cout<<space<<"     * "<<motherPV->GetName()<<" :"<<G4endl;
  G4cout<<space<<"       Parameterised ?     : "<<motherPV->IsParameterised()<<G4endl;
  G4cout<<space<<"       Replicated ?        : "<<motherPV->IsReplicated()<<G4endl;*/
  double doselCubicVolume(doselSV.GetCubicVolume());
  if(motherPV->IsParameterised())
  {
    //G4cout<<space<<"       motherPV->GetName()         : "<<motherPV->GetName() <<G4endl;
    //G4cout<<space<<"       motherPV->GetMultiplicity() : "<<motherPV->GetMultiplicity() <<G4endl;
    G4VPVParameterisation* motherParameterisation(motherPV->GetParameterisation());

    if(matrixFirstTime)
    {
      time_t timer1,timer2;
      time(&timer1);

      G4Box* DABox((G4Box*)DALV->GetSolid()); 
      G4Box* voxelBox((G4Box*)motherParameterisation->ComputeSolid(0,motherPV));
      const int nxVoxel=round(DABox->GetXHalfLength()/voxelBox->GetXHalfLength());
      const int nyVoxel=round(DABox->GetYHalfLength()/voxelBox->GetYHalfLength());
      const int nzVoxel=round(DABox->GetZHalfLength()/voxelBox->GetZHalfLength());

      /*G4cout<<" * Voxels info :"<<G4endl;
      G4cout<<"Total x voxels="<<nxVoxel<<G4endl;
      G4cout<<"Total y Voxels="<<nyVoxel<<G4endl;
      G4cout<<"Total z Voxels="<<nzVoxel<<G4endl;*/

      voxelCubicVolume.resize(nxVoxel);
      voxelMass.resize(nxVoxel);
      for(int x=0;x<nxVoxel;x++)
      {
        voxelCubicVolume[x].resize(nyVoxel);
        voxelMass[x].resize(nyVoxel);
        for(int y=0;y<nyVoxel;y++)
        {
          voxelCubicVolume[x][y].resize(nzVoxel,-1.);
          voxelMass[x][y].resize(nzVoxel,-1.);
        }
      }
      doselReconstructedCubicVolume.resize(mDoseImage.GetValueImage().GetNumberOfValues(),-1.);
      doselReconstructedMass.resize(mDoseImage.GetValueImage().GetNumberOfValues(),-1.);

      //G4cout<<space<<" Computing data form parameterised volume, please wait."<<G4endl;
      matrixFirstTime=false;
      for(signed long int i=0;i<motherPV->GetMultiplicity();i++)
      {
        motherParameterisation->ComputeTransformation(i,motherPV);

        //Computing voxel x y z :
        int xVoxel(round((DABox->GetXHalfLength()+motherPV->GetTranslation().getX()-voxelBox->GetXHalfLength())/(voxelBox->GetXHalfLength()*2.0)));
        int yVoxel(round((DABox->GetYHalfLength()+motherPV->GetTranslation().getY()-voxelBox->GetYHalfLength())/(voxelBox->GetYHalfLength()*2.0)));
        int zVoxel(round((DABox->GetZHalfLength()+motherPV->GetTranslation().getZ()-voxelBox->GetZHalfLength())/(voxelBox->GetZHalfLength()*2.0)));

        /*G4cout<<"xVoxel="<<xVoxel<<G4endl;
        G4cout<<"yVoxel="<<yVoxel<<G4endl;
        G4cout<<"zVoxel="<<zVoxel<<G4endl;*/

        if(xVoxel>=nxVoxel||yVoxel>=nyVoxel||zVoxel>=nzVoxel)
        {
          GateError("!!! ERROR : Too many voxels !!! (xVoxel="<<xVoxel<<",yVoxel="<<yVoxel<<",zVoxel="<<zVoxel<<")"<<Gateendl);
          exit(EXIT_FAILURE);
        }

        voxelCubicVolume[xVoxel][yVoxel][zVoxel]=motherParameterisation->ComputeSolid(i,motherPV)->GetCubicVolume();
        voxelMass[xVoxel][yVoxel][zVoxel]=motherParameterisation->ComputeMaterial(i,motherPV)->GetDensity()*motherParameterisation->ComputeSolid(i,motherPV)->GetCubicVolume();
      }

      double doselReconstructedTotalMass(0.);
      double doselReconstructedTotalCubicVolume(0.);

      for(long int i=0;i<mDoseImage.GetValueImage().GetNumberOfValues();i++)
      {
        long double xDoselmin((DABox->GetXHalfLength()+mDoseImage.GetValueImage().GetVoxelCenterFromIndex(i).getX()-mDoseImage.GetValueImage().GetVoxelSize().getX()/2.0)/(voxelBox->GetXHalfLength()*2.0)),
                    xDoselmax((DABox->GetXHalfLength()+mDoseImage.GetValueImage().GetVoxelCenterFromIndex(i).getX()+mDoseImage.GetValueImage().GetVoxelSize().getX()/2.0)/(voxelBox->GetXHalfLength()*2.0)),
                    yDoselmin((DABox->GetYHalfLength()+mDoseImage.GetValueImage().GetVoxelCenterFromIndex(i).getY()-mDoseImage.GetValueImage().GetVoxelSize().getY()/2.0)/(voxelBox->GetYHalfLength()*2.0)),
                    yDoselmax((DABox->GetYHalfLength()+mDoseImage.GetValueImage().GetVoxelCenterFromIndex(i).getY()+mDoseImage.GetValueImage().GetVoxelSize().getY()/2.0)/(voxelBox->GetYHalfLength()*2.0)),
                    zDoselmin((DABox->GetZHalfLength()+mDoseImage.GetValueImage().GetVoxelCenterFromIndex(i).getZ()-mDoseImage.GetValueImage().GetVoxelSize().getZ()/2.0)/(voxelBox->GetZHalfLength()*2.0)),
                    zDoselmax((DABox->GetZHalfLength()+mDoseImage.GetValueImage().GetVoxelCenterFromIndex(i).getZ()+mDoseImage.GetValueImage().GetVoxelSize().getZ()/2.0)/(voxelBox->GetZHalfLength()*2.0));

        //G4cout<<" * Dosel info :"<<G4endl;
        //G4cout<<"xmin="<<round(xDoselmin)<<", xmax="<<round(xDoselmax)<<G4endl;
        //G4cout<<"ymin="<<round(yDoselmin)<<", ymax="<<round(yDoselmax)<<G4endl;
        //G4cout<<"zmin="<<round(zDoselmin)<<", zmax="<<round(zDoselmax)<<G4endl;

        std::vector<long double> doselMin(3,-1.),doselMax(3,-1.);
        doselMin[0]=xDoselmin; doselMax[0]=xDoselmax;
        doselMin[1]=yDoselmin; doselMax[1]=yDoselmax;
        doselMin[2]=zDoselmin; doselMax[2]=zDoselmax;

        for(int x=round(xDoselmin);x<round(xDoselmax);x++)
          for(int y=round(yDoselmin);y<round(yDoselmax);y++)
            for(int z=round(zDoselmin);z<round(zDoselmax);z++)
            {
              std::vector<bool> isMin(3,false),isMax(3,false);
              std::vector<int>    origCoord(3,-1);
              std::vector<std::vector<int> >    coord(3);
              std::vector<std::vector<double> > coef(3);

              origCoord[0]=x; origCoord[1]=y; origCoord[2]=z;

              // Dimension of the vectors : x = 0, y = 1, z = 2

              for(int dim=0;dim<3;dim++)
              {
                if(origCoord[dim]==round(doselMin[dim])&&fmod(doselMin[dim],1)>1e-8)
                  isMin[dim]=true;
                if(origCoord[dim]==round(doselMax[dim])-1&&fmod(doselMax[dim],1)>1e-8)
                  isMax[dim]=true;

                if(isMin[dim]&&isMax[dim])
                {
                  if(fmod(doselMin[dim],1)>=0.5&&fmod(doselMax[dim],1)<0.5)
                  {
                    coord[dim].push_back(origCoord[dim]-1);
                    coord[dim].push_back(origCoord[dim]);
                    coord[dim].push_back(origCoord[dim]+1);

                    coef[dim].push_back(1-fmod(doselMin[dim],1));
                    coef[dim].push_back(1);
                    coef[dim].push_back(fmod(doselMax[dim],1));
                  }
                  else if(fmod(doselMin[dim],1)>=0.5)
                  {
                    coord[dim].push_back(origCoord[dim]-1);
                    coord[dim].push_back(origCoord[dim]);

                    coef[dim].push_back(1-fmod(doselMin[dim],1));
                    coef[dim].push_back(fmod(doselMax[dim],1));
                  }
                  else if(fmod(doselMax[dim],1)<0.5)
                  {
                    coord[dim].push_back(origCoord[dim]);
                    coord[dim].push_back(origCoord[dim]+1);

                    coef[dim].push_back(1-fmod(doselMin[dim],1));
                    coef[dim].push_back(fmod(doselMax[dim],1));
                  }
                  else
                  {
                    coord[dim].push_back(origCoord[dim]);
                    coef[dim].push_back(abs((1-fmod(doselMin[dim],1))-fmod(doselMax[dim],1))); //FIXME ?
                  }
                }
                else if(isMin[dim])
                {
                  if(fmod(doselMin[dim],1)>=0.5)
                  {
                    coord[dim].push_back(origCoord[dim]-1);
                    coord[dim].push_back(origCoord[dim]);

                    coef[dim].push_back(1-fmod(doselMin[dim],1));
                    coef[dim].push_back(1);
                  }
                  else
                  {
                    coord[dim].push_back(origCoord[dim]);
                    coef[dim].push_back(1-fmod(doselMin[dim],1));
                  }
                }
                else if(isMax[dim])
                {
                  if(fmod(doselMax[dim],1)<0.5)
                  {
                    coord[dim].push_back(origCoord[dim]);
                    coord[dim].push_back(origCoord[dim]+1);

                    coef[dim].push_back(1);
                    coef[dim].push_back(fmod(doselMax[dim],1));
                  }
                  else
                  {
                    coord[dim].push_back(origCoord[dim]);
                    coef[dim].push_back(fmod(doselMax[dim],1));
                  }
                }
                else
                {
                  coord[dim].push_back(origCoord[dim]);
                  coef[dim].push_back(1);
                }
                if(coord[dim].size()!=coef[dim].size())
                  GateError("!!! ERROR : Size of coord and coef are not the same !"<<Gateendl);
              }



              for(size_t xVox=0;xVox<coord[0].size();xVox++)
                for(size_t yVox=0;yVox<coord[1].size();yVox++)
                  for(size_t zVox=0;zVox<coord[2].size();zVox++)
                  {
                    double coefVox(coef[0][xVox]*coef[1][yVox]*coef[2][zVox]);

                    doselReconstructedCubicVolume[i]+=voxelCubicVolume[coord[0][xVox]][coord[1][yVox]][coord[2][zVox]]*coefVox;
                    doselReconstructedMass[i]+=voxelMass[coord[0][xVox]][coord[1][yVox]][coord[2][zVox]]*coefVox;
                  }
            }
        doselReconstructedTotalMass+=doselReconstructedMass[i];
        doselReconstructedTotalCubicVolume+=doselReconstructedCubicVolume[i];

        //if(i!=0&&abs(doselReconstructedCubicVolume[i]-doselCubicVolume)>0.5)
        //  GateWarning("Dosel "<<i<<" has incorrect reconstructed cubic volume : "<<Gateendl
        //      <<"   "<<G4BestUnit(doselReconstructedCubicVolume[i],"Volume")<<" (Original cubic volume : "<<G4BestUnit(mDoseImage.GetVoxelVolume(),"Volume")<<")"<<Gateendl
        //      <<"   "<< "xmin="<<xDoselmin<<" ("<<round(xDoselmin)<<"), xmax="<<xDoselmax<<" ("<<round(xDoselmax)<<")"<<Gateendl
        //      <<"   "<< "ymin="<<yDoselmin<<" ("<<round(yDoselmin)<<"), ymax="<<yDoselmax<<" ("<<round(yDoselmax)<<")"<<Gateendl
        //      <<"   "<< "zmin="<<zDoselmin<<" ("<<round(zDoselmin)<<"), zmax="<<zDoselmax<<" ("<<round(zDoselmax)<<")"<<Gateendl);

      }
      time(&timer2);
      seconds=difftime(timer2,timer1);

      G4cout<<G4endl<<space<<"================================================================"<<G4endl;
      G4cout<<space<<" * SUMMARY : Dosel mass calculation for voxelized volume :"<<G4endl;
      G4cout<<space<<"     Time elapsed : "<<seconds<<" seconds"<<G4endl;
      G4cout<<space<<"     Number of voxels : "<<motherPV->GetMultiplicity()<<G4endl;
      G4cout<<space<<"     Number of dosels : "<<mDoseImage.GetValueImage().GetNumberOfValues()<<G4endl;
      G4cout<<space<<"     Dosels reconstructed total mass : "<<G4BestUnit(doselReconstructedTotalMass,"Mass")<<G4endl;
      G4cout<<space<<"     Dosels reconstructed total cubic volume : "<<G4BestUnit(doselReconstructedTotalCubicVolume,"Volume")<<G4endl;
      G4cout<<space<<"================================================================"<<G4endl<<G4endl;
    }

    motherProgenyMass=doselReconstructedMass[index];
    motherProgenyCubicVolume=doselCubicVolume;
  }
  //////////////////////////////////////////////////////////////////////////
  else
  {
    // Overlap Mother-Dosel
    motherSV=new G4IntersectionSolid(motherSV->GetName()+"∩"+doselSV.GetName(),
                                    motherSV,
                                    &doselSV, 
                                    &motherRotation, // Local rotation
                                    doselTranslation); // Local translation

    //double motherVoxelOverlapCubicVolume(RealZero(motherSV->GetCubicVolume(),motherLV->GetSolid()->GetCubicVolume(),0.00001));
    double motherDoselOverlapCubicVolume(motherSV->GetCubicVolume());

    G4cout<<G4endl<<space<<"================================================================"<<G4endl;
    G4cout<<space<<" * SUMMARY : Mother : "<<motherLV->GetName()<<" (generation n°"<<Generation<<") :"<<G4endl
          <<space<<"    Material           : "<<motherLV->GetMaterial()->GetName()<<G4endl
          <<space<<"    Density            : "<<G4BestUnit(motherDensity,"Volumic Mass")<<G4endl
          <<space<<"    Original volume    : "<<G4BestUnit(motherLV->GetSolid()->GetCubicVolume(),"Volume")<<G4endl
          <<space<<"    Overlap volume (with dosel "<<index<<") : "<<G4BestUnit(motherSV->GetCubicVolume(),"Volume")<<G4endl
          <<space<<"    Diff volume        : "<<(motherSV->GetCubicVolume()/motherLV->GetSolid()->GetCubicVolume())<<G4endl
          <<space<<"    Overlap voxel      : "<<G4BestUnit(motherDoselOverlapCubicVolume,"Volume")<<G4endl;
      G4cout<<space<<"================================================================"<<G4endl<<G4endl;

    // If the mother's intersects the voxel.
    //if(motherVoxelOverlapCubicVolume==0.)
    if(motherDoselOverlapCubicVolume==0.)
    {
      G4cout<<space<<" *** IS NOT IN THE DOSEL N°"<<index<<" ***"<<G4endl;
      daughterVolume=0.; // Only for testing
      if(Generation>0)
        space.resize(space.size()-3);
      return 0.;
    }

    // Calculation for daughter(s) ///////////////////////////////////////////
    if(motherLV->GetNoDaughters()>0) 
    {
      //G4cout<<space<<"    Has "<<motherLV->GetNoDaughters()<<" daughter(s)"<<G4endl;
      
      for(int i=0;i<motherLV->GetNoDaughters();i++)
      {
        G4VPhysicalVolume*  daughterPV(motherLV->GetDaughter(i));
        G4LogicalVolume*    daughterLV(daughterPV->GetLogicalVolume());
        G4VSolid*           daughterSV(daughterLV->GetSolid());

        /*GateDebugMessage("Actor", 5,space<<"   * Daughter n°"<<i<<" : "<<daughterLV->GetName()<<" :"<<Gateendl
            <<space<<"       Has "<<daughterLV->GetNoDaughters()<<" daughter(s)"<<Gateendl);*/

        // Calculation of absolute translation and rotation.
        G4RotationMatrix daughterRotation(daughterPV->GetObjectRotationValue());
        G4ThreeVector    daughterTranslation(daughterPV->GetObjectTranslation());

        if(Generation>0)
        {
          daughterTranslation=motherTranslation+daughterPV->GetObjectTranslation().transform(motherRotation);
          daughterRotation=daughterRotation.transform(motherRotation);
        }


        // Substraction Mother-Daughter
        motherSV=new G4SubtractionSolid(motherSV->GetName()+"-"+daughterSV->GetName(),
                                      motherSV, // Already overlapped with voxel volume
                                      daughterSV,
                                      daughterPV->GetObjectRotation(), // Local rotation
                                      daughterPV->GetObjectTranslation()); // Local translation

        //G4cout<<space<<"  Mother-daughter overlap volume     : "<<G4BestUnit(motherSV->GetCubicVolume(),"Volume")<<G4endl;

        // Calculation of the mass of the Mother Progeny
        motherProgenyMass+=VoxelIteration(daughterPV,Generation+1,daughterRotation,daughterTranslation,index);
        motherProgenyCubicVolume+=daughterVolume; // Only for verification.
      }
    }
    /*else
      G4cout<<space<<"    Has no daughter"<<G4endl;*/
    //////////////////////////////////////////////////////////////////////////

    // Mother mass & volume //////////////////////////////////////////////////
    double motherCubicVolume(motherSV->GetCubicVolume());
    motherMass=motherCubicVolume*motherDensity;

    motherProgenyMass+=motherMass;
    motherProgenyCubicVolume+=motherCubicVolume;
    //////////////////////////////////////////////////////////////////////////
  }

  // Mother information dump ///////////////////////////////////////////////
  /*G4cout<<space<<"  Mother original volume     : "<<G4BestUnit(motherLV->GetSolid()->GetCubicVolume(),"Volume")<<G4endl
        <<space<<"  Mother w/o daughters volume: "<<G4BestUnit(motherCubicVolume,"Volume")<<G4endl
        <<space<<"  Mother  +  daughters volume: "<<G4BestUnit(motherProgenyCubicVolume,"Volume")<<G4endl
        <<space<<"  Mother Mass                : "<<G4BestUnit(motherMass,"Mass")<<G4endl
        <<space<<"  Mother  +  daughters mass  : "<<G4BestUnit(motherProgenyMass,"Mass")<<" (MotherProgenyMass)"<<G4endl;*/
  //////////////////////////////////////////////////////////////////////////

  if(Generation>0)
    space.resize(space.size()-3);

  daughterVolume=motherProgenyCubicVolume;
  return motherProgenyMass;
}

bool GateDoseActor::HasDaughter(const G4LogicalVolume* doseactorLV)
{
    if(doseactorLV->GetNoDaughters()>0)
      return true;
    return false;
}
//-----------------------------------------------------------------------------

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



  if(FirstTime)
  {
    FirstTime=false;
    //bool mIsVoxelised(false);
    VoxelDensity=0.;
    VoxelVolume=0.;
    doselMass.resize(mDoseImage.GetValueImage().GetNumberOfValues(),-1.);

    DAPV=G4PhysicalVolumeStore::GetInstance()->GetVolume(mVolumeName+"_phys");
    DALV=G4LogicalVolumeStore::GetInstance()->GetVolume(mVolumeName+"_log");
    mHasDaughter=HasDaughter(DALV);

    VoxelIteration(DAPV,0,DAPV->GetObjectRotationValue(),DAPV->GetObjectTranslation(),index); //FIXME For testing only !!!!

    /*G4cout<<G4endl
          <<"DEBUG of DoseActor "<< GetObjectName()<<" :"<<G4endl
          <<" DA volume name   : "<<mVolumeName<<G4endl;
    if(mHasDaughter)
      G4cout<<" INFO : Has daughter(s)"<<G4endl;*/
    /*if(mHasDaughter&&!mIsVoxelised)
      voxelMass[0]=VolumeIteration(DAPV,0,DAPV->GetObjectRotationValue(),DAPV->GetObjectTranslation()).first;
    else if(mHasDaughter)*/

    /*time_t timer1,timer2;
    time(&timer1);

      double totalMass(0.);
      for(long int i=0;i<mDoseImage.GetValueImage().GetNumberOfValues();i++)
      {
        doselMass[i]=VoxelIteration(DAPV,0,DAPV->GetObjectRotationValue(),DAPV->GetObjectTranslation(),i);
        totalMass+=doselMass[i];*/


        //if(i%10==0) std::cout<<" * Please wait ... "<<i*100/mDoseImage.GetValueImage().GetNumberOfValues()<<"%\r"<<std::flush;
        //G4cout<<" DEBUG : VoxelNx="<<mDoseImage.GetValueImage().GetVoxelNx()<<", VoxelNy="<<mDoseImage.GetValueImage().GetVoxelNy()<<", VoxelNz="<<mDoseImage.GetValueImage().GetVoxelNz()<<G4endl;
        /*G4cout<<" DEBUG : Index : "<<i<<G4endl;
        G4cout<<" DEBUG : Voxel size   : X="<<G4BestUnit(mDoseImage.GetValueImage().GetVoxelSize().getX(),"Length")<<", Y="<<G4BestUnit(mDoseImage.GetValueImage().GetVoxelSize().getY(),"Length")<<", Z="<<G4BestUnit(mDoseImage.GetValueImage().GetVoxelSize().getZ(),"Length")<<G4endl;
        G4cout<<" DEBUG : Voxel center : X="<<G4BestUnit(mDoseImage.GetValueImage().GetVoxelCenterFromIndex(i).getX(),"Length")<<", Y="<<G4BestUnit(mDoseImage.GetValueImage().GetVoxelCenterFromIndex(i).getY(),"Length")<<", Z="<<G4BestUnit(mDoseImage.GetValueImage().GetVoxelCenterFromIndex(i).getZ(),"Length")<<G4endl;
        G4cout<<" DEBUG : Voxel corner : X="<<G4BestUnit(mDoseImage.GetValueImage().GetVoxelCornerFromIndex(i).getX(),"Length")<<", Y="<<G4BestUnit(mDoseImage.GetValueImage().GetVoxelCornerFromIndex(i).getY(),"Length")<<", Z="<<G4BestUnit(mDoseImage.GetValueImage().GetVoxelCornerFromIndex(i).getZ(),"Length")<<G4endl;
        //G4cout<<" DEBUG : doselMass["<<i<<"] : "<<G4BestUnit(doselMass[i],"Mass")<<G4endl;
      }*/

    /*time(&timer2);
    seconds=difftime(timer2,timer1);

    G4cout<<G4endl;
    G4cout<<" DEBUG : Dosels totalMass : "<<G4BestUnit(totalMass,"Mass")<<G4endl;
    G4cout<<" DEBUG : Total time : "<<seconds<<" s"<<G4endl;

    GateDebugMessage("Actor", 5," DA volume           : "<<G4BestUnit(DALVCubicVolume,"Volume")<<Gateendl
          <<" DA original mass    : "<<G4BestUnit(DALV->GetMass(true,true,0),"Mass")<<" (propagated)"<<Gateendl
          <<" DA original mass    : "<<G4BestUnit(DALV->GetMass(false,false,0),"Mass")<<" (not propagated)"<<Gateendl
          <<" DA original density : "<<G4BestUnit(DALV->GetMass()/DALVCubicVolume,"Volumic Mass")<<Gateendl
          <<" DA corrected mass   : "<<G4BestUnit(VoxelMass,"Mass")<<Gateendl
          <<" DA corrected density: "<<G4BestUnit(VoxelDensity,"Volumic Mass")<<Gateendl);*/
  }

  double dose=0.;
  double density = step->GetPreStepPoint()->GetMaterial()->GetDensity();

  /*if(mHasDaughter)
  {
    if(doselMass[index]==-1.)
    {
      nbofRecVoxels++;
      doselMass[index]=VoxelIteration(DAPV,0,DAPV->GetObjectRotationValue(),DAPV->GetObjectTranslation(),index);
      //if(nbofRecVoxels%100==0) std::cout<<" * Please wait ... "<<(double)nbofRecVoxels*100/(double)mDoseImage.GetValueImage().GetNumberOfValues()<<"% of dosels reconstructed ...\r"<<std::flush;
    }
    density = doselMass[index]/mDoseImage.GetVoxelVolume();
  }*/


  if (mIsDoseImageEnabled) {

    // ------------------------------------
    // Convert deposited energy into Gray

    // OLD version (correct but not clear)
    // dose = edep/density*1e12/mDoseImage.GetVoxelVolume();

    // NEW version (same results but more clear)
    dose = edep/density/mDoseImage.GetVoxelVolume()/gray;
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

  //if (mIsMassImageEnabled) mMassImage.AddValue(index, doselMass[index]);
  if (mMassFirstTime&&mIsMassImageEnabled)// FIXME : For testing only!!!!!
  {
    mMassFirstTime=false;
    for(size_t i=0;i<doselReconstructedMass.size();i++)
      mMassImage.AddValue(i, doselReconstructedMass[i]/6.24151e+21); // FIXME : For testing only!!!!!
  }

  GateDebugMessageDec("Actor", 4, "GateDoseActor -- UserSteppingActionInVoxel -- end\n");
}
//-----------------------------------------------------------------------------
