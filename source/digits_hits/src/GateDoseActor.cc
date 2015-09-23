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
#include <G4SmartVoxelHeader.hh>
#include <G4VPVParameterisation.hh>

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
  matrixFirstTime=true;
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

  G4Box doselSV("Voxel_index",//FIXME convertir i en string
                mDoseImage.GetValueImage().GetVoxelSize().getX()/2.0,
                mDoseImage.GetValueImage().GetVoxelSize().getY()/2.0,
                mDoseImage.GetValueImage().GetVoxelSize().getZ()/2.0);

  // Calculation of voxel's local rotation and translation
  G4RotationMatrix doselRotation(mDoseImage.GetValueImage().GetTransformMatrix());
  G4ThreeVector    doselTranslation(mDoseImage.GetValueImage().GetVoxelCenterFromIndex(index));

  if(Generation>0)
    doselTranslation-=motherTranslation;

  // Parameterisation //////////////////////////////////////////////////////
  //G4cout<<space<<"       Parameterised ?     : "<<motherPV->IsParameterised()<<G4endl;
  //G4cout<<space<<"       Replicated ?        : "<<motherPV->IsReplicated()<<G4endl;
  double doselCubicVolume(doselSV.GetCubicVolume());
  double doselMass(0.);
  if(motherPV->IsParameterised())
  {
    G4cout<<space<<"       motherPV->GetName()         : "<<motherPV->GetName() <<G4endl;
    G4cout<<space<<"       motherPV->GetMultiplicity() : "<<motherPV->GetMultiplicity() <<G4endl;
    G4VPVParameterisation* motherParameterisation(motherPV->GetParameterisation());
    double doselReconstructedCubicVolume(0.); // For control only
    int nbVoxelInsideDosel(0);

    if(matrixFirstTime)
    {
      G4cout<<space<<" Computing data form parameterised volume, please wait."<<G4endl;
      matrixFirstTime=false;
      voxelMatrix.resize(motherPV->GetMultiplicity(),0.);
      voxelAbsoluteTranslation.resize(motherPV->GetMultiplicity());
      voxelCubicVolume.resize(motherPV->GetMultiplicity(),0.);
      voxelMass.resize(motherPV->GetMultiplicity(),0.);
      for(signed long int i=0;i<motherPV->GetMultiplicity();i++)
      {
        voxelMatrix[i]=i;
        motherParameterisation->ComputeTransformation(i,motherPV);
        voxelAbsoluteTranslation[i]=motherPV->GetTranslation();
        voxelCubicVolume[i]=motherParameterisation->ComputeSolid(i,motherPV)->GetCubicVolume();
        voxelMass[i]=motherParameterisation->ComputeMaterial(i,motherPV)->GetDensity()*voxelCubicVolume[i];
      }
    }

    G4cout<<G4endl<<space<<"* Computing the parameterised volume for dosel n°"<<index<<", please wait."<<G4endl;

    //for(long int i=0;i<motherPV->GetMultiplicity();i++)
    for(long int i=0;i<motherPV->GetMultiplicity();i++)
      if(voxelMatrix[i]!=-1)
      {
        //G4cout<<space<<"       voxelMatrix["<<i<<"] : "<<voxelMatrix[i]<<G4endl;
        G4ThreeVector voxelTranslation=voxelAbsoluteTranslation[i]-mDoseImage.GetValueImage().GetVoxelCenterFromIndex(index);

        //G4cout<<space<<"         Dosel->Voxel Translation : X="<<G4BestUnit(voxelTranslation.getX(),"Length")<<",Y="<<G4BestUnit(voxelTranslation.getY(),"Length")<<",Z="<<G4BestUnit(voxelTranslation.getZ(),"Length")<<G4endl;
        //G4cout<<space<<"         Mother->Voxel Translation : X="<<G4BestUnit(motherPV->GetTranslation().getX(),"Length")<<",Y="<<G4BestUnit(motherPV->GetTranslation().getY(),"Length")<<",Z="<<G4BestUnit(motherPV->GetTranslation().getZ(),"Length")<<G4endl;

        if(doselSV.Inside(voxelTranslation)==kInside) //FIXME : This method takes way too much time.
        {
          nbVoxelInsideDosel++;

          //if(i%100000==0) std::cout<<i*100/voxelMatrix.size()<<"%\r"<<std::flush;

          doselReconstructedCubicVolume+=voxelCubicVolume[i];
          doselMass+=voxelMass[i];

          voxelMatrix[i]=-1;
        }

        /*G4cout<<space<<"       Voxel n°"<<i<<" :"<<G4endl;
        G4cout<<space<<"         CubicVolume : "<<G4BestUnit(motherParameterisation->ComputeSolid(i,motherPV)->GetCubicVolume(),"Volume")<<G4endl;
        G4cout<<space<<"         Material    : "<<motherParameterisation->ComputeMaterial(i,motherPV)->GetName()<<G4endl;
        G4cout<<space<<"         Translation : X="<<G4BestUnit(motherPV->GetTranslation().getX(),"Length")<<",Y="<<G4BestUnit(motherPV->GetTranslation().getY(),"Length")<<",Z="<<G4BestUnit(motherPV->GetTranslation().getZ(),"Length")<<G4endl;*/

      }

    G4cout<<G4endl;
    G4cout<<space<<"       Nb of voxels inside dosel : "<<nbVoxelInsideDosel<<"/"<<motherPV->GetMultiplicity()<<" ("<<nbVoxelInsideDosel*100/motherPV->GetMultiplicity()<<"%)"<<G4endl;
    G4cout<<space<<"       Dosel original cubic volume : "<<G4BestUnit(doselCubicVolume,"Volume")<<G4endl;
    G4cout<<space<<"       Dosel reconstructed cubic volume : "<<G4BestUnit(doselReconstructedCubicVolume,"Volume")<<G4endl;
    G4cout<<space<<"       Dosel mass : "<<G4BestUnit(doselMass,"Mass")<<G4endl;
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

    double motherVoxelOverlapCubicVolume(RealZero(motherSV->GetCubicVolume(),motherLV->GetSolid()->GetCubicVolume(),0.00001));

    G4cout<<space<<"* Mother : "<<motherLV->GetName()<<" (generation n°"<<Generation<<") :"<<G4endl
          <<space<<"    Original volume    : "<<G4BestUnit(motherLV->GetSolid()->GetCubicVolume(),"Volume")<<G4endl
          <<space<<"    Overlap  volume    : "<<G4BestUnit(motherSV->GetCubicVolume(),"Volume")<<G4endl
          <<space<<"    Diff volume        : "<<(motherSV->GetCubicVolume()/motherLV->GetSolid()->GetCubicVolume())<<G4endl
          <<space<<"    Overlap voxel      : "<<G4BestUnit(motherVoxelOverlapCubicVolume,"Volume")<<G4endl
          <<space<<"    Density            : "<<G4BestUnit(motherDensity,"Volumic Mass")<<G4endl;

    // If the mother's intersects the voxel.
    //if(motherVoxelOverlapCubicVolume==0.)
    if(motherVoxelOverlapCubicVolume==0.)
    {
      G4cout<<space<<" *** IS NOT IN THE DOXEL N°"<<index<<" ***"<<G4endl;
      daughterVolume=0.; // Only for testing
      if(Generation>0)
        space.resize(space.size()-3);
      return 0.;
    }
  }

  // Voxel information dump ////////////////////////////////////////////////
  /*space+="   ";
  G4cout<<space<<"* Voxel n°"<<index<<" informations :"<<G4endl
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

  //////////////////////////////////////////////////////////////////////////
  if(motherLV->GetNoDaughters()>0) // If the mother's volume has daughter(s).
  {
    G4cout<<space<<"    Has "<<motherLV->GetNoDaughters()<<" daughter(s)"<<G4endl;
    
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

      G4cout<<space<<"  Mother-daughter overlap volume     : "<<G4BestUnit(motherSV->GetCubicVolume(),"Volume")<<G4endl;

      // Calculation of the mass of the Mother Progeny
      motherProgenyMass+=VoxelIteration(daughterPV,Generation+1,daughterRotation,daughterTranslation,index);
      motherProgenyCubicVolume+=daughterVolume; // Only for verification.
    }
  }
  else
    G4cout<<space<<"    Has no daughter"<<G4endl;
  //////////////////////////////////////////////////////////////////////////


  // Mother mass & volume //////////////////////////////////////////////////
  double motherCubicVolume(motherSV->GetCubicVolume());
  if(motherPV->IsParameterised())
  {
    motherProgenyMass=doselMass;
    motherProgenyCubicVolume=doselCubicVolume;
  }
  else
  {
    motherMass=motherCubicVolume*motherDensity;

    motherProgenyMass+=motherMass;
    motherProgenyCubicVolume+=motherCubicVolume;
  }
  //////////////////////////////////////////////////////////////////////////

  // Mother information dump ///////////////////////////////////////////////
  G4cout<<space<<"  Mother original volume     : "<<G4BestUnit(motherLV->GetSolid()->GetCubicVolume(),"Volume")<<G4endl
        <<space<<"  Mother w/o daughters volume: "<<G4BestUnit(motherCubicVolume,"Volume")<<G4endl
        <<space<<"  Mother  +  daughters volume: "<<G4BestUnit(motherProgenyCubicVolume,"Volume")<<G4endl
        <<space<<"  Mother Mass                : "<<G4BestUnit(motherMass,"Mass")<<G4endl
        <<space<<"  Mother  +  daughters mass  : "<<G4BestUnit(motherProgenyMass,"Mass")<<" (MotherProgenyMass)"<<G4endl;
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

    //FIXME : Savoir si le volume associé au DoseActor a des filles. Si oui utiliser l'algo, si non utiliser la formule de la densité originelle. Faire un booléen.
    //

    FirstTime=false;
    bool mIsVoxelised(false);
    VoxelDensity=0.;
    VoxelVolume=0.;
    doselMass.resize(mDoseImage.GetValueImage().GetNumberOfValues(),0.);

    DAPV=G4PhysicalVolumeStore::GetInstance()->GetVolume(mVolumeName+"_phys");
    DALV=G4LogicalVolumeStore::GetInstance()->GetVolume(mVolumeName+"_log");

    // Study of voxelised volumes:
    //G4cout<<"DALV->GetVoxelHeader()->GetNoSlices():"<<DALV->GetVoxelHeader()->GetNoSlices()<<G4endl; // Doesn't work (invalid use of incomplete type ‘struct G4SmartVoxelHeader’)

    mHasDaughter=HasDaughter(DALV);

    G4cout<<G4endl
          <<"DEBUG of DoseActor "<< GetObjectName()<<" :"<<G4endl
          <<" DA volume name   : "<<mVolumeName<<G4endl;
    if(mHasDaughter)
      G4cout<<" INFO : Has daughter(s)"<<G4endl;

    if(mDoseImage.GetValueImage().GetNumberOfValues()>1&&DALV->GetSolid()->GetCubicVolume()>mDoseImage.GetVoxelVolume())
      mIsVoxelised=true;

    if(mIsVoxelised)
      G4cout<<" INFO : Is voxelised"<<G4endl
            <<"        Number of voxels : "<<mDoseImage.GetValueImage().GetNumberOfValues()<<G4endl;

    G4cout        <<"   DALV->GetSolid()->GetCubicVolume() : "<<DALV->GetSolid()->GetCubicVolume()<<G4endl;
    G4cout        <<"   mDoseImage.GetVoxelVolume() : "<<mDoseImage.GetVoxelVolume()<<G4endl;


    /*if(mHasDaughter&&!mIsVoxelised)
      voxelMass[0]=VolumeIteration(DAPV,0,DAPV->GetObjectRotationValue(),DAPV->GetObjectTranslation()).first;
    else if(mHasDaughter)*/
      double totalMass(0.);
      for(long int i=0;i<mDoseImage.GetValueImage().GetNumberOfValues();i++)
      {
        if(i%10==0) std::cout<<" * Please wait ... "<<i<<"/"<<mDoseImage.GetValueImage().GetNumberOfValues()<<"\r"<<std::flush;
        //G4cout<<" DEBUG : VoxelNx="<<mDoseImage.GetValueImage().GetVoxelNx()<<", VoxelNy="<<mDoseImage.GetValueImage().GetVoxelNy()<<", VoxelNz="<<mDoseImage.GetValueImage().GetVoxelNz()<<G4endl;
        /*G4cout<<" DEBUG : Index : "<<i<<G4endl;
        G4cout<<" DEBUG : Voxel size   : X="<<G4BestUnit(mDoseImage.GetValueImage().GetVoxelSize().getX(),"Length")<<", Y="<<G4BestUnit(mDoseImage.GetValueImage().GetVoxelSize().getY(),"Length")<<", Z="<<G4BestUnit(mDoseImage.GetValueImage().GetVoxelSize().getZ(),"Length")<<G4endl;
        G4cout<<" DEBUG : Voxel center : X="<<G4BestUnit(mDoseImage.GetValueImage().GetVoxelCenterFromIndex(i).getX(),"Length")<<", Y="<<G4BestUnit(mDoseImage.GetValueImage().GetVoxelCenterFromIndex(i).getY(),"Length")<<", Z="<<G4BestUnit(mDoseImage.GetValueImage().GetVoxelCenterFromIndex(i).getZ(),"Length")<<G4endl;
        G4cout<<" DEBUG : Voxel corner : X="<<G4BestUnit(mDoseImage.GetValueImage().GetVoxelCornerFromIndex(i).getX(),"Length")<<", Y="<<G4BestUnit(mDoseImage.GetValueImage().GetVoxelCornerFromIndex(i).getY(),"Length")<<", Z="<<G4BestUnit(mDoseImage.GetValueImage().GetVoxelCornerFromIndex(i).getZ(),"Length")<<G4endl;*/

        doselMass[i]=VoxelIteration(DAPV,0,DAPV->GetObjectRotationValue(),DAPV->GetObjectTranslation(),i);

        //FIXME : La densité du voxel ne varie pas selon la position de celui-ci dans le volume.



        G4cout<<" DEBUG : doselMass["<<i<<"] : "<<G4BestUnit(doselMass[i],"Mass")<<G4endl;
        totalMass+=doselMass[i];
      }
      G4cout<<" DEBUG : Dosels totalMass : "<<G4BestUnit(totalMass,"Mass")<<G4endl;


    GateDebugMessage("Actor", 5," DA volume           : "<<G4BestUnit(DALVCubicVolume,"Volume")<<Gateendl
          <<" DA original mass    : "<<G4BestUnit(DALV->GetMass(true,true,0),"Mass")<<" (propagated)"<<Gateendl
          <<" DA original mass    : "<<G4BestUnit(DALV->GetMass(false,false,0),"Mass")<<" (not propagated)"<<Gateendl
          <<" DA original density : "<<G4BestUnit(DALV->GetMass()/DALVCubicVolume,"Volumic Mass")<<Gateendl
          <<" DA corrected mass   : "<<G4BestUnit(VoxelMass,"Mass")<<Gateendl
          <<" DA corrected density: "<<G4BestUnit(VoxelDensity,"Volumic Mass")<<Gateendl);
  }

  double dose=0.;
  double density;
  if(mHasDaughter)
    density = doselMass[index]/mDoseImage.GetVoxelVolume();
  else
    density = step->GetPreStepPoint()->GetMaterial()->GetDensity();


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

  GateDebugMessageDec("Actor", 4, "GateDoseActor -- UserSteppingActionInVoxel -- end\n");
}
//-----------------------------------------------------------------------------
