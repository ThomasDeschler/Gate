/*----------------------
  Copyright (C): OpenGATE Collaboration

  This software is distributed under the terms
  of the GNU Lesser General  Public Licence (LGPL)
  See GATE/LICENSE.txt for further details
  ----------------------*/

/*
  \brief Class GateVoxelizedMass :
  \brief
*/

#include "GateVoxelizedMass.hh"
#include "GateMiscFunctions.hh"

#include <G4IntersectionSolid.hh>
#include <G4UnionSolid.hh>
#include <G4SubtractionSolid.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4VSolid.hh>
#include <G4Box.hh>
#include <G4VPVParameterisation.hh>
#include <time.h>

//-----------------------------------------------------------------------------
GateVoxelizedMass::GateVoxelizedMass()
{
  space="";
}
//-----------------------------------------------------------------------------

void GateVoxelizedMass::Initialize(const G4String mExtVolumeName, const GateImageDouble mExtImage)
{
  mVolumeName=mExtVolumeName;
  mImage=mExtImage;
  //G4cout<<"mImage.GetNumberOfValues()="<<mImage.GetNumberOfValues()<<G4endl;
  /*G4cout<<"mImage.GetVoxelNx()="<<mImage.GetVoxelNx()<<G4endl
          <<"mImage.GetVoxelNy()="<<mImage.GetVoxelNy()<<G4endl
          <<"mImage.GetVoxelNz()="<<mImage.GetVoxelNz()<<G4endl;*/
  mIsParameterised=false;
  mIsVecGenerated=false;

  DAPV=G4PhysicalVolumeStore::GetInstance()->GetVolume(mVolumeName+"_phys");
  DALV=DAPV->GetLogicalVolume();

  doselReconstructedMass.clear();
  doselReconstructedCubicVolume.clear();
  doselReconstructedMass.resize(mImage.GetNumberOfValues(),-1.);
  doselReconstructedCubicVolume.resize(mImage.GetNumberOfValues(),-1.);

  doselSV=new G4Box("DoselSV",
                    mImage.GetVoxelSize().getX()/2.0,
                    mImage.GetVoxelSize().getY()/2.0,
                    mImage.GetVoxelSize().getZ()/2.0);

  if(DALV->GetNoDaughters()==1&&DALV->GetDaughter(0)->IsParameterised())
  {
    mIsParameterised=true;
    GateVoxelizedMass::GenerateVoxels();
  }
}

//-----------------------------------------------------------------------------
double GateVoxelizedMass::GetVoxelMass(const int index)
{
  if(doselReconstructedMass[index]==-1.)
  {
    if(mIsParameterised)
    {
      GateVoxelizedMass::ParameterizedVolume(index);
    }
    else
    {
      doselReconstructedData=VoxelIteration(DAPV,0,DAPV->GetObjectRotationValue(),DAPV->GetObjectTranslation(),index);
      doselReconstructedMass[index]=doselReconstructedData.first;
      doselReconstructedCubicVolume[index]=doselReconstructedData.second;
    }
  }

  return doselReconstructedMass[index];
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
std::vector<double> GateVoxelizedMass::GetVoxelMassVector()
{
  GateVoxelizedMass::GenerateVectors();

  return doselReconstructedMass;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void GateVoxelizedMass::GenerateVectors()
{
  time_t timer1,timer2;
  time(&timer1);

  G4cout<<G4endl<<"================================================================"<<G4endl;
  G4cout<<" * Voxelized mass calculation in progress, please wait ... "<<G4endl;


  for(long int i=0;i<mImage.GetNumberOfValues();i++)
  {
    //std::cout<<round((double)i*100./(double)mImage.GetNumberOfValues())<<"%\r";

    if(mIsParameterised)
      GateVoxelizedMass::ParameterizedVolume(i);
    else
    {
      doselReconstructedData=VoxelIteration(DAPV,0,DAPV->GetObjectRotationValue(),DAPV->GetObjectTranslation(),i);
      doselReconstructedMass[i]=doselReconstructedData.first;
      doselReconstructedCubicVolume[i]=doselReconstructedData.second;
    }
  }

  doselReconstructedTotalCubicVolume=0.;
  doselReconstructedTotalMass=0.;
  for(size_t i=0;i<doselReconstructedMass.size();i++)
  {
    doselReconstructedTotalMass+=doselReconstructedMass[i];
    doselReconstructedTotalCubicVolume+=doselReconstructedCubicVolume[i];
  }

  time(&timer2);
  seconds=difftime(timer2,timer1);

  G4cout<<" * SUMMARY : Mass calculation for voxelized volume :"<<G4endl;
  G4cout<<"     Time elapsed : "<<seconds/60<<"min"<<seconds%60<<"s ("<<seconds<<"s)"<<G4endl;
  if(mIsParameterised)
    G4cout<<"     Number of voxels : "<<DALV->GetDaughter(0)->GetMultiplicity()<<G4endl;
  G4cout<<"     Number of dosels : "<<mImage.GetNumberOfValues()<<G4endl;
  G4cout<<"     Dosels reconstructed total mass : "<<G4BestUnit(doselReconstructedTotalMass,"Mass")<<G4endl;
  G4cout<<space<<"     Dosels reconstructed total cubic volume : "<<G4BestUnit(doselReconstructedTotalCubicVolume,"Volume")<<G4endl;
  G4cout<<"================================================================"<<G4endl<<G4endl;

  mIsVecGenerated=true;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void GateVoxelizedMass::GenerateVoxels()
{
  G4VPhysicalVolume*     daughterPV(DALV->GetDaughter(0));
  G4VPVParameterisation* daughterParameterisation(daughterPV->GetParameterisation());

  // Processing the voxels in the parameterized volume /////
  DABox=(G4Box*)DALV->GetSolid();
  voxelBox=(G4Box*)daughterParameterisation->ComputeSolid(0,daughterPV);

  const int nxVoxel=round(DABox->GetXHalfLength()/voxelBox->GetXHalfLength()),
            nyVoxel=round(DABox->GetYHalfLength()/voxelBox->GetYHalfLength()),
            nzVoxel=round(DABox->GetZHalfLength()/voxelBox->GetZHalfLength());

  //G4cout<<" * Voxels info :"<<G4endl;
  //G4cout<<"Total x voxels="<<nxVoxel<<G4endl;
  //G4cout<<"Total y Voxels="<<nyVoxel<<G4endl;
  //G4cout<<"Total z Voxels="<<nzVoxel<<G4endl;

  const int nxDosel=round(DABox->GetXHalfLength()/(mImage.GetVoxelSize().getX()/2.)),
            nyDosel=round(DABox->GetYHalfLength()/(mImage.GetVoxelSize().getY()/2.)),
            nzDosel=round(DABox->GetZHalfLength()/(mImage.GetVoxelSize().getZ()/2.));

  //G4cout<<" * Dosels info :"<<G4endl;
  //G4cout<<"Total x dosels="<<nxDosel<<G4endl;
  //G4cout<<"Total y dosels="<<nyDosel<<G4endl;
  //G4cout<<"Total z dosels="<<nzDosel<<G4endl;

  if(nxDosel>nxVoxel||nyDosel>nyVoxel||nzDosel>nzVoxel)
      GateError("!!! ERROR : The dosel resolution is smaller than the voxel resolution !!!"<<Gateendl);

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

  for(signed long int i=0;i<daughterPV->GetMultiplicity();i++) // Loop over the voxels
  {
    daughterParameterisation->ComputeTransformation(i,daughterPV);

    //Computing voxel x y z :
    const int xVoxel(round((DABox->GetXHalfLength()+daughterPV->GetTranslation().getX()-voxelBox->GetXHalfLength())/(voxelBox->GetXHalfLength()*2.0))),
              yVoxel(round((DABox->GetYHalfLength()+daughterPV->GetTranslation().getY()-voxelBox->GetYHalfLength())/(voxelBox->GetYHalfLength()*2.0))),
              zVoxel(round((DABox->GetZHalfLength()+daughterPV->GetTranslation().getZ()-voxelBox->GetZHalfLength())/(voxelBox->GetZHalfLength()*2.0)));

    //G4cout<<"xVoxel="<<xVoxel<<G4endl;
    //G4cout<<"yVoxel="<<yVoxel<<G4endl;
    //G4cout<<"zVoxel="<<zVoxel<<G4endl;

    if(xVoxel>=nxVoxel||yVoxel>=nyVoxel||zVoxel>=nzVoxel)
      GateError("!!! ERROR : Too many voxels !!! (xVoxel="<<xVoxel<<",yVoxel="<<yVoxel<<",zVoxel="<<zVoxel<<")"<<Gateendl);

    voxelCubicVolume[xVoxel][yVoxel][zVoxel]=daughterParameterisation->ComputeSolid(i,daughterPV)->GetCubicVolume();
    voxelMass[xVoxel][yVoxel][zVoxel]=daughterParameterisation->ComputeMaterial(i,daughterPV)->GetDensity()*daughterParameterisation->ComputeSolid(i,daughterPV)->GetCubicVolume();
  }
  /////////////////////////////////////////////////////////////////////////

  doselReconstructedCubicVolume.resize(mImage.GetNumberOfValues(),-1.);
  doselReconstructedMass.resize(mImage.GetNumberOfValues(),-1.);
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void GateVoxelizedMass::ParameterizedVolume(const int index)
{
  long double xDoselmin((DABox->GetXHalfLength()+mImage.GetVoxelCenterFromIndex(index).getX()-mImage.GetVoxelSize().getX()/2.0)/(voxelBox->GetXHalfLength()*2.0)),
              xDoselmax((DABox->GetXHalfLength()+mImage.GetVoxelCenterFromIndex(index).getX()+mImage.GetVoxelSize().getX()/2.0)/(voxelBox->GetXHalfLength()*2.0)),
              yDoselmin((DABox->GetYHalfLength()+mImage.GetVoxelCenterFromIndex(index).getY()-mImage.GetVoxelSize().getY()/2.0)/(voxelBox->GetYHalfLength()*2.0)),
              yDoselmax((DABox->GetYHalfLength()+mImage.GetVoxelCenterFromIndex(index).getY()+mImage.GetVoxelSize().getY()/2.0)/(voxelBox->GetYHalfLength()*2.0)),
              zDoselmin((DABox->GetZHalfLength()+mImage.GetVoxelCenterFromIndex(index).getZ()-mImage.GetVoxelSize().getZ()/2.0)/(voxelBox->GetZHalfLength()*2.0)),
              zDoselmax((DABox->GetZHalfLength()+mImage.GetVoxelCenterFromIndex(index).getZ()+mImage.GetVoxelSize().getZ()/2.0)/(voxelBox->GetZHalfLength()*2.0));

  //G4cout<<" * Dosel info :"<<G4endl;
  //G4cout<<"xmin="<<round(xDoselmin)<<", xmax="<<round(xDoselmax)<<G4endl;
  //G4cout<<"ymin="<<round(yDoselmin)<<", ymax="<<round(yDoselmax)<<G4endl;
  //G4cout<<"zmin="<<round(zDoselmin)<<", zmax="<<round(zDoselmax)<<G4endl;

  std::vector<double> doselMin(3,-1.),doselMax(3,-1.);
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

              doselReconstructedCubicVolume[index]+=voxelCubicVolume[coord[0][xVox]][coord[1][yVox]][coord[2][zVox]]*coefVox;
              doselReconstructedMass[index]+=voxelMass[coord[0][xVox]][coord[1][yVox]][coord[2][zVox]]*coefVox;
            }
      }
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
std::pair<double,double> GateVoxelizedMass::VoxelIteration(G4VPhysicalVolume* motherPV,const int Generation,G4RotationMatrix motherRotation,G4ThreeVector motherTranslation,const int index)
{
  //FIXME : Doesn't work with daughter overlapping its mother.

  if(motherPV->IsParameterised())
    GateError("The volume "<<motherPV->GetName()<<" is parameterized, can't compute its mass !"<<Gateendl);
  
  if(Generation>0)
    space+="   ";

  G4LogicalVolume* motherLV(motherPV->GetLogicalVolume());
  double motherMass(0.);
  double motherDensity(motherLV->GetMaterial()->GetDensity());
  double motherProgenyMass(0.);
  double motherProgenyCubicVolume(0.);
  G4VSolid* motherSV(motherLV->GetSolid());



  // Calculation of dosel's local rotation and translation
  G4RotationMatrix doselRotation(mImage.GetTransformMatrix());
  G4ThreeVector    doselTranslation(mImage.GetVoxelCenterFromIndex(index));

  if(Generation>0)
    doselTranslation-=motherTranslation;


  // Overlap Mother-Dosel
  motherSV=new G4IntersectionSolid(motherSV->GetName()+"∩"+doselSV->GetName(),
                                  motherSV,
                                  doselSV, 
                                  &motherRotation, // Local rotation
                                  doselTranslation); // Local translation

  double motherDoselOverlapCubicVolume(motherSV->GetCubicVolume());

  /*G4cout<<G4endl<<space<<"================================================================"<<G4endl;
  G4cout<<space<<" * SUMMARY : Mother : "<<motherLV->GetName()<<" (generation n°"<<Generation<<") :"<<G4endl
        <<space<<"    Material           : "<<motherLV->GetMaterial()->GetName()<<G4endl
        <<space<<"    Density            : "<<G4BestUnit(motherDensity,"Volumic Mass")<<G4endl
        <<space<<"    Original volume    : "<<G4BestUnit(motherLV->GetSolid()->GetCubicVolume(),"Volume")<<G4endl
        <<space<<"    Overlap volume (with dosel "<<index<<") : "<<G4BestUnit(motherSV->GetCubicVolume(),"Volume")<<G4endl
        <<space<<"    Diff volume        : "<<(motherSV->GetCubicVolume()/motherLV->GetSolid()->GetCubicVolume())<<G4endl
        <<space<<"    Overlap voxel      : "<<G4BestUnit(motherDoselOverlapCubicVolume,"Volume")<<G4endl;
    G4cout<<space<<"================================================================"<<G4endl<<G4endl;*/

  // Dosel information dump ////////////////////////////////////////////////
  /*space+="   ";
  G4cout<<space<<"* Doxel n°"<<index<<" informations :"<<G4endl
        <<space<<" --> Translation :"<<G4endl
        <<space<<"       Relative : X="<<G4BestUnit(doselTranslation.getX(),"Length")
                                <<",Y="<<G4BestUnit(doselTranslation.getY(),"Length")
                                <<",Z="<<G4BestUnit(doselTranslation.getZ(),"Length")<<G4endl
        <<space<<"       Absolute : X="<<G4BestUnit(mImage.GetVoxelCenterFromIndex(index).getX(),"Length")
                                <<",Y="<<G4BestUnit(mImage.GetVoxelCenterFromIndex(index).getY(),"Length")
                                <<",Z="<<G4BestUnit(mImage.GetVoxelCenterFromIndex(index).getZ(),"Length")<<G4endl
        <<space<<" --> Rotation :"<<G4endl
        <<space<<"       Relative : Phi="  <<doselRotation.getPhi()
                                <<",Theta="<<doselRotation.getTheta()
                                <<",Psi="  <<doselRotation.getPsi()<<G4endl
        <<space<<"       Absolute : Phi="<<mImage.GetTransformMatrix().getPhi()
                                <<",Theta="<<mImage.GetTransformMatrix().getTheta()
                                <<",Psi="<<mImage.GetTransformMatrix().getPsi()<<G4endl;
  space.resize(space.size()-3);*/
  //////////////////////////////////////////////////////////////////////////

  // If the mother's intersects the voxel.
  if(motherDoselOverlapCubicVolume==0.)
  {
    //G4cout<<space<<" *** IS NOT IN THE DOSEL N°"<<index<<" ***"<<G4endl;
    if(Generation>0)
      space.resize(space.size()-3);
    return std::make_pair(0.,0.);
  }

  // Calculation for daughter(s) ///////////////////////////////////////////
  if(motherLV->GetNoDaughters()>0) 
  {
    for(int i=0;i<motherLV->GetNoDaughters();i++)
    {
      G4VPhysicalVolume*  daughterPV(motherLV->GetDaughter(i));
      G4LogicalVolume*    daughterLV(daughterPV->GetLogicalVolume());
      G4VSolid*           daughterSV(daughterLV->GetSolid());

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
     
      std::pair<double,double> daughterIteration(VoxelIteration(daughterPV,Generation+1,daughterRotation,daughterTranslation,index));
      motherProgenyMass+=daughterIteration.first; 
      motherProgenyCubicVolume+=daughterIteration.second;
    }
  }
  //////////////////////////////////////////////////////////////////////////

  // Mother mass & volume //////////////////////////////////////////////////
  double motherCubicVolume(motherSV->GetCubicVolume());
  motherMass=motherCubicVolume*motherDensity;
  motherProgenyMass+=motherMass;
  motherProgenyCubicVolume+=motherCubicVolume;
  //////////////////////////////////////////////////////////////////////////

  // Mother information dump ///////////////////////////////////////////////
  /*G4cout<<space<<"  Mother original volume     : "<<G4BestUnit(motherLV->GetSolid()->GetCubicVolume(),"Volume")<<G4endl
        <<space<<"  Mother w/o daughters volume: "<<G4BestUnit(motherCubicVolume,"Volume")<<G4endl
        <<space<<"  Mother  +  daughters volume: "<<G4BestUnit(motherProgenyCubicVolume,"Volume")<<G4endl
        <<space<<"  Mother Mass                : "<<G4BestUnit(motherMass,"Mass")<<G4endl
        <<space<<"  Mother  +  daughters mass  : "<<G4BestUnit(motherProgenyMass,"Mass")<<" (MotherProgenyMass)"<<G4endl;*/
  //////////////////////////////////////////////////////////////////////////

  if(Generation>0)
    space.resize(space.size()-3);

  return std::make_pair(motherProgenyMass,motherProgenyCubicVolume);
}
//-----------------------------------------------------------------------------
