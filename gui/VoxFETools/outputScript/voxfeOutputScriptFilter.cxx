/*
 * Vox-FE (Version 2): Voxel-based Finite Element Modelling
 *
 * Interface developed as a ParaView plugin by
 * Richard Holbrey and Michael J. Fagan
 *
 * Dept of Mechanical Engineering
 * University of Hull
 * Cottingham Road
 * HU6 7RX
 *
 * Vox-FE (Version 1) was created by Andreas Bitternas, George Sisias, Jia Liu and others.
 */


#include "vtkObjectFactory.h" //for new() macro
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPolyData.h"
#include "vtkPolyDataWriter.h"
#include "vtkCharArray.h"
#include "vtkIntArray.h"
#include "vtkDoubleArray.h"
#include "vtkFieldData.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkMultiBlockDataSet.h"

#include <iostream>
#include <string>
#include <fstream>
#include <limits>
using std::cout;
using std::cerr;
using std::string;
using std::ofstream;
using std::endl;
using std::fixed;
using std::numeric_limits;

#include "voxfeOutputScriptFilter.h"
#include "../itkReader/ReadVoxFEScript.h"

vtkStandardNewMacro(voxfeOutputScriptFilter);

//-----------------------------------------------------------------------------
voxfeOutputScriptFilter::voxfeOutputScriptFilter()
{
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}


//-----------------------------------------------------------------------------
voxfeOutputScriptFilter::~voxfeOutputScriptFilter()
{

}


//----------------------------------------------------------------------------
int voxfeOutputScriptFilter::FillInputPortInformation( int port, vtkInformation* info )
{

  switch (port)
  {
    //model data
    case 0:
      info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet" );
      return 1;

    //boundary condtions
    case 1:
      info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet" );      
      info->Set( vtkAlgorithm::INPUT_IS_REPEATABLE(), 1 );
      return 1;

  default:
      vtkWarningMacro("Invalid input port " << port << " requested.");
      break;
  }
  
  return 0;
}


//----------------------------------------------------------------------------
int voxfeOutputScriptFilter::RequestData(vtkInformation *vtkNotUsed(request),
                              vtkInformationVector **inputVector,
                              vtkInformationVector *outputVector)
{
  std::cout << "OutputScriptFilter::RequestData\n"
            << "Number of information objects on port 0: "
            << inputVector[0]->GetNumberOfInformationObjects() << "\n";

  // ============================== model data on port 0 ==============================================
  //[1]. get the first info object and model
  vtkInformation *inputModelInfo   = inputVector[0]->GetInformationObject(0);
  vtkMultiBlockDataSet* inputModel = vtkMultiBlockDataSet::SafeDownCast( inputModelInfo->Get(vtkDataObject::DATA_OBJECT()) );

  // .... and retrieve the voxfe file-name attribute if we can
  vtkFieldData* fd0 = inputModel->GetFieldData();
  
  string voxfefile;
  vtkCharArray* ca  = static_cast<vtkCharArray*>(fd0->GetArray( VOXFE_VOXFEFILE )); // "voxfefile"
  if( ca ) {
    for(int p=0; p<ca->GetSize(); p++ ) {
    
      char ch = *(ca->GetPointer(p));
      if( ch != '\0' ) {
        voxfefile.push_back( ch );
      }
      else break;  //otherwise we get the extra junk (prob related to actual array size?)
    }
  }
  cout << "\n OutputScriptFilter::RequestData obtained VoxFE file: \'"  << voxfefile << "\'\n\n";
  
  vtkDoubleArray *dvox = static_cast<vtkDoubleArray*>(fd0->GetArray( VOXFE_VOXELSIZE )); // "voxelsize"
  double voxelsize = -1;
  if( dvox ) 
    voxelsize = *(dvox->GetPointer(0));
  else {
    cerr << "*** Voxel size cannot be set; operation cannot be completed ***" << endl << endl;
    return 0;
  }
  cout << "\n VoxFE voxelsize: "  << voxelsize << "\n\n";
  
  // ============================== BC data on port 1 ==============================================
  //[2]. retrieve the BC info if we can but since we may have a number of
  // conditions (multiple_input=1), we have to run through them
  int numInfo1Obj = inputVector[1]->GetNumberOfInformationObjects();
  cout << "\nNumber of information objects on port 1: " << numInfo1Obj << "\n";

  //create the output file for constraints
  ofstream fout(VOXFE_DEFAULT_CONSTRAINT_FILE);  //assume we're in the correct dir

  //from http://stackoverflow.com/questions/554063/how-do-i-print-a-double-value-with-full-precision-using-cout
  int max_dbl_precision = numeric_limits<double>::digits10 + 1;
  fout.precision( max_dbl_precision ); //set a max precision for display of double

  for( int k=0; k<numInfo1Obj; k++ ) {  //for each defined BC

    vtkInformation* inputBCInfo = inputVector[1]->GetInformationObject(k);
    vtkUnstructuredGrid* ug     = GetBoundaryConditionFromUserSelection( inputBCInfo );

    if( !ug ) continue;  //fixme: why did this fail

    //handle the BC
    vtkPoints* pointSet = ug->GetPoints();
    vtkFieldData* fd1   = ug->GetFieldData();
    vtkIntArray*  ia    = static_cast<vtkIntArray*>(fd1->GetArray( VOXFE_CONSTRAINT_ATTRIBUTE ));
    vtkDoubleArray* da  = static_cast<vtkDoubleArray*>(fd1->GetArray( VOXFE_VECTOR_ATTRIBUTE ));
    vtkIdType n_points  = pointSet->GetNumberOfPoints();

    if( ! (ia && da && n_points) ) continue;

    int conditionType = -1, forceDistrib = -1;
    double pt[3], vec[3], pvec[3];
    double ep_average_mag = 0.0, ep_scaling = -1.0, user_mag = 0.0;
    
    conditionType = *(ia->GetPointer(0));
    forceDistrib  = *(ia->GetPointer(1));
    vec[0]   = *(da->GetPointer(0)); vec[1] = *(da->GetPointer(1)); vec[2] = *(da->GetPointer(2));
    user_mag = *(da->GetPointer(3));

    //[3]. Need to loop over Force-to-point conditions first to get the scaling
    //for the force magnitude
    if( conditionType == VOXFE_BC_FORCE_TO_POINT ) {
      for(vtkIdType p=0; p<n_points; p++ ) {

        pointSet->GetPoint(p, pt);

        //vec is the endpoint of the force vector, so subtract selected point
        for(int q=0; q<3; q++ ) {
          pvec[q] = vec[q] - pt[q];
          //cout << "pvec[" << q << "]: " << vec[q] << " - " << pt[q] << endl;
        }
        ep_average_mag += sqrt( pvec[0]*pvec[0] + pvec[1]*pvec[1] + pvec[2]*pvec[2] );  //we apply this below....
        //debug
        //cout << "    Endpoint condition sum (" << p << "): " << ep_average_mag << endl;
      }

      ep_average_mag /= n_points;  //change into an average magnitude (used for end-pt condition)

      //work out scaling factor for endpoint condition
      ep_scaling = user_mag/ep_average_mag;

      cout << "    Endpoint condition total vector magnitude: " << ep_average_mag
           << " across: " << n_points << " points" << endl << endl;
    }

    //[4]. Second pass to actually output data to file
    double scaling = -1;
    if( conditionType == VOXFE_BC_FORCE_PARALLEL ) {
      //in this case all the forces will be the same so compute here
      scaling = user_mag / sqrt( vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] );  //scaling for vector
      if( forceDistrib > 0 ) scaling /= double(n_points);
      pvec[0] = vec[0]*scaling;
      pvec[1] = vec[1]*scaling;
      pvec[2] = vec[2]*scaling;
    }

    for(vtkIdType p=0; p<n_points; p++ ) {

      //we want integer indices here -- might have to check for rounding errors
      pointSet->GetPoint(p, pt);

      if( conditionType == VOXFE_BC_NO_REMODEL )
        fout << "  PRESERVE_NODE_3D  " << "  ";
      else
        fout << "  SELECT_NODE_3D  " << "  ";

      fout
           << (int)( this->round(pt[0]/voxelsize) + voxfeVoxelOffset ) << " "
           << (int)( this->round(pt[1]/voxelsize) + voxfeVoxelOffset ) << " "
           << (int)( this->round(pt[2]/voxelsize) + voxfeVoxelOffset ) << " ";  //NB 1-based indices

      switch( conditionType ) {
      case VOXFE_BC_NODAL:
        fout << "0 0 0 "
             << vec[0] << " " << vec[1] << " " << vec[2] << "\n";
        break;
      case VOXFE_BC_FORCE_PARALLEL:
        fout << pvec[0] << " " << pvec[1] << " " << pvec[2];
        fout << " 0 0 0" << "\n";
        break;
      case VOXFE_BC_FORCE_TO_POINT:
        //in this case all the forces will scaled using the ep_average_mag computed earlier
        for(int q=0; q<3; q++ ) pvec[q] = vec[q] - pt[q];

        scaling = ep_scaling;  //scaling: user-supplied-magnitude / mean-point-force-magnitude
        if( forceDistrib > 0 ) scaling /= double(n_points);
        pvec[0] *= scaling;
        pvec[1] *= scaling;
        pvec[2] *= scaling;
        fout <<  pvec[0] << " " << pvec[1] << " " << pvec[2] << " 0 0 0\n";
        break;
      case VOXFE_BC_NO_REMODEL:
        fout << "0 0 0  0 0 0\n";
      }
    }
  }  //for each input BC
  
  fout.close();

  // ============================== Output ========================================================
  //[4]. Output
  vtkInformation *outputModelInfo = outputVector->GetInformationObject(0);

  //we use vtkMultiBlockDataSet here to pass the geometry through
  vtkMultiBlockDataSet  *outputModel = vtkMultiBlockDataSet::SafeDownCast( outputModelInfo->Get(vtkDataObject::DATA_OBJECT()) );
  outputModel->ShallowCopy(inputModel);

  return 1;
}


vtkUnstructuredGrid* voxfeOutputScriptFilter::
GetBoundaryConditionFromUserSelection( vtkInformation* inputBC ) {

  if( ! inputBC ) {
    cerr << "GetBC: Input null, returning... " << endl;
    return NULL;
  }

  //Ctrl-clicking in ParaView assembles these into blocks (so now they're blocks of blocks...)
  vtkMultiBlockDataSet* inputBC_mb  = vtkMultiBlockDataSet::SafeDownCast( inputBC->Get(vtkDataObject::DATA_OBJECT()));
  if( !inputBC_mb ) {
    cerr << "GetBC: Input cannot be cast to MB, returning... " << endl;
    return NULL;
  }

  //cout << "vtkMultiBlockDataSet (MB) cast of this info object has " << inputBC_mb->GetNumberOfBlocks() << " blocks" << endl;
  vtkMultiBlockDataSet* inputBC_mb0   = vtkMultiBlockDataSet::SafeDownCast( inputBC_mb->GetBlock( 0 ) );
  vtkUnstructuredGrid*  ug            = NULL;

  if( ! inputBC_mb0 ) {

#if 0
    cerr << "GetBC: Input 0 cannot be re-cast to MB, returning... " << endl;
    return NULL;
#else
    //try this for the heck
    ug =  vtkUnstructuredGrid::SafeDownCast( inputBC_mb->GetBlock( 0 ) );
    if( ug ) {
      cout << " .... cast straight to UG did obtain an unstruc grid!\n";
    }
#endif
  }
  else {
    //cout << "[p] Casting the " << p << "'th block of this MB retrieves another MB\n";
    //    cout << "with : " << inputBC_mb0->GetNumberOfBlocks() << " blocks" << endl;
    ug    =  vtkUnstructuredGrid::SafeDownCast( inputBC_mb0->GetBlock( 0 ) );
  }

  if( ug )
    cout << "Got BC object with: " << ug->GetNumberOfPoints() << endl << endl;

  return ug;
}


////////// External Operators /////////////

void voxfeOutputScriptFilter::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}


void voxfeOutputScriptFilter::SetInput0Connection(vtkAlgorithmOutput* algOutput)
{
  std::cout << "Connecting port: 0" /*<< port <<*/ " on Input0\n";
  this->SetInputConnection(0, algOutput);
}

void voxfeOutputScriptFilter::SetInput1Connection(int port, vtkAlgorithmOutput* algOutput)
{
  std::cout << "Adding/Connecting port: " << port << " on Input1\n";

  //this->SetInputConnection(port, algOutput); this replaces connections without adding, so need....
  this->AddInputConnection(port, algOutput);
}


int voxfeOutputScriptFilter::round(double number)
{
  //return (number >= 0) ? (int)(number + 0.5) : (int)(number - 0.5);
  return (int)(number > 0 ? number + 0.5 : ceil(number - 0.5));
}




