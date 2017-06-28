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

#include <vtkObjectFactory.h>
#include <vtkPolyDataReader.h>
#include <vtkSmartPointer.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkFieldData.h>
#include <vtkPointData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkCharArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkMath.h>
#include <vtkPolyData.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkIdList.h>

#include "../VoxFEDefines.h"
#include "voxfeStrainFilter.h"
#include "voxfeRecomputeStrainFilter.h"


//Try and indicate an error
#define CRAZY_STRAIN_VALUE -99.0


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
vtkStandardNewMacro(voxfeStrainFilter);

//----------------------------------------------------------------------------
voxfeStrainFilter::voxfeStrainFilter() : strainTable(0), StrainFile(0), nodeMap(0)
{
  this->SetNumberOfInputPorts(1);  //extra input for reader -- for the model surface point data
  this->SetNumberOfOutputPorts(1);
  
  strainData = new voxfeStrainData;
}

voxfeStrainFilter::~voxfeStrainFilter() {

  if( strainData ) delete strainData;
  if( strainTable ) delete[] strainTable;
  if( nodeMap ) delete[] nodeMap;
}

//----------------------------------------------------------------------------
void voxfeStrainFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
 
  os << indent << "File Name: "
     << (this->StrainFile ? this->StrainFile : "(none)") << "\n";
     
}

//----------------------------------------------------------------------------
int voxfeStrainFilter::FillInputPortInformation( int port, vtkInformation* info )
{

  switch (port)
  {
    //model data
    case 0:
      info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet" );
      return 1;

  default:
      vtkWarningMacro("Invalid input port " << port << " requested.");
      break;
  }
  
  return 0;
}


//----------------------------------------------------------------------------
void voxfeStrainFilter::SetInput0Connection(vtkAlgorithmOutput* algOutput)
{
  int port = 0;
  
  std::cout << "Setting connection on port: " << port << "\n";
  this->SetInputConnection(port, algOutput); //this replaces connections without adding
}
    

int voxfeStrainFilter::RequestData(
           vtkInformation *vtkNotUsed(request),
           vtkInformationVector **inputVector,
           vtkInformationVector *outputVector) {

  
  //[1]. get the first info object and model
  vtkInformation *modelInfo     = inputVector[0]->GetInformationObject(0);
  vtkMultiBlockDataSet* model   = vtkMultiBlockDataSet::SafeDownCast( modelInfo->Get(vtkDataObject::DATA_OBJECT()) );
  if( !model ) { 
    cerr << "Cannot retrieve the model data " << endl;
    return 0;
  }

  //set this for convenience; we must iterate over each block  
  int numBlocks                 = model->GetNumberOfBlocks();
  
  //we also need the voxel size which was stored as field data
  vtkFieldData* modelFieldData  = model->GetFieldData();    // .... and retrieve the attributes
  vtkDoubleArray *dvox = static_cast<vtkDoubleArray*>(modelFieldData->GetArray( "voxelsize" ));  //VOXFE_VOXELSIZE "voxelsize"
  double voxelsize = -1;
  if( dvox ) 
    voxelsize = *(dvox->GetPointer(0));
  else {
    cerr << "*** Voxel size cannot be set; operation cannot be completed ***" << endl << endl;
    return 0;
  }
  cout << "\nAdd strains acquired VoxFE voxelsize: "  << voxelsize << "\n\n";

  // =================================== For OLD SOLVER Displacement output ================================
  //retrieve the model file name, so we can get to the nodemap (if needed -- ie using old solver...)
  string nodeMapfile;
  if( this->nodeMap )
    delete[] this->nodeMap;
  this->nodeMap = 0;
  rvsIntegerType numMappedNodes = 0;

  if( this->GetBMUSolverOutput() ) {

    vtkCharArray* ca  = static_cast<vtkCharArray*>(modelFieldData->GetArray( VOXFE_VOXFEFILE )); // "voxfefile"
    if( ca ) {
      for(int p=0; p<ca->GetSize(); p++ ) {

        char ch = *(ca->GetPointer(p));
        if( ch != '\0' ) {
          nodeMapfile.push_back( ch );
        }
        else break;  //otherwise we get the extra junk (prob related to actual array size?)
      }
    }
    cout << "\n voxfeStrainFilter::RequestData obtained VoxFE file: \'"  << nodeMapfile << "\'\n\n";

    nodeMapfile.append( VOXFE_NODEMAP_FILE_EXTN );
    this->nodeMap = rvsReadNodeMapKeys( nodeMapfile.c_str(), numMappedNodes );
  }


  // get the output object
  vtkInformation *outInfo      = outputVector->GetInformationObject(0);
  vtkMultiBlockDataSet* output = vtkMultiBlockDataSet::SafeDownCast( outInfo->Get(vtkDataObject::DATA_OBJECT()) );
 
  //[2]. Read in strain/disp data to add point attributes
  ifstream strainfile( this->GetStrainFile() );
  if( !strainfile.is_open() ) {
  
    cerr << "Cannot open the strain file: " << this->GetStrainFile() << std::endl;
    output->ShallowCopy(model);
    return 0;
  }
  
  //test
  cout << "Opened strain file: " << this->GetStrainFile() << "\n";

  unsigned long numNodes, dim[3], dimXY, globalnode, NN = numMappedNodes;
  if( ReadStrainData( strainfile, dim, NN ) ) {

    //if this was used it can now be freed
    if( this->nodeMap ) {
      delete[] this->nodeMap;
      this->nodeMap = 0;
    }
  
    //fixme: dim here is assumed to be as voxfeVoxelOffset elsewhere and always 1,1,1
    dimXY = dim[0] * dim[1];

//#define VOXFE_DEBUG_STRAIN_CALCS
#ifdef VOXFE_DEBUG_STRAIN_CALCS
    ofstream debug_out("debug_strain_out.txt");
#endif
	  
    double pt[3];
    for( int blk=0; blk<numBlocks; ++blk ) {
       
      //multiblock -- need to retrieve/downcast block 0
      //We use PointSet here as the superclass of both vtkPolyData and vtkUnstructuredGrid,
      //but for this filter we're only really interested in the points anyway
      vtkPointSet* inputBlock = vtkPointSet::SafeDownCast( model->GetBlock(blk) );
      vtkPoints*     pointSet = inputBlock->GetPoints();
      vtkIdType numPointsInBlock = pointSet->GetNumberOfPoints();

      cout << "Block: " << blk << " (" << numPointsInBlock << " points)" << endl;       //debug

      //create vtk arrays
      vtkSmartPointer<vtkDoubleArray> displacement = vtkSmartPointer<vtkDoubleArray>::New();
        displacement->SetNumberOfComponents(3);
        displacement->SetNumberOfTuples(numPointsInBlock);  //this will be too many  -- hope to resize later
        displacement->SetName("Displacement");

      //declare smart pointers
      vtkSmartPointer<vtkDoubleArray> normalStrainX = 0;
      vtkSmartPointer<vtkDoubleArray> normalStrainY = 0;
      vtkSmartPointer<vtkDoubleArray> normalStrainZ = 0;
      vtkSmartPointer<vtkDoubleArray> shearStrainXY = 0;
      vtkSmartPointer<vtkDoubleArray> shearStrainYZ = 0;
      vtkSmartPointer<vtkDoubleArray> shearStrainZX = 0;
      vtkSmartPointer<vtkDoubleArray>  P1 = 0;
      vtkSmartPointer<vtkDoubleArray>  P2 = 0;
      vtkSmartPointer<vtkDoubleArray>  P3 = 0;
      vtkSmartPointer<vtkDoubleArray>  SED = 0;

      vtkIdType count = 0;

      //====================================================================================
      //================= read displacements, re-compute strains here ======================

      //... needed to compute the strains
      voxfeRecomputeStrainFilter::BMatrix B;
      voxfeRecomputeStrainFilter recomp;
      voxfeRecomputeStrainFilter::uVector u;
      voxfeRecomputeStrainFilter::eVector e;
      recomp.computeVoxelGradientMatrix( voxelsize, B );  //we now have the B matrix (constant for VoxFE)

#ifdef VOXFE_DEBUG_STRAIN_CALCS
      debug_out << "B:\n" << B << "\n\n";
#endif

      //to access material props
      double YM = -1., PR = 0., MOR, LAM;
      ModelMetaData<long> mmd;
      mmd.ReadMaterialsFile( this->strainData->GetMaterialsFile() );
      int materialsArrayIndex = -1;
      inputBlock->GetCellData()->GetArray( "MaterialType", materialsArrayIndex ); //test if MaterialType is present
      if( materialsArrayIndex < 0 )
        cerr << "-- Cannot find cell array: MaterialType --\n-- So cannot compute SED/stress --\n";
      else {

        //determine material props which should be same for every cell of same materialtype
        int material = static_cast<int>(*(inputBlock->GetCellData()->GetArray(materialsArrayIndex)->GetTuple(0)));
        map< int, MaterialType >::iterator it = mmd.materialMap.find( material );
        if( it !=  mmd.materialMap.end() ) {
          YM = it->second.YoungsModulus;
          PR = it->second.PoissonsRatio;
          MOR = recomp.getModulusOfRigidity( YM, PR );
          LAM = recomp.getLambda( YM, PR );
        }
        else {
          cerr << "\n/***********************************************************************************/\n";
          cerr <<   " Material: " << material << " could not be found; SED values will not likely be valid";
          cerr << "\n/***********************************************************************************/\n\n";

          //put something out to file too....
          ofstream log( "debug_strain_log.txt", std::ios::app );
          log << "** Material: " << material << " could not be found; SED values will not likely be valid **\n";
          log << "material: " << material << " has YM: "  << YM  << "  PR: "  << PR  << "\n";
          log << "MOR: " << MOR << "  LAM: " << LAM << "\n\n";
          log << "B:\n" << B << "\n\n";
          log.close();
        }

#ifdef VOXFE_DEBUG_STRAIN_CALCS
        debug_out << "material: " << material << " has YM: "  << YM  << "  PR: "  << PR  << "\n";
        debug_out << "MOR: " << MOR << "  LAM: " << LAM << "\n\n";
#endif

      }

      //add the displacements as point attribute data
      for(vtkIdType p=0; p<numPointsInBlock; p++ ) {

        pointSet->GetPoint(p, pt);
        GetGlobalNodeNumber( pt, voxelsize, dimXY, dim[0], globalnode );

        map<unsigned long, double*>::iterator it = mapGlobalToStrainData.find( globalnode );
        if( it != mapGlobalToStrainData.end() ) {

          displacement->SetTuple(p, (it->second) );

#ifdef VOXFE_DEBUG_STRAIN_CALCS
          debug_out << count << " : " << it->first << " : "
               << it->second[0] << ","
               << it->second[1] << ","
               << it->second[2] << "  ";
          debug_out << "  (" << pt[0] << "," << pt[1] << "," << pt[2] << ")\n";
#endif

          count++;

        }
        else {
          //indicate an error with a 'crazy' value
          double crazy[3] = { CRAZY_STRAIN_VALUE, CRAZY_STRAIN_VALUE, CRAZY_STRAIN_VALUE };
          displacement->SetTuple(p, crazy);

#ifdef VOXFE_DEBUG_STRAIN_CALCS
          debug_out << "Globalnode/Displ at point: " << p << " not found!\n";
#endif
        }
      }

      inputBlock->GetPointData()->AddArray( displacement );          //Add the array for use with paraview

      //=======================================================================================================
      //debug: issue error
      this->OutputIfCountError( pointSet->GetNumberOfPoints(), count, blk );

      //=======================================================================================================
      //compute strains per cell -- only if internal geometry exists
      vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::SafeDownCast( model->GetBlock(blk) );
      if( ugrid ) {

        vtkIdType numCellsInBlock = ugrid->GetNumberOfCells();
        normalStrainX = vtkSmartPointer<vtkDoubleArray>::New();
        normalStrainY = vtkSmartPointer<vtkDoubleArray>::New();
        normalStrainZ = vtkSmartPointer<vtkDoubleArray>::New();
        normalStrainX->SetNumberOfTuples(numCellsInBlock);
        normalStrainY->SetNumberOfTuples(numCellsInBlock);
        normalStrainZ->SetNumberOfTuples(numCellsInBlock);
        normalStrainX->SetName("NormalStrain_X");
        normalStrainY->SetName("NormalStrain_Y");
        normalStrainZ->SetName("NormalStrain_Z");

        shearStrainXY = vtkSmartPointer<vtkDoubleArray>::New();
        shearStrainYZ = vtkSmartPointer<vtkDoubleArray>::New();
        shearStrainZX = vtkSmartPointer<vtkDoubleArray>::New();
        shearStrainXY->SetNumberOfTuples(numCellsInBlock);
        shearStrainYZ->SetNumberOfTuples(numCellsInBlock);
        shearStrainZX->SetNumberOfTuples(numCellsInBlock);
        shearStrainXY->SetName("ShearStrain_XY");
        shearStrainYZ->SetName("ShearStrain_YZ");
        shearStrainZX->SetName("ShearStrain_ZX");

        P1 = vtkSmartPointer<vtkDoubleArray>::New();
        P1->SetNumberOfTuples(numCellsInBlock); P1->SetName("P1");
        P2 = vtkSmartPointer<vtkDoubleArray>::New();
        P2->SetNumberOfTuples(numCellsInBlock); P2->SetName("P2");
        P3 = vtkSmartPointer<vtkDoubleArray>::New();
        P3->SetNumberOfTuples(numCellsInBlock); P3->SetName("P3");
        SED = vtkSmartPointer<vtkDoubleArray>::New();
        SED->SetNumberOfTuples(numCellsInBlock); SED->SetName("SED");

        for(vtkIdType k=0; k<numCellsInBlock; k++ ) {

          vtkCell* cell = ugrid->GetCell( k );
          vtkIdList* pointIds = cell->GetPointIds();

          //set the displacements around each cell
          //cout << "Cell: " << k << " has point Ids: ";
          for(int m=0; m<pointIds->GetNumberOfIds(); m++ ) {
            //cout << pointIds->GetId(m) << " ";
            displacement->GetTuple( pointIds->GetId(m), pt );  //node order for VTK matches that for VoxFE fortunately
            int m3 = 3*m;
            u(m3  ) = pt[0];
            u(m3+1) = pt[1];
            u(m3+2) = pt[2];
          }
          //cout << "\n";

          //compute basic strains
          recomp.getStrainTensor( B, u, e );
          normalStrainX->SetTuple( k, &(e(0)) );
          normalStrainY->SetTuple( k, &(e(1)) );
          normalStrainZ->SetTuple( k, &(e(2)) );
          shearStrainXY->SetTuple( k, &(e(3)) );
          shearStrainYZ->SetTuple( k, &(e(4)) );
          shearStrainZX->SetTuple( k, &(e(5)) );

          //compute principal strains
          recomp.computePrincipalStrain( e.data(), pt );
          P1->SetTuple(k, pt    );
          P2->SetTuple(k, pt + 1);
          P3->SetTuple(k, pt + 2);

          //compute strain energy density -- this needs Y modulus and Poissons R
          //double strain_tensor[6];
          if( materialsArrayIndex >= 0 ) { //fixme: We haven't found the material type in the data.....

            //e.data() ... may not always get what we want
            //for(int q=0; q<6; q++) strain_tensor[q] = e(q);

#if 0
            recomp.computeStrainEnergyDensity( e.data(), YM, PR, pt[0] );
#else
            recomp.computeStrainEnergyDensity2( e.data(), MOR, LAM, pt[0] );
#endif
            SED->SetTuple(k, pt);
          }
          else {
            pt[0] = CRAZY_STRAIN_VALUE;
            SED->SetTuple(k, pt);
          }

        }

        inputBlock->GetCellData()->AddArray( normalStrainX );
        inputBlock->GetCellData()->AddArray( normalStrainY );
        inputBlock->GetCellData()->AddArray( normalStrainZ );
        inputBlock->GetCellData()->AddArray( shearStrainXY );
        inputBlock->GetCellData()->AddArray( shearStrainYZ );
        inputBlock->GetCellData()->AddArray( shearStrainZX );
        inputBlock->GetCellData()->AddArray( P1 );
        inputBlock->GetCellData()->AddArray( P2 );
        inputBlock->GetCellData()->AddArray( P3 );
        inputBlock->GetCellData()->AddArray( SED );

      } //if (ugrid) add do strain calcs

    }  //for each block

    if(this->strainTable) delete[] this->strainTable;
    this->strainTable = 0;
    mapGlobalToStrainData.clear(); //values will be void now anyway....
    if(this->nodeMap) delete[] this->nodeMap;

#ifdef VOXFE_DEBUG_STRAIN_CALCS
    debug_out.close();
#endif


  } //if( ReadStrainData....)

  output->ShallowCopy( model );
  return 1;
}
  

/*
 *  Fixme:: This routine is getting well contorted, but is trying to do several jobs:
 *
 *  1) we need to read in all the displacements (at least) and store them in a table
 *  so that we can work out later where they go -- via the global map. The problem being
 *  that the data may be in blocks, so we have lost the original order of the points.
 *
 *  2) we're still trying to support the formats of the old solver (all strain data)
 *  and the new petsc-based one (displ only) or some mix of the two. Hopefully,
 *  this can all be cleaned up later...
 *
 *  3) we're also still trying to support the formats of the old solver for just
 *  displacements, useful when the solver has to be run in parallel....
 */
double* voxfeStrainFilter::ReadStrainData( ifstream& fin, unsigned long dim[3], const unsigned long& actualNodeCount ) {

  //set whether disp only data in reader
  this->strainData->SetDisplacementOnly( true );

  if( !strainData->ReadHeader( fin ) ) return 0;
  
  //Kludge for BMU solver output
  if( this->GetBMUSolverOutput() ) {
    if( !this->nodeMap ) return 0;

    //use actual node count instead of proc 0 number
    strainData->SetNumberOfNodes(actualNodeCount);
  }

  //retrieve the header data
  unsigned long nodesRead, nNodes;
  strainData->GetHeaderData( nNodes, nodesRead, dim );
  
  //make space to store the strain data until we can pass onto vtk
  mapGlobalToStrainData.clear();

  if( mapGlobalToStrainData.max_size() < nNodes ) {
    cerr << "\n\n** Map of strain data may exceed system memory **\n\n";
  }

  if( this->strainTable ) delete[] this->strainTable;
  this->strainTable = 0;

  //make space for a single large table of data to be returned to the caller
  this->strainTable = new (std::nothrow) double [ this->strainData->GetNumDoubles() * nNodes ];  //fixme: can we smartpointer this

  if( !this->strainTable ) return 0;

  unsigned long glob_node, dimXY = dim[0]*dim[1];

  if( this->GetBMUSolverOutput() ) {

    unsigned long k = 0;
    while( strainData->ReadLineDispOnly( fin, nodesRead ) > 0 ) {

      //extract the bits we want to the vtk arrays
      unsigned long index = (nodesRead-1) * strainData->GetNumDoubles();
      const voxfeStrainData::SData* imported_data = strainData->GetData();

#ifndef NDEBUG
      if( !imported_data ) {
        for( int p=0; p<3; p++ ) {
          strainTable[ index + p ] =   CRAZY_STRAIN_VALUE;
        }
      }
      else {
        for( int p=0; p<3; p++ ) {
          strainTable[ index + p ] =   imported_data->disp[p];     //displacement at node
        }
      }
#else
      for( int p=0; p<3; p++ ) {
        strainTable[ index + p ] =   imported_data->disp[p];     //displacement at node
      }
#endif

      //use the global node value in nodeMap
      glob_node = this->nodeMap[k];
      k++;

      mapGlobalToStrainData.insert( pair< unsigned long, double* >( glob_node, &(strainTable[index]) ) );
    }
  }
  else {

    while( strainData->ReadLineLocAndDisp( fin, nodesRead ) > 0 ) {

      //extract the bits we want to the vtk arrays
      unsigned long index = (nodesRead-1) * strainData->GetNumDoubles();
      const voxfeStrainData::SData* imported_data = strainData->GetData();

#ifndef NDEBUG
      if( !imported_data ) {
        for( int p=0; p<3; p++ ) {
          strainTable[ index + p ] =   CRAZY_STRAIN_VALUE;
        }
      }
      else {
        for( int p=0; p<3; p++ ) {
          strainTable[ index + p ] =   imported_data->disp[p];     //displacement at node
        }
      }
#else
      for( int p=0; p<3; p++ ) {
        strainTable[ index + p ] =   imported_data->disp[p];     //displacement at node
      }
#endif

      //work out global node num as index to map
      glob_node = (imported_data->location[2] * dimXY)  +
                  (imported_data->location[1] * dim[0]) +
                  imported_data->location[0];

      mapGlobalToStrainData.insert( pair< unsigned long, double* >( glob_node, &(strainTable[index]) ) );
    }
  }
  
  return strainTable;
} 


void voxfeStrainFilter::OutputIfCountError( const vtkIdType& expected, const vtkIdType& count, const int& blk ) {

  if( count != expected ) {
    cerr << "\n************************************\n"
         << "Block: " << blk
         << " has bad count: " << count << "["
         << expected << "]\n"
         << "\n************************************\n" << endl;
  }
}

//we assume here that we should do +1 for the border and -1 to get Node 0 (ie we can ignore offset)
void voxfeStrainFilter::GetGlobalNodeNumber( const double* pt, const double& voxelsize,
                                     const double& dimXY, const double& dimX,
                                     unsigned long& globalnode  ) {

  globalnode = vtkMath::Round( pt[2]/voxelsize ) * dimXY  +
               vtkMath::Round( pt[1]/voxelsize ) * dimX +
               vtkMath::Round( pt[0]/voxelsize );

}


