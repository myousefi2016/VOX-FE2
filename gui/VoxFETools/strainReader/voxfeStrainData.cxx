
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

#include "voxfeStrainData.h"
#include "../VoxFEDefines.h"


#define VOXFE_DISPLACEMENT_DATA_SIZE 3
#define VOXFE_ALL_STRAIN_DATA_SIZE   12

//----------------------------------------------------------------------------
// HELPER CLASS

voxfeStrainData::voxfeStrainData() {

  sdata.ren_node    = -1;
  sdata.glob_node   = -1;

  for( int k=0; k<3; k++ ) {
    sdata.location[k] = 0;
    sdata.disp[k]     = 0.0;
    sdata.nStrain[k]  = 0.0;
    sdata.sStrain[k]  = 0.0;
    sdata.pStrain[k]  = 0.0;
  }

  SetDisplacementOnly( true );

  nNodes    = 0;
  nodesRead = 0;
  dim[0] = dim[1] = dim[2] = 0;
  materialsFile = "0";
}

int voxfeStrainData::ReadHeader( ifstream& fin ) {

  nNodes = 0;
  string sdummy;
  int idummy;

  if( !this->GetDisplacementOnly() ) {

    //skip the first 5 line header
    char dummy[VOXFE_DUMMY_LINE_LENGTH_MAX];
    for( int k=0; k<5; k++ )
      fin.getline( dummy, VOXFE_DUMMY_LINE_LENGTH_MAX );

    fin >> sdummy >> idummy;                                      //Last_RI
  }

  fin >> sdummy >> nNodes;                                      //RN
  fin >> sdummy >> dim[0] >> dim[1] >> dim[2];                  //RX1,RY1,RZ1 (nodes in x,y,z)

  int border[3]; //fixme: remove?
  if( this->GetDisplacementOnly() )
    fin >> sdummy >> materialsFile;
  else
    fin >> sdummy >> border[0] >> border[1] >> border[2];         //Arbitrary border

  //test
  //cout << "nNodes: " << nNodes << " (" << dim[0] << "," << dim[1] << "," << dim[2] << ")\n";

  if( fin.good() ) { //hardline check for all failure bits
    nodesRead = 0;

    cout << "Read header.....\n";
    return 1;
  }

  cout << "Could not read header.....\n";
  return 0;
}

int voxfeStrainData::ReadLineAllStrain( ifstream& fin, unsigned long& _nodesRead ) {

  for( int k=0; k<3; k++ ) fin >> sdata.location[k];
  fin >> sdata.ren_node;
  fin >> sdata.glob_node;
  for( int k=0; k<3; k++ ) fin >> sdata.disp[k];
  for( int k=0; k<3; k++ ) fin >> sdata.nStrain[k];
  for( int k=0; k<3; k++ ) fin >> sdata.sStrain[k];
  for( int k=0; k<3; k++ ) fin >> sdata.pStrain[k];

  if( !fin.eof() ) nodesRead++;
  _nodesRead = nodesRead;

  if( fin.good() ) return 1;
  return 0;
}


int voxfeStrainData::ReadLineDispOnly( ifstream& fin, unsigned long& _nodesRead ) {

  fin >> sdata.ren_node;
  for( int k=0; k<3; k++ ) fin >> sdata.disp[k];

  if( !fin.eof() ) nodesRead++;
  _nodesRead = nodesRead;

  if( fin.good() ) {
    //cout << "Read line of data....\n";
    return 1;
  }

  //cout << "Failed to read data line.....\n";
  return 0;
}


int voxfeStrainData::ReadLineLocAndDisp( ifstream& fin, unsigned long& _nodesRead ) {

  fin >> sdata.ren_node;
  for( int k=0; k<3; k++ ) fin >> sdata.location[k];
  for( int k=0; k<3; k++ ) fin >> sdata.disp[k];

  if( !fin.eof() ) nodesRead++;
  _nodesRead = nodesRead;

  //debug
  /*cout << sdata.ren_node << "  ";
  for( int k=0; k<3; k++ ) cout << sdata.location[k] << "  ";
  for( int k=0; k<3; k++ ) cout << sdata.disp[k] << "  ";
*/
  if( fin.good() ) {
    //cout << "Read line of data....\n";
    return 1;
  }

  //cout << "Failed to read data line.....\n";
  return 0;
}

void voxfeStrainData::GetHeaderData( unsigned long& _nNodes, unsigned long& _nodesRead,
                        unsigned long _dim[3] ) {

  _nNodes    = nNodes;
  _nodesRead = nodesRead;
  _dim[0] = dim[0];  _dim[1] = dim[1];  _dim[2] = dim[2];
  //_border[0] = border[0];  _border[1] = border[1];  _border[2] = border[2];

}

void voxfeStrainData::SetDisplacementOnly( bool disp_only ) {
  sdata.dispOnly = disp_only;
  if( sdata.dispOnly ) nDoubles = VOXFE_DISPLACEMENT_DATA_SIZE; //6 with loc data;
  else nDoubles = VOXFE_ALL_STRAIN_DATA_SIZE;
}



double* voxfeStrainData::ReadBMUSolverDisplacementData( ifstream& fin, const unsigned long& actualNodeCount  ) {

  //set whether disp only data in reader
  this->SetDisplacementOnly( true );
  if( !this->ReadHeader( fin ) ) return 0;
  this->SetNumberOfNodes(actualNodeCount);  //hack to fix num nodes error in file

  //retrieve the header data
  unsigned long nodesRead, numNodes;
  this->GetHeaderData( numNodes, nodesRead, this->dim );  //numNodes should be correct now

  //make space for a single large table of data to be returned to the caller
  double* strainTable = new (std::nothrow) double [ this->GetNumDoubles() * actualNodeCount ];  //fixme: can we smartpointer this
  if( !strainTable ) return 0;

  while( this->ReadLineDispOnly( fin, nodesRead ) > 0 ) {

    //extract the bits we want to the vtk arrays
    unsigned long index = (nodesRead-1) * this->GetNumDoubles();
    const SData* imported_data = this->GetData();

    for( int p=0; p<3; p++ ) {
      strainTable[ index + p ] =   imported_data->disp[p];     //displacement at node
    }
  }

  return strainTable;
}

double* voxfeStrainData::ReadPETSCSolverDisplacementData( ifstream& fin, map< unsigned long, double* >& mapGlobalToStrainData ) {

  //set whether disp only data in reader
  this->SetDisplacementOnly( true );
  if( !this->ReadHeader( fin ) ) return 0;

  //retrieve the header data
  unsigned long nodesRead, numNodes;
  this->GetHeaderData( numNodes, nodesRead, this->dim );  //numNodes should be correct now

  //make space for a single large table of data to be returned to the caller
  double* strainTable = new (std::nothrow) double [ this->GetNumDoubles() * numNodes ];  //fixme: can we smartpointer this
  if( !strainTable ) return 0;

  mapGlobalToStrainData.clear();
  unsigned long dimXY = dim[0]*dim[1], glob_node;
  while( this->ReadLineLocAndDisp( fin, nodesRead ) > 0 ) {

    //extract the bits we want to the vtk arrays
    unsigned long index = (nodesRead-1) * this->GetNumDoubles();
    const voxfeStrainData::SData* imported_data = this->GetData();


    for( int p=0; p<3; p++ ) {
      strainTable[ index + p ] =   imported_data->disp[p];     //displacement at node
    }

    //work out global node num as index to map
    glob_node = (imported_data->location[2] * dimXY)  +
                (imported_data->location[1] * dim[0]) +
                imported_data->location[0];

    mapGlobalToStrainData.insert( pair< unsigned long, double* >( glob_node, &(strainTable[index]) ) );
  }

  return strainTable;

}

