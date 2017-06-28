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


//==================================================================================================

#include "convertVoxFEScript.hxx"

//===================================== MAIN =======================================================

//#define VOXFE_CREATE_NODEMAP
#define VOXFE_CREATE_ELEMENTMAP

//typedef rvsIntegerType     ConvertScriptIntType;
typedef rvsLongIntegerType ConvertScriptIntType;


int main( int argc, char* argv[] ) {

  if( argc != 3 ) {
    //Assume the script file points to the correct model file
    cout << "Usage: " << argv[0] << "  Solver-script-file  TYPE" << endl;
    cout << "where TYPE is either A (or ANSYS), M (METIS),";
    cout << " F (FELT) or V (VTK)" << endl;
    exit(-1);
  }

  int border[3];
  border[0] = border[1] = border[2] = 0;

  //1. set the conversion type
  string conversion( argv[2] );
  ConversionType conversionType = NO_CONVERSION;
  if( !conversion.compare("a") || !conversion.compare("A") || !conversion.compare("ANSYS") )
    conversionType = ANSYS;
  if( !conversion.compare("m") || !conversion.compare("M") || !conversion.compare("METIS") )
    conversionType = METIS;
  if( !conversion.compare("f") || !conversion.compare("F") || !conversion.compare("FELT") )
    conversionType = FELT;
  if( !conversion.compare("v") || !conversion.compare("V") || !conversion.compare("VTK") )
    conversionType = VTK;

  if( conversionType == NO_CONVERSION ) {
    cerr << "Could not interpret the input parameters.... please re-enter" << endl;
    exit(-1);
  }

  //2. get the model meta data
  ModelMetaData<ConvertScriptIntType> model;
  if( !rvsReadModelData( argv[1], border, model ) ) {
    cerr << "Could not read the input script file...." << endl;
    exit(-2);
  }

  //3. create and fill voxel table
  Element<ConvertScriptIntType>* element = 0;
  ConvertScriptIntType  n_elements = 0;
  ConvertScriptIntType  minmaxZ[2];
  /*bool retainElementData = true;*/
  map<ConvertScriptIntType, ConvertScriptIntType> nodes;   //mapping from unique node index to position in map (latter used in output file)
  if( !rvsReadVoxelFile( model.voxelFile.c_str(), element, n_elements, model,  minmaxZ, nodes, border ) ) {
    cerr << "Could not read voxel file...." << endl;
    if( element ) delete[] element;
    exit(-3);
  }

  //4. Do conversion
  if( conversionType == ANSYS ) {   //Write out the ansys file

    string ansysfile = string(argv[1]) + string(".ansys");
    if( !WriteAnsysSolid45File( ansysfile.c_str(), nodes, element, n_elements, model ) ) {

      cout << "Could not write out the ansys file...." << endl;
    }
    else cout << "Written out the ansys file...." << endl;

  }
  else if( conversionType == METIS ) {  //Write out a Metis mesh file format

    string ansysfile = string(argv[1]) + string(".metis");
    if( !WriteMetisMeshFile( ansysfile.c_str(), nodes, element, n_elements, model ) ) {

      cout << "Could not write out the Metis file...." << endl;
    }
    else cout << "Written out the Metis file...." << endl;
    
#ifdef VOXFE_CREATE_NODEMAP
    //create a list of global IDs for nodes present
    string nm_file = string(argv[1]) + string(".nodemap");
    rvsWriteNodeMapKeys( nm_file.c_str(), nodes ); 
    
    //debug
    nm_file.append(".txt");
    rvsWriteNodeMapKeysASCII( nm_file.c_str(), nodes ); 
#endif
#ifdef VOXFE_CREATE_ELEMENTMAP
    //do same for elements
    string el_file = string(argv[1]) + string(".elmap");
    rvsWriteElementIDs( el_file.c_str(), element, n_elements, model );
    
    //debug
    el_file.append(".txt");
    rvsWriteElementIDsASCII( el_file.c_str(), element, n_elements, model );
#endif    
  }
  else if( conversionType == FELT ) {  //Write out a FELT mesh file

    string ansysfile = string(argv[1]) + string(".flt");
    if( !WriteFeltMeshFile( ansysfile.c_str(), nodes, element, n_elements, model ) ) {

      cout << "Could not write out the felt file...." << endl;
    }
    else cout << "Written out the felt file...." << endl;
  }
  else if( conversionType == VTK ) {  //Write out a VTK voxel mesh file

    string ansysfile = string(argv[1]) + string(".vtk");
    if( !WriteVTKMeshFile( ansysfile.c_str(), nodes, element, n_elements, model ) ) {

      cout << "Could not write out the vtk file...." << endl;
    }
    else cout << "Written out the vtk file...." << endl;
  }



#ifdef _DEBUG
  string nodefile = string(argv[1]) + string(".nodes");
  ofstream fnodes( nodefile.c_str() );
  PrintNodes( model, nodes, fnodes );
  fnodes.close();
#endif

  //======================== wrap up ===========================
  if( element ) delete[] element;

  return 0;
}



