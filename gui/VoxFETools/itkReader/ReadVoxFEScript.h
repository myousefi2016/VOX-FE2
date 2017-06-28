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

#ifndef READ_VOXFE_SCRIPT_H
#define READ_VOXFE_SCRIPT_H

#include "math.h"
#include "string.h"
#include "stdlib.h"

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <limits>
#include <iomanip>
#include <algorithm>

//#include <memory>
using std::vector;
using std::map;
using std::set;
using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::ifstream;
using std::ofstream;
using std::fstream;
using std::ostream;
using std::skipws;
using std::pair;
using std::ios;
using std::setprecision;


#include "../VoxFEDefines.h"


//--------------------------------------------------------------------------------------------------
//------------------------------------------ Structs -----------------------------------------------
//--------------------------------------------------------------------------------------------------

//FIXME: could convert to template types here?
typedef unsigned int rvsIntegerType;
typedef long int     rvsLongIntegerType;

/** Basic element data */
template<class T=rvsIntegerType>
struct Element {
#ifdef VOXFE_READ_DEBUG
  T  elNum;  //element number, dummied out for now
#endif
  unsigned short elType;  //element material type (only <= 256 expected)
  T  index[3]; //position
};

/** Material properties data */
struct MaterialType {
  double YoungsModulus;
  double PoissonsRatio;
  bool   RemodelFlag;
  double* LSM;

  MaterialType() {
    YoungsModulus =  0.0;
    PoissonsRatio = -1.0;
    RemodelFlag   = true;
    LSM           =  0;  //this is crucial as not always set, but will destroyed in Model
  }

  MaterialType(double ym, double pr, bool flag=true) {
    YoungsModulus = ym;
    PoissonsRatio = pr;
    RemodelFlag   = flag;
    LSM           =  0;
  }
};

/** Boundary condition data */
template<class T=rvsIntegerType>
struct Constraint {
  T nodeIndex[3];
  double force[3];
  unsigned short axisFixed[3];
  char type; //F=force, N=nodal, P=preserve(ie no remodelling)
};

/** Hold metadata for the model */
template<class T=rvsIntegerType>
struct ModelMetaData {

  ModelMetaData() {
    materialMap.clear();
  }

  //need destructor as we have a pointer member in materialMap
  ~ModelMetaData() {

     map< int, MaterialType >::iterator it = materialMap.begin();
     for( ; it != materialMap.end(); ++it )
       if( it->second.LSM ) delete [] it->second.LSM;
  }

  map< int, MaterialType > materialMap;  ///< Map from material id to material props

  double voxSize[4];   ///< We expect cubic voxels, but allow that each axis might be scaled. The 4th parameter is a universal scaling (not used).
  string algorithm;    ///< Name of solver algorithm type
  unsigned int maxIterations;   ///< Specify the max number of iterations for gradient solvers
  double tolerance;             ///< Specify the min tolerance for gradient solvers

  string scanType;  ///< Name of scan-type (not used)
  int scanDim[3];   ///< Specify the size of the 3D scan data block in voxels
  T X, XY, XYZ;     ///< Axis dimensions X, XY and XYZ are maintained for nodes
  string voxelFile; ///< The name of the model data file

  vector< Constraint<T> > constraint;  ///< Store constraint (boundary condition data)

  /** Get node indices from the unique node index */
  inline void GetIndexXYZ(const T& index, T xyz[3]) const {

    T rem;
    xyz[2] = index / XY; // "z"
    rem    = index % XY;
    xyz[1] = rem / X;     // "y"
    xyz[0] = rem % X;     // "x"
  }

  /** Read the material properties file */
  void ReadMaterialsFile( const char* MatFileName ) {

    fstream fin( MatFileName );
    if( !fin.is_open() )
      return;

    int m, flag;
    double ym, pr;
    while( !fin.eof() ) {

      fin >> m >> ym >> pr >> flag;
      materialMap.insert(  std::pair<int, MaterialType>( m, MaterialType(ym,pr,flag) ) );
    }

    fin.close();
  }

  /** Read constraint data from file
   *
   *  We allow 3 types:
   *  N - nodal constraint
   *  F - Force boundary condition
   *  P - Preserve the node at this location (used for controlling remodelling growth)
   *
   * */
  void AddConstraint( ifstream& fin, bool preserve=false ) {

    if( (!fin.is_open()) || (fin.eof()) )
      return;

    Constraint<T> ct;
    for( int p=0; p<3; p++ ) {
      fin >> ct.nodeIndex[p];
      ct.nodeIndex[p] -= 1;   //switch from 1-based to zero based indices
    }

    int nodal_con    = 0; //ct.axisFixed[0] + ct.axisFixed[1] + ct.axisFixed[2];
    double force_con = 0.0;
    for( int p=0; p<3; p++ ) {
      fin >> ct.force[p];
      force_con += fabs(ct.force[p]);
    }
    for( int p=0; p<3; p++ ) {
      fin >> ct.axisFixed[p];
      nodal_con += ct.axisFixed[p];
    }

    //decide type   FIXME: Is this sufficient to distinguish?????
    if( preserve )
      ct.type = 'P';
    else if( (nodal_con > ZERO_TOL) && (force_con < ZERO_TOL) )
      ct.type = 'N';
    else if( (nodal_con < ZERO_TOL) && (force_con > ZERO_TOL) )
      ct.type = 'F';
    else {
      cout << "===================================\n"
              "** Error reading constraint type **\n"
              "===================================\n";

      ct.type = '0';
   }

#ifdef VOXFE_CONSTRAINT_DEBUG
   //debug =================================================
   for( int p=0; p<3; p++ ) cout << ct.nodeIndex[p] + 1 << "  ";
   cout << "    ";

   for( int p=0; p<3; p++ ) cout << ct.force[p] << "  ";
   cout << "    ";

   for( int p=0; p<3; p++ ) cout << ct.axisFixed[p] << "  ";
   cout << "    " << ct.type << "\n";
   //debug =================================================
#endif

    constraint.push_back( ct );
  }

};


/** Used for reading image data via a '.voxfe' file.
 *
 *  Groups add the capability of defining material ranges (eg. bone 1-50) and
 *  border widths. The .voxfe file also fixes the voxel size.
 *
 *  .voxfe files/image data are converted to VoxFE script files and then re-read.
 */
struct ImageGroupData {

  /** Constants to define how materials are to be grouped */
  enum {

    LABEL_PRESENT                   = 1,
    LABEL_IN_GROUP                  = 2,
    LABEL_PRESENT_IN_GROUP_NO_REMOD = 3,
    LABEL_REMOD                     = 4,
    LABEL_PRESENT_IN_GROUP_REMOD    = 7,

    MAX_MATERIALS = 256
  };

  struct CMP_STR {
     bool operator()(char const *a, char const *b) const
     {
        return ::strcmp(a, b) < 0;
     }
  };

  /** Define a group of image labels as one material type */
  struct GroupDef {
    string MaterialName;  ///< Group material name eg 'bone'
    short Start;          ///< Where the group numbering starts (inclusive)
    short Stop;           ///< Where the group numbering stops (inclusive)

    short RemodelFlag;    ///< if 0, cannot remodel this group

    GroupDef(const string& name,
             const short& begin,
             const short& end,
             const short& flag):
      MaterialName(name), Start(begin), Stop(end), RemodelFlag(flag) {}
  };

  typedef map<const char*, MaterialType, CMP_STR> RefMaterialMap;
  typedef RefMaterialMap::iterator       RefMaterialMapIterator;
  typedef RefMaterialMap::const_iterator RefMaterialMapConstIterator;

  vector<GroupDef> Group;  ///< The list of available groups
  string ImageFile;        ///< The name of the imag file to be read
  double VoxelSize;        ///< Voxel size (cubic voxels assumed)
  bool TallyDone;          ///< Does the supplied range of materials match what was read
  int BorderOffset[3];     ///< Allow a an image border to be defined (for remodelling growth)

  /** Constructor */
  ImageGroupData();

  /** Destructor */
  ~ImageGroupData();

  /** Maintain a list of labels present (and add single groups if any) */
  void TallyGroups( const set<short>& labels_present );

  /** Write out a materials file, if none present */
  void WriteMaterialsFile( const char* MatFileName, const RefMaterialMap& RefMaterials );


protected:

  /** Keep a tally of which materials are present.
      0 indicates not present, >=1 otherwise
  */
  unsigned char Tally[MAX_MATERIALS];

  /** Find labels which don't appear to have been grouped -- maybe they indicate an error   */
  void GetNonGroups( set<short>& non_groups );
};

//--------------------------------------------------------------------------------------------------
//------------------------------------------ Funcs -------------------------------------------------
//--------------------------------------------------------------------------------------------------

/** Read Voxfe script file of model metadata.
 *
 */
template<class T>
bool rvsReadModelData( const char* file, const int _border[3], ModelMetaData<T>& model );

/** Read Voxfe model/voxel file. NB element is reallocated (but see detail below).
    \param element will be reallocated to store element index and type
    \param readElements determines whether elements are retained.
    \param nodes is used to store an in-memory mapping from unique node number to the
           sequence 1...n, for n nodes
    \param startValue sets the start value for nodes to be mapped to (0 or 1
           probably being preferable)
    \return true if no error was encountered
*/
template<class T>
bool rvsReadVoxelFile( const char* file,
                    Element<T>*& element,
                    T& n_elements,
		            const ModelMetaData<T>& model,
                    T minmaxZ[2],
                    map<T, T>& nodes,
                    const int _border[3],
                    T startValue=1 );


// --------------------------------------------------------------------------------------------------------
// //fixme:  Havent added in the border thing to routines below  ------------------------------------------
// --------------------------------------------------------------------------------------------------------


/** Read Voxfe model/voxel file to create node mapping, either to map ('nodes') or file.
    \param nodes is used to store a mapping from unique node number to the sequence startValue...n, for n nodes
           if nodeMapToFile != NULL
    \param startValue sets the start value for nodes to be mapped to (0 or 1 probably being preferable)
    \param nodeMapToFile If readElements is *false*, setting this parameter to non-NULL sends the node-mapping
           to binary file (rather than being stored in 'nodes' map). All potential global nodes will be
           mapped: the binary file will contain -1, where the node does not exist and will be >= 0
           where it does. The mapping is then random access by default since every possible node on the
           grid is represented and we can look up via: global node number = z*NX*NY + y*NX + x
           NB. If not NULL, the node count is return via this arg.
    \param nodeMapFile names the file to be used if nodeMapToFile > 0
    \return true if no error was encountered
*/
template<class T>
bool rvsCreateNodeMap( const char* file,
                    T& n_elements,
		            const ModelMetaData<T>& model,
                    T minmaxZ[2],
                    map<T, T>& nodes,
                    T startValue,
                    T* nodeMapToFile=NULL,
                    const char* nodeMapFile=NULL );


/** Helper func for rvsCreateNodeMap, read an element */
template<class T>
bool rvsCreateNodeMap_ReadElement( std::ifstream& fin,
                                   Element<T>& element,
                                   const ModelMetaData<T>& model,
                                   T minmaxZ[2],
                                   T gnode[NODESPERELEMENT] );

/** Helper func for rvsCreateNodeMap, process element with layer swapping if needed
    \param bin_file is the output binary map file

*/
template<class T>
void rvsCreateNodeMap_SwapLayers( std::fstream* bin_out,
                                  set<T>*& topSet,
                                  set<T>*& botSet,
                                  T& currentLayer,
                                  T& nodeCount );

/** Helper func for rvsCreateNodeMap, check existing map-file entry before writing
    \param bin_out is the output binary map file

*/
template<class T>
bool rvsCreateNodeMap_WriteMapEntry( std::fstream* bin_out,
                                     const T& location,
                                     const T& new_value );


/** Helper func for rvsCreateNodeMap, remap the node numbering into sequence
    \param bin_out is the output binary map file
    \param size gives the total number of entries (each of sizeof(T)) in the map file

*/
template<class T>
T rvsCreateNodeMap_RemapNodes( std::fstream* bin_out, const T& size );


/** Open the voxel file and read in chunk elements.
    \param element will be reallocated to store element index and type for chunk elements.
    Use alternate form (which takes ifstream as an arg) to read next chunk.
    If a file read failure is encountered, the stream will simply be returned.
*/
template<class T>
T rvsReadVoxelFileInChunks( ifstream& fin, Element<T>*& element,
                            T& n_elements, const unsigned int& chunk );


/** Read in chunk elements from the supplied voxel file stream.
    \param element is assumed to be allocated (and will not be reallocated).
    If a file read failure is encountered, the number of entries read will be returned.
*/
template<class T>
T rvsReadVoxelFileInChunks( ifstream& fin, Element<T>*& element, const unsigned int& chunk );



/** Gather node indices for the given voxel at index */
template<class T>
void rvsGetVoxelNodes( const ModelMetaData<T>& model, const T index[3],
		    T a[NODESPERELEMENT] );

/** After the unique node set has been collected, they will generally be
    an ordered but non-sequential and starting at some arbitrary number.
    Map these to the range 1... N by default.
*/
template<class T>
void rvsMapNodes( map<T, T>& nodes, T startValue );

/** Debug check on node numbering */
template<class T>
void rvsPrintNodes( const ModelMetaData<T>& model, const map<T, T>& nodes,
                 ofstream& fout );

/** Read the local stiffness matrices from a text file into the supplied model.
    Only entries in the material map will be read in.
*/
template<class T>
bool rvsReadLSMs( const char* lsmFile, const int& lsmSize, ModelMetaData<T>& model );

/** Find min and max indices.
    \param min_index should be initialised to the current known min value.
    \param max_index should be initialised to known max value.
    \param initialise if true, initialise min_index and max_index first.
*/
template<class T>
void rvsFindMinMax( const Element<T>* element, const T& n_elements, T min_index[3], T max_index[3], bool initialise=true );



/** Read a voxfe image-group file (which allows different materials
    to be considered as one group eg bone).
    \param filepath The complete file path to the .voxfe file.
    \return NULL if the file cannot be opened, or an allocated ImageGroupData object (user to delete) otherwise.
*/
ImageGroupData* rvsReadImageGroupFile( const char* filepath );

/** Convenience function to obtain the set of element material types present.
    \param element will be iterated through to collect the set
    \param label_set will contain the set of materials
*/
template<class T>
void rvsGetMaterialSet( const Element<T>* element, const T& n_elements, set<short>& label_set );

/** Write the global node numbers out to (binary) file in the order of the map. Note
    that the first element stored is the number of nodes.
    \param nodes is used to store an in-memory mapping from unique node number to the
           sequence 1...n, for n nodes
    \return true if no error was encountered
*/
template<class T>
bool rvsWriteNodeMapKeys( const char* bin_file,
                      const map<T, T>& nodes );

/**  Ascii version of above.
 *
 */
template<class T>
bool rvsWriteNodeMapKeysASCII( const char* asc_file,
                      const map<T, T>& nodes );

/** Read the global node numbers from a binary file.
    \param numNodes is read in from the file if successful
    \return array of size numNodes or NULL if error was encountered
*/
template<class T>
T* rvsReadNodeMapKeys( const char* bin_file, T& numNodes );



/** Write the global node numbers (type T) relating to element N0s out to (binary) file in numeric order, with the
 *  associated material (unsigned char). Note that the first element stored is the number of elements.
    \param element is the unordered array of elements.
    \return true if no error was encountered
*/
template<class T>
bool rvsWriteElementIDs( const char* bin_file, const Element<T>* element, const T& n_elements, const ModelMetaData<T>& model );

/**  Ascii version of above.
 *
 */
template<class T>
bool rvsWriteElementIDsASCII( const char* asc_file, const Element<T>* element, const T& n_elements, const ModelMetaData<T>& model );


 /** Read the global node number IDs (equivalent to N0) and materials of elements from a binary file.
    \param numNodes is read in from the file if successful
    \return array of size numNodes or NULL if error was encountered
*/
template<class T>
void rvsReadElementIDs( const char* bin_file, T& numNodes, T*& IDs, unsigned char*& materials );


//===================================================================================================
//================================ Write funcs ======================================================
//===================================================================================================

/**  Write out a basic VoxFE script file, nominating the model
 *   \param counter is added to the output files described in
 *          scriptfile if > 0.
 */
void rvsWriteVoxFEScriptStub( const char* scriptfile,
                              const char* modelfile,
                              double spacing[3],
                              int size[3],
                              const char* algorithm,
                              bool displacement_only=false,
                              const char* counter=0 );

/** Write out the voxel file with the given border.
 *
 */
template<class T>
bool rvsWriteVoxelFile( const char* file,
                        const Element<T>*& element,
                        const T& n_elements,
                        const int _border[3] );



#endif

