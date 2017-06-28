
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

#ifndef __voxfeStrainData_h
#define __voxfeStrainData_h


#include "../itkReader/ReadVoxFEScript.h"

#define VOXFE_DUMMY_LINE_LENGTH_MAX 256

#include <iostream>
#include <fstream>
#include <map>
#include <string>
using std::ostream;
using std::istream;
using std::cout;
using std::endl;
using std::ifstream;
using std::map;
using std::string;
using std::pair;


/** \class voxfeStrainData
    \brief This class performs the reading of displacement/strain data table for voxfeStrainFilter.
*/
class voxfeStrainData {

  public:

    /** Intended to read a row of displacement and/or strain data (to construct a table) */
    typedef struct _SData {
      long location[3]; //node loc
      unsigned long ren_node;    //renumbered node #
      unsigned long glob_node;   //global node #
      double disp[3];   //displacement at node
      double nStrain[3]; //normal strain
      double sStrain[3]; //shear strain
      double pStrain[3]; //principal strain [1-3]

      //prob not best place for this, but while testing......
      bool dispOnly;  //if reading displacements only

    } SData;

    /** Constructor, set deafults */
    voxfeStrainData();

    /** Allow access to read data */
    const SData* GetData() { return &sdata; }
    
    /** Read header data from the supplied file */
    int ReadHeader( ifstream& fin );

    /** Read a row of pre-computed strain data from file */
    int ReadLineAllStrain( ifstream& fin, unsigned long& _nodesRead );

    /** Read a row of displacement data (and ID -- not used) */
    int ReadLineDispOnly ( ifstream& fin, unsigned long& _nodesRead );

    /** Read a row of displacement and 3D location data (and ID -- not used) */
    int ReadLineLocAndDisp( ifstream& fin, unsigned long& _nodesRead );
    
    /** Access the header data after reading */
    void GetHeaderData( unsigned long& _nNodes, unsigned long& _nodesRead,
                        unsigned long _dim[3] );
                        
    /** Signify that the data will be displacement only */
    void SetDisplacementOnly( bool disp_only );

    /** Get the displacement-only flag setting */
    bool GetDisplacementOnly() { return sdata.dispOnly; }

    /** How many doubles to read in and save (per line/node) */
    const int& GetNumDoubles() { return nDoubles; }

    /** Get the name of the materials file */
    const char* GetMaterialsFile() { return materialsFile.c_str(); }

    /** Fix for BMU solver output -- which may have incorrect setting in header */
    void SetNumberOfNodes(const unsigned long& N) { nNodes = N; }

    /** Read displacement data form ParaBMU solver output.
       \param actualNodeCount needed as that in the file is in error.
    */
    double* ReadBMUSolverDisplacementData( ifstream& fin, const unsigned long& actualNodeCount );

    /** Read displacement data from PetSc solver output. */
    double* ReadPETSCSolverDisplacementData( ifstream& fin,  map< unsigned long, double* >& mapGlobalToStrainData );

  protected:
  
    unsigned long nNodes;    ///< total number of nodes
    unsigned long nodesRead; ///< number of nodes read in (check)
    unsigned long dim[3];    ///< 3D size of block in terms of nodes
    string materialsFile;    ///< materials file for ym/pr
    
    SData sdata;    ///< struct type to hold a row of disp/strain data
    int  nDoubles;  ///< How much data we're keeping track of per table row

};


#endif

