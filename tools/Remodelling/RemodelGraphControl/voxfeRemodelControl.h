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

#ifndef VOXFE_REMODEL_CONTROL_GRAPH_H
#define VOXFE_REMODEL_CONTROL_GRAPH_H

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
#include <sstream>

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
using std::stringstream;

#include "strainReader/voxfeRecomputeStrainFilter.h"
#include "itkReader/ReadVoxFEScript.h"
#include "RunningStat.h"


#define VOXFE_METIS_MAX_BUFF_SIZE 4096  // needs to be able to contain up to 26 long integers

#define VOXFE_GET_N26_INDEX(vox_diff) ( (vox_diff[2]+1)*9 + (vox_diff[1]+1)*3 + (vox_diff[0]+1) )

#define VOXFE_GET_N26_XYZ(index,xyz) {\
    unsigned long rem;         \
    xyz[2] = (index / 9) - 1;  \
    rem    = index % 9;        \
    xyz[1] = (rem / 3) - 1;    \
    xyz[0] = (rem % 3) - 1;    \
}

#define VOXFE_NUM_HISTOGRAM_BINS 100000

/** \class voxfeRemodelControl
    \brief This class maintains a graph of voxel connectivity by global
    node ID (also giving voxel indices) and bits identifying the
    adjacent 26 neighbours.


*/
class voxfeRemodelControl {

  public:

    typedef unsigned int uint;
    typedef rvsLongIntegerType voxfeRemodelIntType;      ///< The type used to write the binary nodemap file
    typedef uint LabelType;
    //typedef rvsIntegerType voxfeRemodelIntType;

    struct VertexData {
      uint connectivity;       ///< Connectivity of this voxel, each bit = one neighbour (max 26 bits)
      short status;            ///< Hsas voxel been added or removed
      unsigned char material;  ///< Material ID
      unsigned char remodel;   ///< Allow this voxel to be removed or added to (1 = yes)
      double SED;              ///< Strain energy density of this voxel
      LabelType component;     ///< For connected component label

      VertexData();
    };

    typedef unsigned long IDType;
    typedef map<IDType, VertexData*> GraphType;   ///< Mapping from global ID# (equivalent to Node 0 of this voxel) to data,
    typedef map<IDType, VertexData*>::iterator GraphTypeIterator;

    typedef set<IDType> IDSet;   ///< Set of element IDs
    typedef set<IDType>::iterator IDSetIterator;

    typedef set<LabelType> LabelSet;   ///< Set of Connected Component labels
    typedef set<LabelType>::iterator LabelSetIterator;

    typedef map<IDType, double*> NodeToDisplacementMapType;
    typedef map<IDType, double*>::iterator NodeToDisplacementMapIterator;

    voxfeRemodelControl(const char* scriptFile,
                        const char* elementMapFile,
                        const char* metisGraphFile );

    ~voxfeRemodelControl();

    bool IsOnSurface( const IDType& vertexID );  ///< yes if < 26 neighbours

    /** Get neighbours (existing) and non-neighbours (not present) of this voxel */
    bool GetNeighbours( const IDType& vertexID, IDSet& neighbours, IDSet& non_neighbours );

    /** Get only face-connected neighbours of this voxel (6-connectivity) */
    bool GetFaceNeighbours( const IDType& vertexID, IDSet& neighbours );

    /** Get the existing neighbours surrounding the ID Node (N0) of this voxel */
    bool GetNodeNeighbours( const IDType& vertexID, IDSet& neighbours );

    /** Input diplacements and compute strains/energy densities.
        If successful, remodelling files will also be generated.
    */
    bool SetDisplacements( const char* dispFile, bool withLocation,
                           bool removeDisconnectedElements=true, bool use26NeighConnectivity=false );

    /** Set thresholds used to determine remodelling growth */
    void SetThresholds( const double& lowerSEDThreshold, const double& upperSEDThreshold, bool adaptive=false );

    /** Dump the graph contents */
    void PrintGraph( std::ostream& out, bool withNeighbours );

    /** Write out a new model file with the graph at the current update step */
    void WriteModel(bool removeDisconnectedElements);

    /** Perform connected components labelling, keeping only the largest component
        if selectLargest is true (other components are set to zero)
    */
    void SetConnectedComponents( bool use26NeighConnectivity, bool selectLargest );

    /** Dump the Equivalence table (connected components) */
    void PrintEquivTable( const char* tableFile );


  protected:

    voxfeRemodelControl();  //not implemented

    GraphType Graph;                             ///< The graph of elements
    ModelMetaData<voxfeRemodelIntType>* Model;   ///< Link to VoxFE model metadata
    uint* NeighbourData;           ///< Store the connectivity bit pattern for neighbours
    double* Displacement;          ///< To store diplacement data
    double Threshold[2];           ///< Thresholds (lower & upper) to determine remodelling
    bool AdaptiveThresholding;     ///< Are the thresholds supplied hard limits ie particular strain values
                                   ///< or are they 'adaptive' eg. 0.05 and 0.95 would indicate the bottom 5%
                                   ///< and the highest 5% for remodelling
    int UpdateStep;                ///< Number of updates to the graph
    const char* ScriptFile;        ///< Retain the name of the original script file

    RunningStat SED_stats;        ///< Compute running statistics for SED
    IDSet voxels_to_add_onto;     ///< Remodelling growth voxels
    IDSet voxels_to_remove;       ///< Remodelling shrinkage

    //fixme: should we zero these at each step?
    IDType NumVoxelsRemoved;      ///< Count the voxels actually removed (status<0) to get the 'actual' size of the graph
    IDType NumVoxelsAdded;        ///< Count the voxels actually added, not those we're adding onto

    vector< LabelSet > EquivTable;  ///< Table of equivalences for Connected Components

    NodeToDisplacementMapType NodeToDisplacementMap;  ///< Map from Node IDs to displacements (to compute strains)

    /** Add graph connections and create vertex if not extant .
        Inputs here are metis 1-based element numbers, 1,2,..., nEl.
        Called from constructor.
    */
    bool CreateConnections( const IDType& metisElementID,
                         const vector<IDType>& metisNeighIDs,
                         const voxfeRemodelIntType* elementIDList,
                         const unsigned char* materialList );

    /** Load the element (dual) graph obtained from "m2gmetis -gtype=dual"
        conversion of node connectivity.

        Called from constructor.
    */
    bool CreateGraphFromMetis( const char* metisGraphFile,
                          const voxfeRemodelIntType* element,
                          const unsigned char* material  );

    /** Forbid elements with boundary conditions or loads attached from remodelling
        Called from constructor.
    */
    void CreateConstraintRestrictions();


    /** Helper func to assist processing metis graph data */
    void SplitLine(const string& s, char delim, vector<IDType>& elems);

    /** Set the bits defining neighbouring voxels */
    void CreateNeighbourDataBits();

    /** Initialize the mapping from nodes to displacements.
        UpdateNodeToDisplacementMap() must be called to finalize the map, once
        displacements have been read.
    */
    void InializeNodeToDisplacementMap();

    /** Finalize the NodeToDisplacementMap() */
    void UpdateNodeToDisplacementMap();

    /** Compute strain energies across the graph.
        \return false, if nodes cannot be mapped
    */
    bool ComputeSED(bool surface_only);

    /** Determine the voxels to add onto or remove, depending on SED thresholds. */
    void SelectBySED(bool surface_only);

    /** Assemble voxel displacement vector */
    void SetVoxelDisplacements( GraphTypeIterator& git, voxfeRecomputeStrainFilter::uVector& u );

     /** Use thresholded voxels to determine remodelling */
    void Remodel();

    /** Initial pass to assign component labels with neighbour filling to lowest label */
    void AssignComponentLabels( bool use26NeighConnectivity );

    /** Work recursively through the EquivTable to find the lowest (parent) label,
        effectively flattening the table
     */
    void FlattenEquivalenceTable();

    /** Work recursively through the EquivTable to find all parent labels */
    void FindAllParentLabels( const LabelType& label, LabelSet& parents );

    /** Work through the component labels of elems, and assign the elements to the
        local parent/lowest label, with updates to EquivTable (highest to lowest)
        \param elems should include the current and neighbouring voxels
    */
    bool SetToLocalParentLabel( IDSet& elems );

    /** Reset the equivalence table, to initialize connected components gathering */
    void ResetComponents();

    void GetNextLabel( LabelType& label );               ///< Get the next label available

    //const LabelType& SetToLowestLabel( IDSet& labels );  ///< Find lowest label and set within the EquivTable

    /** Check whether a voxel element is in the graph and has status >= 0 (ie has
        not been 'removed' by remodelling/connections
    */
    inline bool ElementIsValid( const IDType& vertexID, GraphTypeIterator& it ) {
      it = this->Graph.find( vertexID );
      if(  ( it != this->Graph.end() ) && ( it->second->status >= 0 ) )
        return true;

      return false;
    }

    enum NeighbourDataConstants {
      //Assume : L=left      (x-1), N=null (0), R=right    (x+1),
      //         P=posterior (y-1), N=null (0), A=anterior (y+1),
      //         I=inferior  (z-1), N=null (0), S=superior (z+1)

      ND_NO_NEIGH   = 0x00000000,
      ND_NUM_FACES  = 6,
      ND_NODE_NEIGH = 13,
      ND_NUM_NEIGH  = 26,
      ND_BLOCK_SIZE = 27,
      ND_ARRAY_PRESERVE = 30,  ///< Identify voxels around this vertex N0
      ND_ARRAY_LAST = 31,  ///< Collect the largest 'full' neighbour value
      ND_ARRAY_SIZE = 32,

#if 0
      ix_LPI         = 0x000000001,   // -1, -1, -1
      ix_NPI         = (ix_LPI <<  1),   //  0, -1, -1
      ix_RPI         = (ix_LPI <<  2),   //  1, -1, -1

      ix_LNI         = (ix_LPI <<  3),   // -1,  0, -1
      ix_NNI         = (ix_LPI <<  4),   //  0,  0, -1
      ix_RNI         = (ix_LPI <<  5),   //  1,  0, -1

      ix_LAI         = (ix_LPI <<  6),   // -1,  1, -1
      ix_NAI         = (ix_LPI <<  7),   //  0,  1, -1
      ix_RAI         = (ix_LPI <<  8),   //  1,  1, -1
      //========================================
      ix_LPN         = (ix_LPI <<  9),   // -1, -1,  0
      ix_NPN         = (ix_LPI << 10),   //  0, -1,  0
      ix_RPN         = (ix_LPI << 11),   //  1, -1,  0

      ix_LNN         = (ix_LPI << 12),   // -1,  0,  0
   /* NNN         = (LPI << 13),   //  0,  0,  0    not needed!  */
      ix_RNN         = (ix_LPI << 14),   //  1,  0,  0

      ix_LAN         = (ix_LPI << 15),   // -1,  1,  0
      ix_NAN         = (ix_LPI << 16),   //  0,  1,  0
      ix_RAN         = (ix_LPI << 17),   //  1,  1,  0
      //========================================
      ix_LPS         = (ix_LPI << 18),   // -1, -1,  1
      ix_NPS         = (ix_LPI << 19),   //  0, -1,  1
      ix_RPS         = (ix_LPI << 20),   //  1, -1,  1

      ix_LNS         = (ix_LPI << 21),   // -1,  0,  1
      ix_NNS         = (ix_LPI << 22),   //  0,  0,  1
      ix_RNS         = (ix_LPI << 23),   //  1,  0,  1

      ix_LAS         = (ix_LPI << 24),   // -1,  1,  1
      ix_NAS         = (ix_LPI << 25),   //  0,  1,  1
      ix_RAS         = (ix_LPI << 26),   //  1,  1,  1
#endif

      ND_CENTRE_VOXEL  = ( 1 << 13 ),
      ND_ALL_NEIGH26   = ( ( 1 << 27 ) - 1 ) - ( 1 << 13 ),
      ND_ALL_NEIGH     = 0xFFFFFFFF
    };










};

#endif
