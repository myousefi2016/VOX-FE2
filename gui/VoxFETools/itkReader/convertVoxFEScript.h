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

/*
 *
 * Script to convert to ansys input file from vox-fe data
 *
 * g++ -c ReadVoxfeScript.cpp
 * g++ -o convertVoxFEScript convertVoxFEScript.cpp ReadVoxFEScript.o
 *
 * //some useful awk scripts for Snap output (edited for num voxels)
 * awk '{ if( NF > 1 ) print NR-2, " 1 ", $0; else print $0; }' biomed-cube_all.xyz > biomed-cube_all.mod
 * awk '{ if( ( NF > 1 ) && ( $3 == 0 ) ) print "    SELECT_NODE_3D ", $0, " 0 0 0 1 1 1"; }' biomed-cube_all.xyz > constraints.txt
 *
 * //or to examine the ansys script
 * awk 'BEGIN{ FS=","; }{ if( ($1=="n") && ($3==0) && ($5==0) ) print $0; }' circ.script.ansys > test.txt
 */



#ifndef CONVERT_TO_ANSYS_SCRIPT_H
#define CONVERT_TO_ANSYS_SCRIPT_H

#include "convertVoxFEScriptDefines.h"

#include "ReadVoxFEScript.h"

#include "stdlib.h"
#include "math.h"
#include <set>
using std::set;

/** \brief A set of routines to convert VoxFE script files to other input types:
 * Ansys, Felt, VTK and Metis
 *
 */

/** Conversion type */
enum ConversionType {
  NO_CONVERSION,
  ANSYS,
  METIS,
  FELT,
  VTK
};


//================================ FUNCTION HEADERS =================================================

/** Write out Ansys mesh file */
template<class T>
bool WriteAnsysSolid45File( const char* opfilename,
			    const map<T, T>& nodes,
          const Element<T>* element,
			    const T& n_elements,
			    const ModelMetaData<T>& cmdata );


/** Write out .flt file for input to FELT library */
template<class T>
bool WriteFeltMeshFile( const char* opfilename,
                        const map<T, T>& nodes,
                        const Element<T>* element,
                        const T& n_elements,
			const ModelMetaData<T>& model );

/** Write out .metis file for input to Metis graph partition */
template<class T>
bool WriteMetisMeshFile( const char* opfilename,
                        const map<T, T>& nodes,
                        const Element<T>* element,
                        const T& n_elements,
			const ModelMetaData<T>& model );

/** Write out .vtk (legacy) file for input to VTK/paraview */
template<class T>
bool WriteVTKMeshFile( const char* opfilename,
                        const map<T, T>& nodes,
                        const Element<T>* element,
                        const T& n_elements,
			const ModelMetaData<T>& model );



/** Collect nodes (and elements) at the top and bottom of the mesh for partitioning with ParFE */
template<class T>
void CollectTopAndBottomNodes( const T  minmaxZ[2],
                               const map<T, T>& nodes,
                               const Element<T>* element,
                               const T& n_elements,
                               const ModelMetaData<T>& cmdata,
                               set<T>& topNodes, set<T>& topElements,
                               set<T>& botNodes, set<T>& botElements );

/** Helper function to output top node set etc */
template<class T>
void PrintNodeSet( ofstream& fout, const set<T>& elements );

#endif

