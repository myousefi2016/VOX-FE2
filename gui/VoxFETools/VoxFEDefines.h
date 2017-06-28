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


#ifndef __voxfeDefines_h
#define __voxfeDefines_h

/**
  Include/resolve all global defines here
*/

//file constants
#define VOXFE_VOXFEFILE "voxfefile"
#define VOXFE_VOXELSIZE "voxelsize"
#define VOXFE_DEFAULT_CONSTRAINT_FILE   "constraints.txt"
#define VOXFE_MATERIALS_FILE            "materials.txt"
#define VOXFE_NODEMAP_FILE_EXTN ".nodemap"

const int voxfeVoxelOffset=1;  //Voxel indices assumed to be 1-based
const unsigned int voxfeDimension         = 3;

//FE constants
#define LSM_BUFF_LENGTH 32
#define ZERO_TOL 1e-6
#define VOXFE_COMMENT_LINE_LENGTH 256
#define VOXFE_YOUNGSMODULUS  1.7e4   //assuming N/mm^2
#define VOXFE_POISSONSRATIO  0.3
#define VOXFE_MAX_MATERIALS 256
#define NODESPERELEMENT 8
#define HALFNODESPERELEMENT 4

#include "addBoundaryCondition/voxfeAddConstraintFilterDefines.h"

#endif
