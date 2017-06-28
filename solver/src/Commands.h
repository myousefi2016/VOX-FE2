//
//  Commands.h
//  
//
//  Created by N. Banglawala
//
//

#ifndef GUARD_COMMANDS_H
#define GUARD_COMMANDS_H


// filelengths etc.
#define MAX_COMMAND_LENGTH 256
#define MAX_DATA_LENGTH 256
#define MAX_FILENAME_LENGTH 256

// Rank of the master process
enum { MASTER=0 };

//======================== COMMANDS =====================================

// compatibility with old solver
enum { CMD_SET_VOXEL_SIZE, CMD_LOAD_MATERIALS, CMD_LOAD_MODEL, CMD_SET_TOLERANCE, CMD_SET_MAX_ITER, CMD_SET_ALGORITHM_FEA, CMD_SELECTION_OF_NODES, CMD_LOAD_CONSTRAINTS, CMD_SELECT_NODE_3D, CMD_PRESERVE_NODE_3D, CMD_COMPUTE_SED, CMD_SOLVE, CMD_PRINT_DISPLACEMENTS, CMD_FINISH };


#endif /*( GUARD_COMMANDS_H )*/
