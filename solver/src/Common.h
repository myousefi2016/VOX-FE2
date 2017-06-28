//
//  Common.h
//  
//
//  Created by N. Banglawala on 05/02/2015.
//
//

// Common definitions etc.

#ifndef GUARD_COMMON_H
#define GUARD_COMMON_H

//======================= MPI / Petsc HEADER FILES =====================

#include "mpi.h"
#include "petscdmda.h"
#include "petscsys.h"
#include "petscksp.h"
#include "petsctime.h"

//=======================  DEFINITIONS  =================================

// SIGN function
#define SIGN(x) ((x > 0) - (x < 0))

#define DOF_3D              3   // Number of degrees of freedom per node
#define NUM_NEIGHBRS        27  // Max number of neighbours (element or nodes)
#define NUM_TERMS           81  // DOF_3D x NUM_NEIGHBRS
#define NODES_PER_ELEMENT   8

//======================== TYPEDEF UNITS ================================

typedef unsigned int        xyzType;
typedef unsigned int        midxType;
typedef uint64_t            idxType;

//======================= ENUMS =========================================

#endif /**  (defined GUARD_COMMON_H) */

