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


#ifndef __voxfeAddConstraintFilterDefines_h
#define __voxfeAddConstraintFilterDefines_h

#define VOXFE_CONSTRAINT_ATTRIBUTE "ConditionType"
#define VOXFE_VECTOR_ATTRIBUTE     "ConditionVector"
#define VOXFE_PANEL_PARAMETER      "BoundaryCondition"

enum BoundaryCondition {
  VOXFE_BC_NONE,
  VOXFE_BC_NODAL,
  VOXFE_BC_FORCE_PARALLEL,
  VOXFE_BC_FORCE_TO_POINT,
  VOXFE_BC_NO_REMODEL
};

#endif
