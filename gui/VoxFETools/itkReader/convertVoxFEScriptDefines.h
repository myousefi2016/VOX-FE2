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

#ifndef CONVERT_VOXFE_DEFINES_H
#define CONVERT_VOXFE_DEFINES_H

//define some offsets for standard "1mm wide" ansys models
#define ANSYS_OFFSET_X 0.5
#define ANSYS_OFFSET_Y 0.5
#define ANSYS_OFFSET_Z 0.0

#define ANSYSHEADER_LINES 10
const char* ANSYSHEADER[] = {
  "!  This is ANSYS script format, which can be read in directly by ANSYS.",
  "!  Everything after  !  on a line is ignored ",
  "          ",
  "fini	! exits from any processor and returns to the start level",
  "/clear,start	! resets the ANSYS database",
  "          ",
  "/prep7               ! starts the pre-processor",
  "          ",
  "et,1,SOLID45         ! defines element type 1 as ANSYS-type  SOLID45 (8 node brick)",
  "keyopt,1,2,1         ! reduced integration"
};


#define FELTHEADER_LINES 18
const char* FELTHEADER[] = {

"problem description\n",
"title=\"Brick FE problem\" nodes=",
"  elements=",
"  analysis=static\n",
"\nnodes\n",
"\nbrick elements\n",
"\nmaterial properties\n",
"\nforces\n",
"\nconstraints",
"fixed    Tx=c Ty=c Tz=c",
"free     Tx=u Ty=u Tz=u",
"roll_x   Tx=u Ty=c Tz=c",
"roll_y   Tx=c Ty=u Tz=c",
"roll_z	  Tx=c Ty=c Tz=u",
"roll_xy  Tx=u Ty=u Tz=c",
"roll_yz  Tx=c Ty=u Tz=u",
"roll_zx  Tx=u Ty=c Tz=u",
"\nend\n"

};


#define VTK_VOXEL_CELLTYPE 11
#define VTKHEADER_LINES 8  //last entries here are if we want to add strains etc
const char* VTKHEADER[] = {

"# vtk DataFile Version 3.1\n",
"Created by convertVoxFEScript\n",
"ASCII\n\n",
"DATASET UNSTRUCTURED_GRID\n",
"POINTS ",
" float \n",
"CELLS ",
"CELL_TYPES ",
"CELL_DATA ",
"SCALARS ",
" int\n",
"LOOKUP_TABLE default\n"

};




#endif

