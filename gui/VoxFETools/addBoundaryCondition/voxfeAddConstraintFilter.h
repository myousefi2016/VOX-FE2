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

#ifndef __voxfeAddConstraintFilter_h
#define __voxfeAddConstraintFilter_h


#include <iostream>
#include <vtkMultiBlockDataSetAlgorithm.h>

class vtkAppendFilter;
class vtkAppendPolyData;

/** \class voxfeAddConstraintFilter
 *  \brief voxfeAddConstraintFilter allows the specification of boundary condition (BC) information
 *  for the input point data.
 *
 *  The chosen point data are assumed to have been selected and extracted by the user, so that
 *  they are vtkUnstructuredGrid data. Point data created by some other means may not be
 *  supported.
 *
 */
class voxfeAddConstraintFilter : public vtkMultiBlockDataSetAlgorithm /*vtkUnstructuredGridAlgorithm*/
{
 public:

  //fixme: debug test
  static int objectcounter;

  static voxfeAddConstraintFilter *New();
  vtkTypeMacro(voxfeAddConstraintFilter, vtkMultiBlockDataSetAlgorithm);  //vtkUnstructuredGridAlgorithm

  void PrintSelf(ostream& os, vtkIndent indent);

  //Set the constraint mode
  vtkGetMacro(BoundaryCondition, int);
  vtkSetMacro(BoundaryCondition, int);  
  
  //0. VOXFE_BC_NODAL Set the axis/axes to constrain
  vtkGetVector3Macro(AxisConstraint, int);
  vtkSetVector3Macro(AxisConstraint, int);
  void SetAxisConstraintX(int _b) { AxisConstraint[0] = _b; this->Modified();  } 
  int  GetAxisConstraintX()       { return AxisConstraint[0]; }
  void SetAxisConstraintY(int _b) { AxisConstraint[1] = _b; this->Modified();  } 
  int  GetAxisConstraintY()       { return AxisConstraint[1]; }
  void SetAxisConstraintZ(int _b) { AxisConstraint[2] = _b; this->Modified();  } 
  int  GetAxisConstraintZ()       { return AxisConstraint[2]; }


  //@ 1. VOXFE_BC_FORCE_PARALLEL Set the force vector directly
  vtkGetVector3Macro(ForceVector, double);
  vtkSetVector3Macro(ForceVector, double);

  //@ 2. VOXFE_BC_FORCE_TO_POINT Set the force endpoint, and compute the vector later
  vtkGetVector3Macro(ForceEndpoint, double);
  vtkSetVector3Macro(ForceEndpoint, double);
  
  //@ St Force magnitude
  vtkGetMacro(ForceMagnitude, double);
  vtkSetMacro(ForceMagnitude, double);  

  //@ Set how the force is distributed
  vtkGetMacro(ForceDistribution, int);
  vtkSetMacro(ForceDistribution, int);  

  //@ Set the glyph scale representation
  vtkGetMacro(GlyphScale, double);
  vtkSetMacro(GlyphScale, double);
  
  //@ Set the glyph offset step (for use with slider)
  vtkGetMacro(OffsetStep, double);
  vtkSetMacro(OffsetStep, double);

  //@ Set the glyph offset, to make more visible
  vtkGetMacro(OffsetX, double);
  vtkSetMacro(OffsetX, double);

  //@ Set the glyph offset, to make more visible
  vtkGetMacro(OffsetY, double);
  vtkSetMacro(OffsetY, double);

  //@ Set the glyph offset, to make more visible
  vtkGetMacro(OffsetZ, double);
  vtkSetMacro(OffsetZ, double);

protected:

  voxfeAddConstraintFilter();   ///< Protected constructor
  ~voxfeAddConstraintFilter();  ///< Protected destructor

  int    BoundaryCondition;  ///< BC type based on SM def: 0=none, 1=nodal, 2=force/parallel, 3=force/to point
  int    AxisConstraint[3];  ///< Nodal constraint 0=free, 1=fixed
  double ForceVector[3];     ///< Force vector (unscaled)
  double ForceEndpoint[3];   ///< User-defined force end-point
  double ForceMagnitude;     ///< Scaling for force vector
  int ForceDistribution;     ///< Spread of force based on SM def: 0=node, 1=area(divided over all nodes)
  double GlyphScale;         ///< Scaling of glyph representation
  double OffsetStep;         ///< Offset step (used as multiple with slider), to nudge glyph into view
  double OffsetX;            ///< Offset slider X value
  double OffsetY;            ///< Offset slider Y value
  double OffsetZ;            ///< Offset slider Z value
  
  vtkAppendFilter* AppendUGridData;    ///< Pointer to add unstructured grid data
  vtkAppendPolyData* AppendGlyphData;  ///< Pointer to append glyph data
  bool IsDefined;            ///< Has the BC been previously defined... we might just want to change the glyphs

  /** \brief Make sure the pipeline knows what type we expect as input */
  int FillInputPortInformation( int port, vtkInformation* info );
  
  /** Generate boundary condition data */
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:

  voxfeAddConstraintFilter(const voxfeAddConstraintFilter&);  ///< Copy constructor, not implemented.
  void operator=(const voxfeAddConstraintFilter&);            ///< Assignment, not implemented.
};

#endif
