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

#ifndef __voxfeOutputScriptFilter_h
#define __voxfeOutputScriptFilter_h

#include "vtkSmartPointer.h" // compiler errors if this is forward declared
#include <vtkMultiBlockDataSetAlgorithm.h>
#include "vtkExecutive.h"

#include "../VoxFEDefines.h"

class vtkUnstructuredGrid;
class vtkInformation;
class vtkInformationVector;

#include "vtkMath.h"

/** \class voxfeOutputScriptFilter
 *  \brief voxfeOutputScriptFilter adds in boundary conditions (BCs) for the voxfe solver script.
 *
 *  The BCs, added on port 1 via the change input dialog, are actually written to the file:
 *  "constraints.txt". The supplied vtkMultiBlockDataSet object is shallow copied through.
 */
class voxfeOutputScriptFilter : public vtkMultiBlockDataSetAlgorithm
{
 public:
 
  /** VTK standard New call */
  static voxfeOutputScriptFilter *New();

  vtkTypeMacro(voxfeOutputScriptFilter, vtkMultiBlockDataSetAlgorithm);

  void PrintSelf(ostream &os, vtkIndent indent);
  
  /** Set the main mesh object input */
  void SetInput0Connection(vtkAlgorithmOutput* algOutput);

  /** Set one or more boundary condition input */
  void SetInput1Connection(int port, vtkAlgorithmOutput* algOutput);
  
 protected:

  /** Constructor defines port settings */
  voxfeOutputScriptFilter();

  ~voxfeOutputScriptFilter();

  /** Make sure the pipeline knows what type we expect as input */
  int FillInputPortInformation( int port, vtkInformation* info );

  /** Generate output (constraint file) */
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *); //the function that makes this class 
                                                                                      //work with the vtk pipeline

  /**  Get the BC data from the supplied multiblock data.
   *
   *   The input data is wrapped in various multiblock structures - this is a helper function
   *   to get to the actual boundary condition.
   */
  vtkUnstructuredGrid* GetBoundaryConditionFromUserSelection( vtkInformation* inputBC );

  int round(double number); ///< Rounding function as per http://www.codeproject.com/Articles/58289/C-Round-Function

};

#endif

