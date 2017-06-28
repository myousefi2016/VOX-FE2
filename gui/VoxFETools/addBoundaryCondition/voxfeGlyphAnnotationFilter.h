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

#ifndef __voxfeGlyphAnnotationFilter_h
#define __voxfeGlyphAnnotationFilter_h

#include "vtkSmartPointer.h" // compiler errors if this is forward declared
#include <vtkPolyDataAlgorithm.h>
#include "vtkExecutive.h"

class vtkPolyData;
class vtkInformation;
class vtkInformationVector;

#include "vtkMath.h"

/** \class voxfeGlyphAnnotationFilter
 *  \brief voxfeGlyphAnnotationFilter adds glyph annotation to points in the boundary condition (BC).
 *
 *  This has to be separated from the input BC data as adding the glyph geometry 
 *  may get mixed up with real model geometry. Also, the user may wish to be able to
 *  switch the glyphs off separately.
 */
class voxfeGlyphAnnotationFilter : public vtkPolyDataAlgorithm
{
 public:
 
  static voxfeGlyphAnnotationFilter *New();
  vtkTypeMacro(voxfeGlyphAnnotationFilter, vtkPolyDataAlgorithm);
  
  //glyph scale representation
  vtkGetMacro(GlyphScale, double);
  vtkSetMacro(GlyphScale, double);
  
 protected:
  voxfeGlyphAnnotationFilter();
  ~voxfeGlyphAnnotationFilter();
  
  double GlyphScale;         ///< Scaling of glyph representation  

  /** Make sure the pipeline knows what type we expect as input */
  int FillInputPortInformation( int port, vtkInformation* info );

  /** Generate output data */
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *); //the function that makes this class 
                                                                                      //work with the vtk pipeline

};

#endif

