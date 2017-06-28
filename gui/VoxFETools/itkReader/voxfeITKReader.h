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

#ifndef __voxfeITKReader_h
#define __voxfeITKReader_h

#include "../VoxFEDefines.h"
#include "vtkUnstructuredGrid.h"
#include "vtkMultiBlockDataSetAlgorithm.h"

#include <memory>

/** \class voxfeITKReader
 *  \brief voxfeITKReader reads in grey-scale image data, presumed to be segmented, and generates a
 *  vtkMultiBlockDataSet object.
 *
 *  The ITK library is used to read in the images and sanity checks are performed for (1) the sequence of
 *  images (results are output to a log file), (2) connected components, to remove small objects (assumed
 *  to be noise) and (3) all detected labels are allocated to a group (if no group is given, the label is assumed
 *  to be the group).
 *
 *  If image data are to be read, then this is specified in a ".voxfe" group file.
 *  The image data are converted to a voxfe script file format (essentially an unstructured grid)
 *  and this is read in to create vtk voxel data. If a voxfe script (".script") is initially supplied,
 *  this is used as is, without conversion of any image data.
 *
 *  The GenerateSurfaceMesh option is enabled by default and this removes
 *  internal geometry from each group as defined in the ".voxfe" file (if this is supplied).
 *  That is, each group will be represented as a block with a surface mesh within the resulting
 *  multi-block. This allows point selections on the surface of each block (to add boundary conditions).
 *  However, for visualization, turning GenerateSurfaceMesh off may be preferable, so that internal
 *  strain etc can be examined.
 *
 *  Currently the user must write the ".voxfe" file.
 */
class voxfeITKReader : public vtkMultiBlockDataSetAlgorithm
{
public:

  enum _SolverAlgorithm {
    KSPCG_PCJACOBI,
    BMU_JAPCG
  };

  vtkTypeMacro(voxfeITKReader,vtkMultiBlockDataSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  static voxfeITKReader *New();

  //@ Specify file name of the .voxfe file.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  //@ Specify file name of the strain data file.
  vtkSetStringMacro(DisplacementDataFile);
  vtkGetStringMacro(DisplacementDataFile);

  //@ Whether to remove internal geometry (needed for strain display though!)
  vtkGetMacro(GenerateSurfaceMesh, bool);
  vtkSetMacro(GenerateSurfaceMesh, bool);

  //@ Define the input image type
  vtkGetMacro(InputImageType, int);
  vtkSetMacro(InputImageType, int);

  //@ Whether to reverse the stack of images according to number ing of 'z' axis
  vtkGetMacro(ImageFlip, bool);
  vtkSetMacro(ImageFlip, bool);

  //@ A printf style image file specifier
  vtkSetStringMacro(FileSpecifier);
  vtkGetStringMacro(FileSpecifier);

  //@ The range of files to be read in (All matching files will be read if not matching; not needed for 3D image types)
  vtkSetVector3Macro(FileRange, int);
  vtkGetVector3Macro(FileRange, int);

  //@ Should we create a metis file, defining elements by their node numbers
  vtkGetMacro(CreateMetisFile, bool);
  vtkSetMacro(CreateMetisFile, bool);

  //@ Should we create a Felt file, for use with Felt library (felt.sourceforge.net)
  vtkGetMacro(CreateFeltFile, bool);
  vtkSetMacro(CreateFeltFile, bool);

  //@Set the solver algorithm
  vtkGetMacro(SolverAlgorithm, int);
  vtkSetMacro(SolverAlgorithm, int);

  //@ Define default material props
  vtkSetVector2Macro(DefaultBoneMaterial, double);
  vtkGetVector2Macro(DefaultBoneMaterial, double);
  
  //@ Define default material props
  vtkSetVector2Macro(DefaultToothMaterial, double);
  vtkGetVector2Macro(DefaultToothMaterial, double);

  //@ Define default material props
  vtkSetVector2Macro(DefaultMarrowMaterial, double);
  vtkGetVector2Macro(DefaultMarrowMaterial, double);

  //@ Define default material props
  vtkSetVector2Macro(DefaultChitinMaterial, double);
  vtkGetVector2Macro(DefaultChitinMaterial, double);

  //@ Define default material props
  vtkSetVector2Macro(DefaultMetalMaterial, double);
  vtkGetVector2Macro(DefaultMetalMaterial, double);

protected:

  /** Constructor -- define vtk port properties */
  voxfeITKReader();

  /** Destructor */
  ~voxfeITKReader(){}

  /** Assign vtk port info */
  int FillOutputPortInformation(int port, vtkInformation* info);

  /** Produce output data, in this case mainly the VoxFE script/model file is created and read
   *  to create a VTK unstructured grid (in a multiblock)
   *  */
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  /** How the data is created */
  enum voxfe_vtk_file_created {
    VOXFE_VTK_FILE_NONE,
    VOXFE_VTK_FILE_INMEMORY
  };

  /** Possible input image types */
  enum voxfe_input_image_type {
    VOXFE_ITK_MHD,
    VOXFE_BMP_SERIES,
	VOXFE_TIFF_SERIES,
	VOXFE_PNG_SERIES,
    VOXFE_SCRIPT_ONLY
  };

private:

  voxfeITKReader(const voxfeITKReader&);  ///< Private constructor, not implemented.
  void operator=(const voxfeITKReader&);  ///< Assignment operator, not implemented.

  char*  FileName;                 ///< Script/.voxfe file name
  char*  DisplacementDataFile;     ///< Displacement data file name
  bool   GenerateSurfaceMesh;      ///< Remove internal geometry flag
  int    InputImageType;           ///< The type of input image
  bool   ImageFlip;                ///< Reverse the z-stack order of images?
  char*  FileSpecifier;            ///< printf specifier for importing 2D images
  int    FileRange[3];             ///< The range of 2D images to import
  double DefaultBoneMaterial[2];   ///< Youngs Modulus/Poissons Ratio for this material
  double DefaultToothMaterial[2];  ///< Youngs Modulus/Poissons Ratio for this material
  double DefaultMarrowMaterial[2];  ///< Youngs Modulus/Poissons Ratio for this material
  double DefaultChitinMaterial[2];  ///< Youngs Modulus/Poissons Ratio for this material
  double DefaultMetalMaterial[2];   ///< Youngs Modulus/Poissons Ratio for this material
  bool   CreateMetisFile ;          ///< Create a metis file (towards remodelling)
  bool   CreateFeltFile  ;          ///< Create a Felt mesh file
  int SolverAlgorithm    ;          ///< Choose the solver algorithm

};


#endif
