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
 *
 * This class was adapted from ConnectedComponentImageFilter ITK Wiki example
 *  Revisions:
 *  Dec 2013 output to 16bit image format for display in ITK-Snap
 *  Jan 2014 modified to accommodate models with multiple materials (no thresholds defined)
 */


#ifndef VOXFE_CONNECTED_COMPONENTS_IMAGE_FILTER_H
#define VOXFE_CONNECTED_COMPONENTS_IMAGE_FILTER_H


#include <string>
using std::string;

#include "../VoxFEDefines.h"
#include "ReadVoxFEScript.h"

//ITK may not be available on some machines (eg. Archer), so build without if this is not set
//in CMake. Without ITK, you will not be able to import (and clean up) image data directly.
#ifdef VOXFE_USE_ITK

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkRelabelComponentImageFilter.h"

#include "itkImageSeriesReader.h"
#include "itkNumericSeriesFileNames.h"

#include "itkBMPImageIO.h"
#include "itkPNGImageIO.h"
#include "itkTIFFImageIO.h"

#include "itkRGBToLuminanceImageFilter.h"
#include "itkFixedArray.h"
#include "itkFlipImageFilter.h"

/** Read voxfe 256 grey image using ITK metaheader (mhd) file.
   \param filepath The full path and filename for the input file.
   \param imageType The image type we expect to obtain.
   \param voxelsize is used to obtain the voxel for output/idexing.
   \param bool To flip the image twice to rotate by 180 deg (usu back to original bmp orientation).

   ITK may not be available on some machines (eg. Archer), so build without if this is not set
   in CMake. Without ITK, you will not be able to import (and clean up) image data directly.
 */
template<class TImage>
itk::SmartPointer< TImage >
voxfeITKReadImageFileMHD( const char* filepath, const TImage* image, double& voxelsize, const bool flip=true );

/** Read voxfe 256 grey image using ITK from png or bmp series. These files require an rgb stream
    to be read in -- the returned image will be cast to that supplied as the secong arg.
    [Note: only tested so far on bmp images.]
   \param filepath The full path and C-type printf file specifier for the input file (eg. /path/to/files/file%03d.bmp)
   \param imageType The image type we expect to obtain.
   \param voxelsize The size of each voxel (voxfe expects cubes!)
   \param first The start of the bmp file series.
   \param last The last image of the bmp file series.
   \param incr The increment in which to read files - for voxfe this probably should always be 1.
   \param bool To flip the image twice to rotate by 180 deg (usu back to original bmp orientation).
 */
template<class TImage>
itk::SmartPointer< TImage >
voxfeITKReadImageFilePNG( const char* filepath, 
                          const TImage* image,
                          const double& voxelsize,
                          const int& first, const int& last, const int incr=1,
                          const bool flip=true );


/** Write the image to file using ITK
   \param inputImage Pointer to itkImage.
   \param filepath The full path and filename for the output file.
 */
template<class TImage>
void voxfeITKWriteToImageFileMHD( const TImage* inputImage, const char* filepath );


/** Create an image of the same size as the one supplied
   \param inputImage Pointer to itkImage.
   \param pxl The output image pixel type.
   \param init If true, set all pixels to zero
 */
template<class InputImage, class OutputImagePixelType>
itk::SmartPointer< itk::Image<OutputImagePixelType, voxfeDimension> >
voxfeITKCreateEmptyImage( const InputImage* inputImage, OutputImagePixelType pxl, bool init );


/** Perform ITK connected components (CC) routine on supplied image. Either an image file OR
    VoxFE model/script files will be generated, depending upon CH option selected.
   \brief The largest object from running CC is kept (other objects will be removed).
          The output, however, is that portion of the input image which matches.
          the largest object (label 1) in the CC image.
   \param inputImage Pointer to itkImage.
   \param clabel The intermediate CC image pixel type (unsigned, but note: an 16bit image will only
                 will only support ~64K labelled objects).
   \param CH If 'I' is given, an image is generated, 'V' will produce voxfe script/model files.
   \param filename The full path and filename for the input file ("_voxel.mhd" or "_voxel.model" will be appended).
   \param algorithm The name of the algorithm for the solver
   \param _border The border added for remodelling growth
   \return The set of labels contained in the connected image (user is responsible for deletion).
 */
template<class TImage, class CCLabelType>
set<short>*
voxfeITKProcessImage( TImage* inputImage, CCLabelType clabel,
                     const char& CH, const char* filename,
                     const char* algorithm,
                     const int _border[3]);

//some general typedefs
typedef itk::Image<unsigned char, voxfeDimension> voxfeGreyUCharImage;
typedef itk::ImageRegionIteratorWithIndex<voxfeGreyUCharImage> voxfeGreyUCharImageIterator;
typedef voxfeGreyUCharImage::Pointer voxfeGreyUCharImagePointer;

typedef itk::Image<short, voxfeDimension> voxfeGreyShortImage;
typedef itk::ImageRegionIteratorWithIndex<voxfeGreyShortImage> voxfeGreyShortImageIterator;

typedef itk::Image<unsigned short, voxfeDimension> voxfeGreyUShortImage;
typedef itk::ImageRegionIteratorWithIndex<voxfeGreyShortImage> voxfeGreyUShortImageIterator;

typedef itk::Image<unsigned long, voxfeDimension> voxfeGreyULongImage;
typedef itk::ImageRegionIteratorWithIndex<voxfeGreyULongImage> voxfeGreyULongImageIterator;

typedef itk::Index<voxfeDimension>                 voxfeIndex3D;
typedef	itk::ImageRegion<voxfeDimension>           voxfeRegion3D;
typedef itk::ImageRegion<voxfeDimension>::SizeType voxfeSize3D;

#endif //#ifdef VOXFE_USE_ITK


/** Helper split string func, courtesy of cplusplus.com */
void voxfeSplitFilename (const string& str, string& path, string& file, const char* sep="/\\" );


#endif

