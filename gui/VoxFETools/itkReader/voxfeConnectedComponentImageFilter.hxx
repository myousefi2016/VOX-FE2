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

#ifndef VOXFE_CONNECTED_COMPONENTS_IMAGE_FILTER_HXX
#define VOXFE_CONNECTED_COMPONENTS_IMAGE_FILTER_HXX

#include "voxfeConnectedComponentImageDefines.h"
#include "voxfeConnectedComponentImageFilter.h"

#include <fstream>
using std::ofstream;

#ifdef VOXFE_USE_ITK

#include "itksys/SystemTools.hxx"
#include "itksys/Directory.hxx"
#include "itksys/Glob.hxx"

// =====================================================================================================

template<class TImage>
itk::SmartPointer< TImage >
voxfeITKReadImageFileMHD( const char* filepath, const TImage* image,
                          double& voxelsize, const bool flip ) {

  //create a log file to record errors -- windows users in partic will not see cerr
  string logFile = string(filepath) + string(voxfeCCLogFileAppendString);
  ofstream logOut(logFile.c_str());

  typedef itk::ImageFileReader<TImage> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(filepath);

  typename TImage::Pointer newImage = 0;

  if( flip ) {

    //flip the image twice to rotate it back upwards
    itk::FixedArray<bool, voxfeDimension> flipAxes;

    flipAxes[0] = false;
    flipAxes[1] = true;
    flipAxes[2] = true;

    typedef itk::FlipImageFilter<TImage> FlipImageFilterType;
    typename FlipImageFilterType::Pointer flipFilter = FlipImageFilterType::New ();
    flipFilter->SetInput( reader->GetOutput() );
    flipFilter->SetFlipAxes(flipAxes);

    try
    {
      flipFilter->Update();
    }
    catch (itk::ExceptionObject &ex)
    {
      logOut << ex << endl;
      logOut.close();
      return newImage;
    }

    newImage = flipFilter->GetOutput();
  }
  else {

    try
    {
      reader->Update();
    }
    catch (itk::ExceptionObject &ex)
    {
      logOut << ex << endl;
      logOut.close();
      return newImage;
    }

    newImage = reader->GetOutput();
  }

  logOut.close();

  //obtain the voxelsize of our new image -- check cubic??
  const typename TImage::SpacingType& spacing = newImage->GetSpacing();
  voxelsize = spacing[0];

  return newImage;
}


template<class TImage>
itk::SmartPointer< TImage >
voxfeITKReadImageFilePNG( const char* filepath, 
                          const TImage* image,
                          const double& voxelsize,
                          const int& first, const int& last, const int incr,
                          const bool flip ) {

  //pull out the file extension //fixme: check...?
  string imageFName, imageFExt;
  voxfeSplitFilename( string(filepath), imageFName, imageFExt, "." );

  //create a log file to record errors -- windows users in partic will not see cerr
  string logFile(filepath);

  //above may well contain '%', so replace
  for( string::iterator sit=logFile.begin(); sit!=logFile.end(); ++sit )
    if( *sit == '%' ) *sit = '_';

  logFile += string(voxfeCCLogFileAppendString);
  ofstream logOut(logFile.c_str());

  //code here much borrowed from ImageSeriesReadWrite.cxx
  typedef unsigned char ComponentType;

#ifdef VOXFE_USE_RGB_IMAGE
  typedef itk::RGBPixel<ComponentType> RGBPixelType;
  typedef itk::Image< RGBPixelType, voxfeDimension >  RGBImageType;
  typedef itk::ImageSeriesReader< RGBImageType >  RGBReaderType;
  typedef itk::ImageSeriesReader< RGBImageType >  RGBReaderType;
  RGBReaderType::Pointer reader = RGBReaderType::New();
#else
  typedef itk::Image< ComponentType, voxfeDimension >  bit8ImageType;
  typedef itk::ImageSeriesReader< bit8ImageType >  bit8ReaderType;
  bit8ReaderType::Pointer reader = bit8ReaderType::New();
#endif

  typedef itk::NumericSeriesFileNames    NameGeneratorType;
  NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();

  string path, specifier;
  voxfeSplitFilename (string(filepath), path, specifier);
  nameGenerator->SetSeriesFormat( specifier.c_str() );
  nameGenerator->SetStartIndex( (unsigned int)first );
  nameGenerator->SetEndIndex( (unsigned int)last );
  nameGenerator->SetIncrementIndex( (unsigned int)incr );

  //=================== debug check =========================
  logOut << "\n\n***  voxfeITKReadImageFilePNG has spec: " << specifier
	     << " from: "       << filepath
	     << " (with path: " << path
	     << ")  ***\n\n";

  itksys::Glob fileGlobber;
  fileGlobber.RecurseOff();
  string globExpression = path;
  globExpression += string("/*.") + imageFExt; //"/*.bmp";
  itksys::SystemTools::ConvertToUnixSlashes(globExpression);
  fileGlobber.FindFiles(globExpression);

  //construct a set to compare with the supplied list
  set<string,voxfeStringLessThan> fileSet;
  vector<string>& globbedFiles = fileGlobber.GetFiles();
  for(unsigned int i = 0; i < globbedFiles.size(); ++i)
  {
    string s = itksys::SystemTools::GetFilenameName(globbedFiles[i]);
    fileSet.insert(s);

    logOut << globbedFiles[i] << "  [" << s << "]" << endl;   //debug
  }

  int fileCount = 0;
  vector< string > fileNames = nameGenerator->GetFileNames();
  vector< unsigned char > fileFound;
  bool seqError = false;
  for(unsigned int i = 0; i < fileNames.size(); ++i)
  {
    if( fileSet.find( fileNames[i] ) == fileSet.end() ) {
      logOut << "*** FILE: " << fileNames[i] << " WAS NOT FOUND ***" << endl;
      fileFound.push_back(0);
      seqError = true;
    }
    else {
      fileFound.push_back(1);
      fileCount++;
    }
  }

  //if there's an error use the globbed files
  vector<string> sortedGlobbedFiles;
  if( seqError ) {

    logOut << "\n\nSelecting files for inclusion:\n\n";  //debug

    set<string,voxfeStringLessThan>::iterator nameSetIt =  fileSet.begin();
    for( ; nameSetIt!=fileSet.end(); ++nameSetIt ) {
      sortedGlobbedFiles.push_back( *nameSetIt );

      logOut << *nameSetIt << endl;  //debug
    }
  }

  //determine image extn
  if( (imageFExt == "bmp") || (imageFExt == "BMP") )
	reader->SetImageIO( itk::BMPImageIO::New() );
  else if( (imageFExt == "png") || (imageFExt == "PNG") )
    reader->SetImageIO( itk::PNGImageIO::New() );
  else if( (imageFExt == "tif") || (imageFExt == "TIF") || (imageFExt == "tiff") || (imageFExt == "TIFF") )
	reader->SetImageIO( itk::TIFFImageIO::New() );

  if( seqError )
    reader->SetFileNames( sortedGlobbedFiles );  //try to guess the correct files
  else
    reader->SetFileNames( fileNames );           //use as intended

  try {
    reader->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    logOut << "\n================================================"   << endl;
    logOut << "*** Exception caught for reading image series ****" << endl;
    logOut << "ERROR: " << err << endl;
    logOut << "\n================================================"   << endl;
    logOut.close();

    return 0;  //still null here
  }

#ifdef VOXFE_USE_RGB_IMAGE
  //unsigned char pxl; voxfeGreyUCharImage::Pointer
  typename TImage::PixelType pxl = 0;
  typename TImage::Pointer newImage = voxfeITKCreateEmptyImage( reader->GetOutput(), pxl, false );
  typedef itk::ImageRegionIteratorWithIndex<TImage> TImageIterator;

  //extract streams ourselves
  typedef itk::ImageRegionIteratorWithIndex<RGBImageType> RGBImageTypeIterator;
  RGBImageTypeIterator itRGB( reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion());

  TImageIterator itTarget( newImage, newImage->GetLargestPossibleRegion() );

  for( itRGB.GoToBegin(), itTarget.GoToBegin(); !itRGB.IsAtEnd(); ++itRGB, ++itTarget )  {

#if 0
    const RGBPixelType& rgbpixel = itRGB.Value();
    unsigned char& tpixel        = itTarget.Value();
    /*voxfeIndex3D idx = itRGB.GetIndex();
    cout << "RGB pixel: " << idx[0] << "  " << idx[1] << "  " << idx[2] << "      "
         << (short)rgbpixel.GetRed()   << "  "
         << (short)rgbpixel.GetGreen() << "  "
         << (short)rgbpixel.GetBlue()  << endl;*/
    tpixel = (unsigned char)rgbpixel.GetRed();
#else
    const RGBPixelType& rgbpixel       = itRGB.Value();
    typename TImage::PixelType& tpixel = itTarget.Value();
    tpixel = rgbpixel;
#endif
  }

#else
  typename TImage::Pointer newImage = reader->GetOutput();
#endif

  double spacing[3];
  spacing[0] = spacing[1] = spacing[2] = voxelsize;
  newImage->SetSpacing( spacing );

  if( flip ) {

    //flip the image twice to rotate it back upwards
    itk::FixedArray<bool, voxfeDimension> flipAxes;

    flipAxes[0] = false;
    flipAxes[1] = true;
    flipAxes[2] = true;

    typedef itk::FlipImageFilter<TImage> FlipImageFilterType;
    typename FlipImageFilterType::Pointer flipFilter = FlipImageFilterType::New ();
    flipFilter->SetInput( newImage );
    flipFilter->SetFlipAxes(flipAxes);

    try {
      flipFilter->Update();
    }
    catch( itk::ExceptionObject & err )
    {
      logOut << "\n===============================================" << endl;
      logOut << "*** Exception caught for flip conversion ****" << endl;
      logOut << "(Were the correct file limits set?" << first << " -> " << last << ")" << endl;
      logOut << "ERROR: " << err << endl;
      logOut << "\n===============================================" << endl;
      logOut.close();

      newImage = 0;
      return newImage;
    }

    newImage = flipFilter->GetOutput();  //copy over
  }

  logOut.close();
  return newImage;
}


template<class TImage>
void
voxfeITKWriteToImageFileMHD( const TImage* inputImage, const char* filepath ) {

  typedef itk::ImageFileWriter< TImage >  WriterType;
	itk::SmartPointer< WriterType > writer = WriterType::New();

	writer->SetInput( inputImage );
	writer->SetFileName( filepath );
	//std::cout  << "Writing the image as " << std::endl << std::endl;
	//std::cout  << argv[2] << std::endl << std::endl;

	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject &ex)
	{
		std::cout << ex << std::endl;
	}

}


template<class InputImage, class OutputImagePixelType>
itk::SmartPointer< itk::Image<OutputImagePixelType, voxfeDimension> >
voxfeITKCreateEmptyImage( const InputImage* inputImage, OutputImagePixelType pxl, bool init )
{
  itk::ImageRegion<voxfeDimension> imgRegion = inputImage->GetLargestPossibleRegion();

	// have to use a smart pointer here or the next few lines (setting the region) will fail
  typedef itk::Image<OutputImagePixelType, voxfeDimension> OutputImageType;

	typename OutputImageType::Pointer newImage = OutputImageType::New();
	newImage->SetRegions( imgRegion );
	newImage->CopyInformation( inputImage );
	newImage->Allocate();
	if( newImage.IsNotNull() && init ) newImage->FillBuffer((OutputImagePixelType)0);

	return newImage;
}


template<class TImage, class CCLabelType>
set<short>*
voxfeITKProcessImage( TImage* inputImage, CCLabelType clabel,
                      const char& CH, const char* filename,
                      const char* algorithm,
                      const int _border[3] ) {

  if( !inputImage ) {
    cerr << "\n*** Image not defined (cannot perform CC) ***\n\n";
    return NULL;
  }

  //set up connected component tool
  typedef itk::Image<CCLabelType, voxfeDimension> CCImage;
  typedef itk::ImageRegionIteratorWithIndex<CCImage> CCImageIterator;

  typedef itk::ConnectedComponentImageFilter <voxfeGreyUCharImage, CCImage > CCImageFilterType;
  typename CCImageFilterType::Pointer connected = CCImageFilterType::New ();
  connected->SetInput(inputImage);
  try
  {
    connected->Update();
  }
  catch (itk::ExceptionObject &ex)
  {
	  std::cerr << ex << std::endl;
  }
  cout << "Number of objects: " << connected->GetObjectCount() << endl;
  if( connected->GetObjectCount() >= numeric_limits<CCLabelType>::max() ) {
    std::cerr << "\n*** Warning: Connected component image cannot handle this number of objects ***\n";
    std::cerr << "\n*** Warning: Some part of the object may have been lost.....                ***\n";
  }

  //relabel the components to get the largest
  typedef itk::RelabelComponentImageFilter<CCImage, CCImage> RelabelImageType;
  typename RelabelImageType::Pointer relabel = RelabelImageType::New();
  relabel->SetInput( connected->GetOutput() );

  relabel->InPlaceOn();

  try
  {
    relabel->Update();
  }
  catch (itk::ExceptionObject &ex)
  {
	  std::cerr << ex << std::endl;
  }
  cout << "Largest object: " << relabel->GetSizeOfObjectInPixels(1) << endl;

  //   DEBUG
  //voxfeITKWriteToImageFile( connected->GetOutput(), "/tmp/test-CC-output.mhd" );
  //voxfeITKWriteToImageFile( relabel->GetOutput(),   "/tmp/test-CC-relabel.mhd" );

  //we can free this now.....
  //connected = 0;  //fixme: would be better to do relabelling inplace ie set relabel->InPlaceOn()

  //output -----------------------------------------------------------
  set<short>* valid_labels = new set<short>;
  string voxfile( filename );

  if( CH == 'I' ) {

    //convert to short image and remove unwanted labels
    CCImageIterator itSource( relabel->GetOutput(),
                              relabel->GetOutput()->GetLargestPossibleRegion() );

    voxfeGreyUCharImageIterator itTarget( inputImage, inputImage->GetLargestPossibleRegion() );

    for( itSource.GoToBegin(),itTarget.GoToBegin(); !itSource.IsAtEnd(); ++itSource, ++itTarget )  {

	    assert( !itTarget.IsAtEnd() );

	    if( itSource.Get() != 1 ) {
        itTarget.Set( 0 );   //remove the pixel if it's not in the largest con compoent
      }
 	}

    voxfile.append( voxfeCCImageAppendString );
    voxfeITKWriteToImageFileMHD( inputImage, voxfile.c_str() );
  }
  else if( CH == 'V' ) {

    voxfile.append( voxfeCCModelAppendString ); //needed to put in script file

    string vxfPath, vxfFile;
    voxfeSplitFilename( voxfile, vxfPath, vxfFile );  //redone here to get the appended file

    //1. -------------- create a simple script file -----------------------
    string scriptfile( filename );
    scriptfile.append( voxfeCCScriptAppendString );
    //retrieve image & voxel sizes
    const voxfeRegion3D& region = inputImage->GetLargestPossibleRegion();
    const voxfeSize3D& size     = region.GetSize();
    const typename TImage::SpacingType& spacing = inputImage->GetSpacing();
    int _size[3]; double _spacing[3];
    for( int p=0; p<3; p++ ) {
      _size[p] = size[p] + (2 * _border[p]); //add in borders if any
      _spacing[p] = spacing[p];
    }

    rvsWriteVoxFEScriptStub( scriptfile.c_str(), vxfFile.c_str(), _spacing, _size, algorithm );

    //2. --------------- create the model file ----------------------------
    //modify input image to remove unconnected components
    CCImageIterator itCC( relabel->GetOutput(),
                          relabel->GetOutput()->GetLargestPossibleRegion() );

    voxfeGreyUCharImageIterator itTarget( inputImage, inputImage->GetLargestPossibleRegion() );

    unsigned long voxel_count = 0;
    for( itCC.GoToBegin(),itTarget.GoToBegin(); !itCC.IsAtEnd(); ++itCC, ++itTarget )  {

    	assert( !itTarget.IsAtEnd() );
	    const CCLabelType& pixelCC = itCC.Value();

  		//reject any pixels that where pixelCC != 1
	    if( pixelCC == 1 ) {
	      voxel_count++; //maintain a count for output of voxels
	      valid_labels->insert( itTarget.Get() );
	    }
    	else
        itTarget.Set( 0 );

	    //debug
	    //cout << "Pixel (relabelled): " << itCC.Value() << " set to: "  << (short)(itTarget.Get()) << endl;
    }

    cout << "Outputting to: " << voxfile.c_str() << "  Border: "
         << _border[0] << " "
         << _border[1] << " "
         << _border[2] << "\n";

    //output a model script file
    ofstream fout( voxfile.c_str() );
    fout << "\n" << voxel_count << "\n";

    //echo to screen
    std::cout << "\nVoxel count: " << voxel_count << "\n";

    unsigned long k = 1; //file counter
    for( itTarget.GoToBegin(); !itTarget.IsAtEnd(); ++itTarget )  {

      const unsigned char& pixel = itTarget.Value();
      if( pixel ) {
        voxfeIndex3D idx = itTarget.GetIndex();

        //fixme: not sure if the indices should start at 1??
        fout << k      << " " << (int)pixel  << " "
             << (idx[0] + _border[0] + voxfeVoxelOffset) << " "
             << (idx[1] + _border[1] + voxfeVoxelOffset) << " "
             << (idx[2] + _border[2] + voxfeVoxelOffset) << "\n";
        k++;
      }
    }

    fout.close();
  }
  
  return valid_labels;
}

#endif
#endif

