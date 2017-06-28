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

#include "voxfeITKReader.h"
#include "voxfeConnectedComponentImageFilter.hxx"

#include "ReadVoxFEScript.hxx"
#include "convertVoxFEScript.hxx"

// ================== additional VTK functions ===================================================

#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkDataObject.h>
#include <vtkSmartPointer.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkGeometryFilter.h>
#include <vtkPointData.h>

#include <vtkVoxel.h>
#include <vtkCellData.h>
#include <vtkIntArray.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPolyDataWriter.h>
#include <vtkThreshold.h>
#include <vtkCharArray.h>
#include <vtkDoubleArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkSphereSource.h>

// ================== additional Reader functions =================================================

  /** Create a vtk unstructured grid from voxfe data (file).
     \brief From data read in using ReadVoxelScript functions create a single vtk file. Materials
            are described as voxel cell data.
     \param nodes Stores an in-memory mapping from unique node number to the
            sequence 1...n, for n nodes.
     \param element Element meta data.
     \param n_elements The total number of elements.
     \return A pointer to the created vtkUnstructuredGrid.
   */
  template<class T>
  vtkUnstructuredGrid* CreateVTKMeshFromVoxFEModel(
                         const map<T, T>& nodes,
                         const Element<T>* element,
                         const T& n_elements,
                         const ModelMetaData<T>& model );

  /** Create surface vtk polydata meshes according to groups given. A geometry filter
      is also applied so that only the surface geometry is retained.
      \brief The group names (from ImageGroupData) are used to name the output files.
      \param groupData Describes the groups within the voxel data.
      \param unstructuredGrid The vtk data to split.
      \param v_polydata gathers the filtered polydata into a std vector.
   */
  void CreateMultiBlockDataGroups( const ImageGroupData* groupData,
                            vtkUnstructuredGrid* unstructuredGrid,
                            vector< vtkDataObject* >& v_polydata,
                            bool removeInternalGeometry=true);
  
  /** Create a listing common materials from defaults and user-defined input */
  void AssembleRefMaterials(voxfeITKReader* reader, ImageGroupData::RefMaterialMap& refMap);

// ================== additional VTK functions ===================================================


vtkStandardNewMacro(voxfeITKReader);

voxfeITKReader::voxfeITKReader() :
    FileName(NULL),
    DisplacementDataFile(NULL),
    FileSpecifier(NULL)
{
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
  SetFileSpecifier( "default%03d.bmp" );
}

int voxfeITKReader::FillOutputPortInformation(int vtkNotUsed (port),
                                             vtkInformation *info )
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
  return 1;
}


int voxfeITKReader::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  // get the info object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the ouptut object
  vtkMultiBlockDataSet *output =
      vtkMultiBlockDataSet::SafeDownCast (
        outInfo->Get ( vtkMultiBlockDataSet::DATA_OBJECT() ) );            

  // Change the working directory so we can load voxfe image data (which assumes the same
  // directory for the raw data)
  string vxfPath, vxfFile;
  voxfeSplitFilename( string(this->FileName), vxfPath, vxfFile );
  if( chdir( vxfPath.c_str() ) )
    cerr << "\n\n *** Failed to change working directory\n" << endl << endl;
  else
    cout << "\n\n*** Changing working directory to: " << vxfPath << endl << endl;

  vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = 0;
  voxfe_vtk_file_created fileCreated = VOXFE_VTK_FILE_NONE;

  //read the image group file, unless the user gave the script file directly (check first!)
  ImageGroupData* groupDataPtr = 0;
  string scriptFName, scriptFExt;
  voxfeSplitFilename( vxfFile, scriptFName, scriptFExt, "." );
  if( scriptFExt == "script" ) {

    //in this case treat as VOXFE_SCRIPT_ONLY with no groups
    groupDataPtr = new ImageGroupData;
    groupDataPtr->ImageFile = vxfFile; //set this as same to select VOXFE_SCRIPT_ONLY
  } 
  else
    groupDataPtr = rvsReadImageGroupFile( this->FileName );
    
  std::auto_ptr<ImageGroupData> groupData;
  std::auto_ptr< set<short> >   labelsPresent;
  vector< vtkDataObject* > v_polydata;

  if( groupDataPtr ) {

    groupData.reset( groupDataPtr );

    //determine the input file type    
    string imageFName, imageFExt;
    voxfeSplitFilename( groupDataPtr->ImageFile, imageFName, imageFExt, "." );
    cout << "Found file extension: " << imageFExt << "\n";

    if     ( imageFExt == "mhd" ) 
      SetInputImageType( VOXFE_ITK_MHD );
    else if( (imageFExt == "bmp") || (imageFExt == "BMP") )
      SetInputImageType( VOXFE_BMP_SERIES );
	else if( (imageFExt == "png") || (imageFExt == "PNG") )
      SetInputImageType( VOXFE_PNG_SERIES );
	else if( (imageFExt == "tif") || (imageFExt == "TIF") || (imageFExt == "tiff") || (imageFExt == "TIFF") )
      SetInputImageType( VOXFE_TIFF_SERIES );
    else if( imageFExt == "script" )
      SetInputImageType( VOXFE_SCRIPT_ONLY );

#ifdef VOXFE_USE_ITK
    //read the image filethis->GetBorderOffset()
    voxfeGreyUCharImagePointer image = 0;  //should be allocated below
    if( this->InputImageType == VOXFE_ITK_MHD ) {

      image = voxfeITKReadImageFileMHD( groupDataPtr->ImageFile.c_str(), image.GetPointer(),
                                        groupDataPtr->VoxelSize, this->ImageFlip );
    }
    else if( (this->InputImageType == VOXFE_BMP_SERIES) || 
		     (this->InputImageType == VOXFE_PNG_SERIES) || 
			 (this->InputImageType == VOXFE_TIFF_SERIES) ) {

      string subpath, _spec; //user may have given a valid bmp filespec with "./", so remove
      voxfeSplitFilename(  groupDataPtr->ImageFile, subpath, _spec );

      // Fixme: Assume the supplied spec (in the voxfe file) is correct for now
      SetFileSpecifier(_spec.c_str());

      if( groupDataPtr->VoxelSize > 0.0 ) {

        string pathWithFileSpec(vxfPath + "/" + this->FileSpecifier);
        image = voxfeITKReadImageFilePNG( pathWithFileSpec.c_str(),
                                          image.GetPointer(),
                                          groupDataPtr->VoxelSize,
                                          this->FileRange[0],
                                          this->FileRange[1],
                                          this->FileRange[2],
                                          this->ImageFlip
                                         );
      }
      else {
        cerr << "\n\n *** Voxel size has not been set in the .voxfe file -- please check specification.\n" << endl << endl;
      }
    }
#endif

    string scriptFile;
    if( this->InputImageType == VOXFE_SCRIPT_ONLY ) {

      //read a script file directly
      string subpath, _spec; //user may have given a filename with "./", so remove
      voxfeSplitFilename(  groupDataPtr->ImageFile, subpath, _spec );

      scriptFile = _spec;
    }
#ifdef VOXFE_USE_ITK    
    else
    if( image ) {      //==== ITK stuff may throw an exception which will result in a null image ====

      char CH = 'V';  //generate voxfe model file  -- OR --   use: char CH = 'I';  //to generate image

      //do connected components (64 bit version) and save to model file
      unsigned long  ccLabel = 0;
      const char* algorithm = (this->GetSolverAlgorithm() == KSPCG_PCJACOBI) ? "KSPCG PCJACOBI" : "JAPCG";
      labelsPresent.reset( voxfeITKProcessImage( image.GetPointer(),
                                                 ccLabel, CH,
                                                 this->FileName,
                                                 algorithm,
                                                 groupDataPtr->BorderOffset
                                                 ) );
    
      if( labelsPresent.get() ) {
        scriptFile = this->FileName;
        scriptFile.append( voxfeCCScriptAppendString );

        //debug
        cout << "Found labels:\n";
        for( set<short>::iterator it=labelsPresent->begin(); it!= labelsPresent->end(); ++it) {
          cout << *it << "  ";
        }
        cout << endl;
      }
        
      //at this point image can be released
      image = 0;      
    }
#endif    
    
    //If scriptFile is set then we can continue, else an error will be returned
    //read the model file to create a .vtk
    if( scriptFile.length() > 0 ) {

      //borders -- we now set this in group file, *only* applies to image data, set this to zero for script/model
      int border[3];
      border[0] = border[1] = border[2] = 0;

      cout << "Reading script file: " << scriptFile << endl;
      ModelMetaData<rvsIntegerType> model;
      if( rvsReadModelData( scriptFile.c_str(), border, model ) ) {

        //create and fill voxel table
        Element<rvsIntegerType>* element = 0;
        rvsIntegerType  minmaxZ[2], n_elements = 0;
        map<rvsIntegerType, rvsIntegerType> nodes;   //mapping from unique node index to position in map 
                                                     // (latter used in output file)
        if( rvsReadVoxelFile( model.voxelFile.c_str(), element, n_elements, model,  minmaxZ, nodes, border ) ) {
        
          //sort the groups into a coherent bunch
          if( this->InputImageType == VOXFE_SCRIPT_ONLY ) { //we can only now get the labels
            
            set<short>* label_set = new set<short>;
            rvsGetMaterialSet( element, n_elements, *label_set );
            labelsPresent.reset( label_set );
            
            //we also need to retrieve the voxel size to set as field data
            groupDataPtr->VoxelSize = model.voxSize[0];
          }
          groupData->TallyGroups( *labelsPresent );  //add in all groups (even singles not previously listed)

          //store the global node order for reading displacements/strains in later
          string nodeMapFile( scriptFile.c_str() ); //this->FileName);
          nodeMapFile.append(VOXFE_NODEMAP_FILE_EXTN);
          rvsWriteNodeMapKeys( nodeMapFile.c_str(), nodes );  //fixme: prob better always to create -- model might have changed
                                                              // if slices being read.. but ... make it so user can select??

          if( this->GetCreateMetisFile() ) {

            string metisfile = string(scriptFile) + string(".metis");
            if( !WriteMetisMeshFile( metisfile.c_str(), nodes, element, n_elements, model ) ) {

              cerr << "@@@ Could not write out the Metis file: " <<  metisfile << " @@@" << endl;
            }
            else cout << "Written out the Metis file: " << metisfile << endl;

          }

          if( this->GetCreateFeltFile() ) {

            string feltfile = string(scriptFile) + string(".flt");
            if( !WriteFeltMeshFile( feltfile.c_str(), nodes, element, n_elements, model ) ) {

              cerr << "@@@ Could not write out the Felt file: " <<  feltfile << " @@@" << endl;
            }
            else cout << "Written out the Metis file: " << feltfile << endl;

          }


          //======================= debug =========================
          //nodeMapFile.append(".txt");
          //rvsWriteNodeMapKeysASCII( nodeMapFile.c_str(), nodes );
          //======================= debug =========================


          //output a materials file if it does not exist
#ifdef VOXFE_USE_ITK
          if( ! itksys::SystemTools::FileExists( VOXFE_MATERIALS_FILE ) ) {
#else
          //lots of possible alternatives here using boost, stat etc. I just do something simple
          fstream foo;
          foo.open( VOXFE_MATERIALS_FILE );
          if( foo.is_open()) {
            foo.close();
          }
          else {
#endif
            //assemble default values for a handful of pre-defined materials
            //fixme: the strings could be user-defined for greater flexibility
            ImageGroupData::RefMaterialMap refMap;
	        AssembleRefMaterials(this, refMap);
            groupData->WriteMaterialsFile( VOXFE_MATERIALS_FILE, refMap );
          }

          //vtkFile  = scriptFile + string(".vtk");

          //gather the voxel data into groups suggested by voxfe file....
          unstructuredGrid = CreateVTKMeshFromVoxFEModel( nodes, element, n_elements, model );

          //... and generate grouped surface polydata
          CreateMultiBlockDataGroups( groupData.get(), unstructuredGrid, v_polydata, this->GetGenerateSurfaceMesh() );

          fileCreated = VOXFE_VTK_FILE_INMEMORY;

          //clean up
          if( element ) delete[] element;
        }
        else {
          cerr << "Could not read voxel file...." << endl;
          if( element ) delete[] element;
        }
      }
      else cerr << "Could not read the input script file: " << scriptFile << endl;

    } //if( scriptfile.length()  )

  } //if( groupDataPtr )
  else {
    // we will still have fileCreated == VOXFE_VTK_FILE_NONE, but issue warning
    vtkWarningMacro( "Cannot access image group file -- rendering default" );
  }



  if( fileCreated != VOXFE_VTK_FILE_INMEMORY ) {   //fixme: fileCreated could now be a bool? -- meshCreated??
  
    //---- output this trinket as a failsafe ------------------------------------------
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->InsertNextPoint(0.0, 0.0, 0.0);
    polydata->SetPoints(points);

    vtkSmartPointer<vtkVertexGlyphFilter> glyphFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    glyphFilter->SetInputData(polydata);
    glyphFilter->Update();

    vtkPolyData* poly = vtkPolyData::New();
    poly->DeepCopy( glyphFilter->GetOutput() );
    v_polydata.push_back( poly );
  }


  //create the multiblock
  if( v_polydata.size() > 0 ) {

    cout << "Assembling multiblock with: " << v_polydata.size() << " objects" << endl;
    for( unsigned int k=0; k<v_polydata.size(); ++k) {

      output->SetBlock ( k, v_polydata[k] );
      
      //identify the data block with a name
      if( fileCreated == VOXFE_VTK_FILE_NONE ) 
        output->GetMetaData(k)->Set(vtkCompositeDataSet::NAME(), "NO_VOXEL_DATA");
      else {
        if( this->GetGenerateSurfaceMesh() )
          output->GetMetaData(k)->Set(vtkCompositeDataSet::NAME(), groupDataPtr->Group[k].MaterialName.c_str());
        else
          output->GetMetaData(k)->Set(vtkCompositeDataSet::NAME(), "All_data");
      }

      
      //mark for deletion to be safe
      v_polydata[k]->Delete();
    }

    //1. Adds some field data to the geometry (here the original file name)
    vtkFieldData *fd = output->GetFieldData();
    vtkCharArray *ca = vtkCharArray::New();

    string scriptFile;
    if( this->InputImageType == VOXFE_SCRIPT_ONLY )
      scriptFile = groupDataPtr->ImageFile;  //use the existing script file name
    else {
      scriptFile = this->FileName;
      scriptFile.append( voxfeCCScriptAppendString );  //recreate the script file name
    }


    //We could set this as the name, but better to access it by "voxfefile" instead
    for( int p=0; p<scriptFile.length(); p++ ) {

      ca->InsertNextValue(scriptFile[p]);
    }
    ca->InsertNextValue('\0');  //insert null to signify end
    ca->SetName(VOXFE_VOXFEFILE);
    fd->AddArray(ca);
    ca->Delete();

    //2. Add the input voxel size, so we can retrieve the node indices
    vtkDoubleArray *da = vtkDoubleArray::New();
    da->SetName(VOXFE_VOXELSIZE);
    da->InsertNextValue( groupDataPtr->VoxelSize );
    fd->AddArray(da);
    da->Delete();
  }

  return 1;
}


void voxfeITKReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "File Name: "
      << (this->FileName ? this->FileName : "(none)") << "\n";
}


// =============================================================================================================
//                              Additional VTK helper functions
// =============================================================================================================

template<class T>
vtkUnstructuredGrid* CreateVTKMeshFromVoxFEModel(
                       const map<T, T>& nodes,
                       const Element<T>* element,
                       const T& n_elements,
                       const ModelMetaData<T>& model ) {

  //create list of nodes
  vtkSmartPointer<vtkPoints> voxelPoints = vtkSmartPointer<vtkPoints>::New();
  voxelPoints->SetNumberOfPoints( (vtkIdType)nodes.size() );

  T i[3], rem;
  typename map<T,T>::const_iterator mit = nodes.begin();
  for( vtkIdType k=0 ;  mit != nodes.end(); mit++, k++ ) {

    //get node indices from unique node index
    i[2] = (mit->first) / model.XY; // "z"
    rem  = (mit->first) % model.XY;
    i[1] = rem / model.X;     // "y"
    i[0] = rem % model.X;     // "x"

    //fixme: vtkIdType is currently int (?)
    voxelPoints->InsertPoint( k,
                              i[0] * model.voxSize[0],
                              i[1] * model.voxSize[1],
                              i[2] * model.voxSize[2]
                            );
  }

  //======================================================================================

  vtkUnstructuredGrid* unstructuredGrid = vtkUnstructuredGrid::New();
  unstructuredGrid->Allocate( (vtkIdType)n_elements );
  unstructuredGrid->SetPoints( voxelPoints );

  //output elements -- Nodes: VTK same order as VOXFE, no need to reorder
  T a[NODESPERELEMENT];

  for( T k = 0; k < n_elements; k++ ) {

    //get the order of nodes for this element
    rvsGetVoxelNodes( model, element[k].index, a );

    //create a vtk voxel
    vtkSmartPointer<vtkVoxel> aVoxel = vtkSmartPointer<vtkVoxel>::New();

    for( T m = 0; m < NODESPERELEMENT; m++ ) {

      typename map<T,T>::const_iterator mit = nodes.find( a[m] );
      if( mit != nodes.end() )
        aVoxel->GetPointIds()->SetId( m, vtkIdType(mit->second - 1) ) ; //0-based indices in VTK!
    }

    unstructuredGrid->InsertNextCell( aVoxel->GetCellType(), aVoxel->GetPointIds() );
  }

  //mimic different material types with attached cell attribute data
  vtkSmartPointer<vtkIntArray> intArray = vtkSmartPointer<vtkIntArray>::New();
  intArray->SetNumberOfTuples( n_elements );
  for( T k = 0; k < n_elements; k++ ) {

    float material = element[ k ].elType;
    intArray->SetTuple( vtkIdType(k), &material );
  }

  intArray->SetName( "MaterialType" );
  unstructuredGrid->GetCellData()->SetScalars( intArray );

  return unstructuredGrid;
}



void CreateMultiBlockDataGroups( const ImageGroupData* groupData,
                          vtkUnstructuredGrid* unstructuredGrid,
                          vector< vtkDataObject* >& v_polydata,
                          bool removeInternalGeometry ) {

  //to do: add pvtk output to gather the vtk files
  //
#ifdef ORIG_MB_GROUP
  //loop through groups to output surface geometries
  for( short p=0; p<groupData->Group.size(); p++ ) {

    vtkSmartPointer<vtkThreshold> thresh = vtkSmartPointer<vtkThreshold>::New();
    thresh->SetInputData(unstructuredGrid);

    thresh->SetAllScalars(0);
    thresh->ThresholdBetween( groupData->Group[p].Start, groupData->Group[p].Stop );
    thresh->SetInputArrayToProcess( 1, 0, 0, 0, "MaterialType" );
    thresh->UseContinuousCellRangeOn();
    thresh->Update();

    if( thresh->GetOutput()->GetNumberOfCells() == 0 ) {

      cerr << groupData->Group[p].MaterialName << " has no cells here\n\n";
      continue;
    }
    else
      cout << "\nSetting start: " << groupData->Group[p].Start << " stop: "
           << groupData->Group[p].Stop << " for " << groupData->Group[p].MaterialName
           << "\n************************************************************\n";

    if( removeInternalGeometry ) {

      vtkSmartPointer<vtkGeometryFilter> geometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
      geometryFilter->SetInputData(thresh->GetOutput());

      cout << "\nRemoving geometry for group: " << groupData->Group[p].MaterialName << "....\n";
      geometryFilter->Update();
      cout << " done!\n\n";

      //set up output files
      //vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
      vtkPolyData* polydata = vtkPolyData::New();  //deleted when multi-block assembled
      polydata->DeepCopy( geometryFilter->GetOutput() ); //need to copy to retain
      v_polydata.push_back( polydata );
    }
    else {

      cout << "\nAdding group: " << groupData->Group[p].MaterialName << " without surface extraction....\n";

      //store the thresholded unstructuredGrid data
      //vtkSmartPointer<vtkUnstructuredGrid> ugriddata = vtkSmartPointer<vtkUnstructuredGrid>::New();
      vtkUnstructuredGrid* ugriddata = vtkUnstructuredGrid::New();  //deleted when multi-block assembled
      ugriddata->DeepCopy( thresh->GetOutput() ); //need to copy to retain
      v_polydata.push_back( ugriddata );
    }


    

#ifdef DUMP_TO_VTK_FILES

    vtkSmartPointer<vtkPolyDataWriter> pwriter = vtkSmartPointer<vtkPolyDataWriter>::New();
    pwriter->SetInputData(polydata);

    //dump to file
    string polyfile(groupData->Group[p].MaterialName);
    polyfile.append("_group.vtk");
    pwriter->SetFileName( polyfile.c_str() );

    //pwriter->SetFileTypeToASCII();
    pwriter->SetFileTypeToBinary();

    cout << "\nWriting...." << polyfile.c_str();
    pwriter->Update();
    cout << " ..... written!\n\n";
#endif
    
  }

#else

  /* In this version, we only create multi-blocks *if* we are also removing internal geometry
   *
   * Otherwise we pass everything as one lump, to be the first block......
   */
  if( removeInternalGeometry ) {

    for( short p=0; p<groupData->Group.size(); p++ ) {

      vtkSmartPointer<vtkThreshold> thresh = vtkSmartPointer<vtkThreshold>::New();
      thresh->SetInputData(unstructuredGrid);

      thresh->SetAllScalars(0);
      thresh->ThresholdBetween( groupData->Group[p].Start, groupData->Group[p].Stop );
      thresh->SetInputArrayToProcess( 1, 0, 0, 0, "MaterialType" );
      thresh->UseContinuousCellRangeOn();
      thresh->Update();

      if( thresh->GetOutput()->GetNumberOfCells() == 0 ) {

        cerr << groupData->Group[p].MaterialName << " has no cells here\n\n";
        continue;
      }
      else
        cout << "\nSetting start: " << groupData->Group[p].Start << " stop: "
             << groupData->Group[p].Stop << " for " << groupData->Group[p].MaterialName
             << "\n************************************************************\n";

      vtkSmartPointer<vtkGeometryFilter> geometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
      geometryFilter->SetInputData(thresh->GetOutput());

      cout << "\nRemoving geometry for group: " << groupData->Group[p].MaterialName << "....\n";
      geometryFilter->Update();
      cout << " done!\n\n";

      //set up output files
      vtkPolyData* polydata = vtkPolyData::New();  //deleted when multi-block assembled
      polydata->DeepCopy( geometryFilter->GetOutput() ); //need to copy to retain
      v_polydata.push_back( polydata );
    }
  }
  else {

    cout << "\nAdding without grouping to first block (without surface extraction....)\n";

    //store the thresholded unstructuredGrid data
    //vtkUnstructuredGrid* ugriddata = vtkUnstructuredGrid::New();  //deleted when multi-block assembled
    //ugriddata->DeepCopy( thresh->GetOutput() ); //need to copy to retain
    v_polydata.push_back( unstructuredGrid );
  }


#endif
}

void AssembleRefMaterials(voxfeITKReader* reader, ImageGroupData::RefMaterialMap& refMap) {
  
  double data[2];
  reader->GetDefaultBoneMaterial(data);
  refMap.insert( std::pair<const char*, MaterialType>("bone",  MaterialType(data[0],data[1]) ) );
  
  reader->GetDefaultToothMaterial(data);
  refMap.insert(std::pair<const char*, MaterialType>("teeth",  MaterialType(data[0],data[1]) ) );
  
  reader->GetDefaultMarrowMaterial(data);
  refMap.insert(std::pair<const char*, MaterialType>("marrow", MaterialType(data[0],data[1]) ) );
  
  reader->GetDefaultChitinMaterial(data);
  refMap.insert(std::pair<const char*, MaterialType>("chitin", MaterialType(data[0],data[1]) ) );
  
  reader->GetDefaultMetalMaterial(data);
  refMap.insert(std::pair<const char*, MaterialType>("metal",  MaterialType(data[0],data[1]) ) );
  
}



/*
 *  //write an unstructured grid
 *

    //no separation of materials, dump out everything
    vtkSmartPointer<vtkUnstructuredGridWriter> pwriter = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    pwriter->SetInputData(unstructuredGrid);

    //dump to file
    string polyfile(groupData->Group[p].MaterialName);
    polyfile.append("_group.vtk");
    pwriter->SetFileName( polyfile.c_str() );

    //pwriter->SetFileTypeToASCII();
    pwriter->SetFileTypeToBinary();
    pwriter->Update();


 *
 */

