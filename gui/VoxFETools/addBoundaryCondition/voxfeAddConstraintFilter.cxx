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

#include "voxfeAddConstraintFilterDefines.h"
#include "voxfeAddConstraintFilter.h"

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPolyData.h"
#include "vtkCharArray.h"
#include "vtkIntArray.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkFieldData.h"
#include "vtkSmartPointer.h"
#include "vtkCellArray.h"
#include "vtkVertex.h"
#include "vtkUnstructuredGrid.h"
#include "vtkMultiBlockDataSet.h"

#include "vtkGlyphSource2D.h"
#include "vtkGlyph3D.h"
#include "voxfeGlyphAnnotationFilter.h"
#include "vtkAppendFilter.h"
#include "vtkAppendPolyData.h"

#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkCompositeDataSet.h"

//debug
int voxfeAddConstraintFilter::objectcounter = 0;

//vtkStandardNewMacro(voxfeAddConstraintFilter);
voxfeAddConstraintFilter* voxfeAddConstraintFilter::New()
{
  objectcounter++;
  return new voxfeAddConstraintFilter;
}
extern vtkObject* vtkInstantiatorvoxfeAddConstraintFilterNew();
vtkObject* vtkInstantiatorvoxfeAddConstraintFilterNew()
{
    return voxfeAddConstraintFilter::New();
}

//----------------------------------------------------------------------------
voxfeAddConstraintFilter::voxfeAddConstraintFilter()
{
  //ResetDefaultConstraint(); //fixme: works but doesn't update display widgets
  cout << "\n*******Number of input ports: " << this->GetNumberOfInputPorts() << "*******\n";

  AppendUGridData = vtkAppendFilter::New();
  AppendGlyphData = vtkAppendPolyData::New();
  IsDefined = false;
}

//----------------------------------------------------------------------------
voxfeAddConstraintFilter::~voxfeAddConstraintFilter()
{
  AppendUGridData->Delete();
  AppendGlyphData->Delete();
}

//----------------------------------------------------------------------------
void voxfeAddConstraintFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
int voxfeAddConstraintFilter::RequestData(vtkInformation *request,
                              vtkInformationVector **inputVector,
                              vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkDataObject *inputDataObj = NULL;
  if( inInfo ) {
    inputDataObj = vtkDataObject::GetData(inInfo);
    cout << "inputDataObj is a: " << inputDataObj->GetClassName() << "\n";
  }

#ifdef DEBUG
  //==== debug ======================================================================================
  int numInfoObj0 = inputVector[0]->GetNumberOfInformationObjects();
  cout << "\n===============================================================================";
  cout << "\n@[" << objectcounter << "] Num info obj on inputVector[0]: " << numInfoObj0 << "\n";

  cout << "\nNumber of pieces: " << inInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES())
       << "\nPiece number:     " << inInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER())
       << "\n";

  vtkCompositeDataSet* test_cd = vtkCompositeDataSet::GetData(inputVector[0], 0);
  vtkUnstructuredGrid* test_ug = vtkUnstructuredGrid::GetData(inputVector[0], 0);
  if( test_cd ) cout << "\nRetrieved Composite Data\n";
  if( test_ug ) cout << "\nRetrieved Unstruc Grid Data\n";
  //==== debug ======================================================================================
#endif

  // verify the input, selection and output
  if ( ! inputDataObj )
    {
    vtkErrorMacro(<<"No input specified");
    return 0;
    }
    
  int retVal = 1; //success
 
  //we assume surface points have been selected and data must be unstructured
  if( ! inputDataObj->IsA("vtkUnstructuredGrid") )
    {
    vtkErrorMacro(<<"No points selected");  //fixme: QMessage!?
    return 0;
    }

  //reject data pushes if previously called for this obj
  vtkUnstructuredGrid*  input  = 0;
//  if( this->IsDefined ) {
//    cout << "\n\n@[" << objectcounter << "] is already defined\n";
//    input  = AppendUGridData->GetOutput();
//  }
//  else {
    // get the input
    input  = vtkUnstructuredGrid::SafeDownCast(  inInfo->Get(vtkDataObject::DATA_OBJECT()) );
//  }

  vtkPoints*    pdata = input->GetPoints();
  vtkIdType numPoints = pdata->GetNumberOfPoints();

#ifdef VOXFE_ADDBC_DEBUG
  //debug=========================
  cout << "voxfeAddConstraintFilter::ReqestData -- Adding attributes to: "
       << numPoints << " points \n";
  //cout << "( cf. "  << n_pdata << " data items)\n";
  for( vtkIdType np=0; np<numPoints; np++ ){
    double* pt = pdata->GetPoint(np);
    cout << np << " : " << pt[0] << "  " << pt[1] << "  " << pt[2] << "\n";
  }
  //debug=========================
#endif

  if( input && (numPoints > 0) ) {
    cout << "AddBC request cast to: Unstructured Grid\n";
  }
  else {
    cout << "** AddBC Input is void or empty... **\n";
    return 0;
  }

  vtkFieldData* fdata = input->GetFieldData();

  //create field arrays and add them -- vtk replaces the named array if it already exists
  //1. Create an int array giving the drop-down parameters (defined in server xml)
  vtkSmartPointer<vtkIntArray> ia = vtkSmartPointer<vtkIntArray>::New();

  //2. Create another array, giving the axis constraints (0=free,1=fixed) or force vector (parallel) or force
  // endpoint (if to a point) depending upon BoundaryCondition
  vtkSmartPointer<vtkDoubleArray> da = vtkSmartPointer<vtkDoubleArray>::New();
  if( this->BoundaryCondition == VOXFE_BC_NONE ) {
    retVal = 0;  //don't set anything
  }
  else if( this->BoundaryCondition == VOXFE_BC_NODAL ) {

    for( int k=0; k<3; k++ ) {
      double dbl = this->AxisConstraint[k];
      da->InsertNextValue( dbl );  // ? use int array ?
      //cout << "Adding axis[" << k << "] = " << dbl << "\n";
    }
  }
  else if( this->BoundaryCondition == VOXFE_BC_FORCE_PARALLEL ) {

    for( int k=0; k<3; k++ ) { 
      da->InsertNextValue( this->ForceVector[k] );
      //cout << "Adding force_paral[" << k << "] = " << ForceVector[k] << "\n";
    }
  }
  else if( this->BoundaryCondition == VOXFE_BC_FORCE_TO_POINT ) {

    for( int k=0; k<3; k++ ) {
      da->InsertNextValue( this->ForceEndpoint[k] );
      //cout << "Adding force_endpt[" << k << "] = " << ForceEndpoint[k] << "\n";
    }
  }
  //  else if( this->BoundaryCondition == VOXFE_BC_NO_REMODEL ) {
  //
  //  }


  if( retVal ) {
    
    //update part 1.
    ia->InsertNextValue( this->BoundaryCondition );
    ia->InsertNextValue( this->ForceDistribution );

    //3. Add in user-supplied force magnitude
    da->InsertNextValue( this->ForceMagnitude );

    if( this->IsDefined ) {
      fdata->RemoveArray( VOXFE_CONSTRAINT_ATTRIBUTE );
      fdata->RemoveArray( VOXFE_VECTOR_ATTRIBUTE );
    }

    //4. Add the fields to the input data and delete the orig arrays
    ia->SetName( VOXFE_CONSTRAINT_ATTRIBUTE );
    fdata->AddArray(ia);
    da->SetName( VOXFE_VECTOR_ATTRIBUTE );
    fdata->AddArray(da);
  }

  // ===================================== output =====================================
  vtkInformation*      outInfo = outputVector->GetInformationObject(0);
  vtkMultiBlockDataSet* output = vtkMultiBlockDataSet::SafeDownCast( outInfo->Get(vtkDataObject::DATA_OBJECT()) );

  if( (this->BoundaryCondition >= VOXFE_BC_NODAL) && (this->BoundaryCondition <= VOXFE_BC_NO_REMODEL)  ) {

      //============================== ADD GLYPHS ==========================================
      //don't want to mark every point - too much, so find the bounds and mark a/the corners
      double bounds[6],center[3],normal[3];
      input->ComputeBounds();

      //build a small grid (after BuildUGrid.tcl)
      vtkSmartPointer<vtkUnstructuredGrid> ug = vtkSmartPointer<vtkUnstructuredGrid>::New();
      vtkSmartPointer<vtkPoints> glyphpoints = vtkSmartPointer<vtkPoints>::New();

      if( this->BoundaryCondition == VOXFE_BC_NODAL ) {
        //input->GetBounds(bounds);
        input->GetCenter(center);
        ug->Allocate( vtkIdType(1), vtkIdType(1) );
        glyphpoints->SetNumberOfPoints( 1 );
        //glyphpoints->InsertPoint(0,bounds[0],bounds[2],bounds[4]); // min x,y,z
        glyphpoints->InsertPoint(0,
                                 center[0] + (this->GetOffsetX() * this->GetOffsetStep()),
                                 center[1] + (this->GetOffsetY() * this->GetOffsetStep()),
                                 center[2] + (this->GetOffsetZ() * this->GetOffsetStep())
                                 );
      }
      else if( this->BoundaryCondition < VOXFE_BC_NO_REMODEL ) {
        input->GetCenter(center);
        ug->Allocate( vtkIdType(1), vtkIdType(1) );
        glyphpoints->SetNumberOfPoints( 1 );
        glyphpoints->InsertPoint(0,
                                 center[0] + (this->GetOffsetX() * this->GetOffsetStep()),
                                 center[1] + (this->GetOffsetY() * this->GetOffsetStep()),
                                 center[2] + (this->GetOffsetZ() * this->GetOffsetStep())
                                 );

        if( this->BoundaryCondition == VOXFE_BC_FORCE_PARALLEL ) {

          normal[0] = this->ForceVector[0];
          normal[1] = this->ForceVector[1];
          normal[2] = this->ForceVector[2];
        }
        else
        if( this->BoundaryCondition == VOXFE_BC_FORCE_TO_POINT ) {
          normal[0] = this->ForceEndpoint[0] - center[0];
          normal[1] = this->ForceEndpoint[1] - center[1];
          normal[2] = this->ForceEndpoint[2] - center[2];
        }
      }

      ug->SetPoints(glyphpoints);

      vtkSmartPointer<vtkVertex> vertexA = vtkSmartPointer<vtkVertex>::New(); // the 'cells'
      //vtkSmartPointer<vtkVertex> vertexB = vtkSmartPointer<vtkVertex>::New();
      vertexA->GetPointIds()->SetId(0, 0); // cell 0 relates to point 0;
      //vertexB->GetPointIds()->SetId(1, 1); // cell 1 relates to point 1;

      ug->InsertNextCell( vertexA->GetCellType(), vertexA->GetPointIds() );
      //ug->InsertNextCell( vertexB->GetCellType(), vertexB->GetPointIds() );

      vtkSmartPointer<vtkDoubleArray> vnormal = vtkSmartPointer<vtkDoubleArray>::New();
      vnormal->SetNumberOfComponents(3);
      vnormal->SetNumberOfTuples(1);
      vnormal->SetTuple(0,normal);
      ug->GetPointData()->SetNormals(vnormal);


      //need to transfer the ancillary data too, so that our glyph class knows how to act
      vtkFieldData*   ug_fdata = ug->GetFieldData();
      ug_fdata->AddArray(ia);
      ug_fdata->AddArray(da);


      //merge adding array data with glyph annotation
      vtkSmartPointer<voxfeGlyphAnnotationFilter> glyph = vtkSmartPointer<voxfeGlyphAnnotationFilter>::New();
      glyph->SetGlyphScale( this->GetGlyphScale() );
      glyph->SetInputData( ug );
      glyph->Update();

      //add to the output
      if( this->IsDefined ) {

        //if the data has already input, we might still want to update the glyph
        AppendGlyphData->RemoveAllInputs();
        AppendGlyphData->AddInputData( glyph->GetOutput() );
        AppendGlyphData->Update();
      }
      else {
        AppendUGridData->AddInputData( input );
        AppendUGridData->Update();
        AppendGlyphData->AddInputData( glyph->GetOutput() );
        AppendGlyphData->Update();
      }

      cout << "Added 1st block to this condition\n";
      //set the MB blocks
      unsigned int k = 0;
      output->Initialize();
      output->SetNumberOfBlocks(2);

      output->SetBlock(k,AppendUGridData->GetOutput());
      output->GetMetaData(k)->Set(vtkCompositeDataSet::NAME(), "BoundaryCondition");

      cout << "Added 2nd block to this condition\n";
      k = 1;
      output->SetBlock(k,AppendGlyphData->GetOutput() );
      output->GetMetaData(k)->Set(vtkCompositeDataSet::NAME(), "GlyphAnnotation");

      this->IsDefined = true;
      cout << "Added condition.... \n";
  }
  
  
  return retVal;
}


//----------------------------------------------------------------------------
int voxfeAddConstraintFilter::FillInputPortInformation( int port, vtkInformation* info )
{
  /*if (!this->Superclass::FillInputPortInformation(port, info))
    {
    return 0;
    }*/

  //The source dataset must be unstructured points
  if ( port == 0 )
    {
    info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid" );
    return 1;
    }
  return 0;
}


