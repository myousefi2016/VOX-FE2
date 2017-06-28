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

#include "vtkObjectFactory.h" //for new() macro
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkFieldData.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkCharArray.h"
#include "vtkIntArray.h"
#include "vtkDoubleArray.h"
#include "vtkUnstructuredGrid.h"

#include "vtkGlyphSource2D.h"
#include "vtkGlyph3D.h"
#include "vtkConeSource.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkAppendPolyData.h"

#include <iostream>
#include <string>
#include <fstream>
#include <limits>
using std::cout;
using std::cerr;
using std::string;
using std::ofstream;
using std::endl;
using std::fixed;
using std::numeric_limits;

#include "../VoxFEDefines.h"
#include "voxfeGlyphAnnotationFilter.h"





vtkStandardNewMacro(voxfeGlyphAnnotationFilter);

//-----------------------------------------------------------------------------
voxfeGlyphAnnotationFilter::voxfeGlyphAnnotationFilter()
{
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}


//-----------------------------------------------------------------------------
voxfeGlyphAnnotationFilter::~voxfeGlyphAnnotationFilter()
{

}


//----------------------------------------------------------------------------
int voxfeGlyphAnnotationFilter::FillInputPortInformation( int port, vtkInformation* info )
{

  switch (port)
  {
    //point data
    case 0:
      info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid" );
      return 1;

  default:
      vtkWarningMacro("Invalid input port " << port << " requested.");
      break;
  }
  
  return 0;
}


//----------------------------------------------------------------------------
int voxfeGlyphAnnotationFilter::RequestData(vtkInformation *vtkNotUsed(request),
                              vtkInformationVector **inputVector,
                              vtkInformationVector *outputVector)
{
  //[1]. get the first info object and model
  vtkInformation *inputModelInfo  = inputVector[0]->GetInformationObject(0);
  vtkUnstructuredGrid* inputBC    = vtkUnstructuredGrid::SafeDownCast( inputModelInfo->Get(vtkDataObject::DATA_OBJECT()) );

  if( !inputBC ) return 0; // no deal !
        
  vtkPoints* pointSet = inputBC->GetPoints();
  vtkIdType numPts    = pointSet->GetNumberOfPoints();
  
  cout << "voxfeGlyphAnnotationFilter found: " << numPts << " points\n";

  vtkFieldData* fd   = inputBC->GetFieldData();
  vtkIntArray*  ia   = static_cast<vtkIntArray*>(fd->GetArray( VOXFE_CONSTRAINT_ATTRIBUTE ));
  vtkDoubleArray* da = static_cast<vtkDoubleArray*>(fd->GetArray( VOXFE_VECTOR_ATTRIBUTE ));
  
  int conditionType = -1, forceDistrib  = -1;
  if( ia ) {
    conditionType = *(ia->GetPointer(0));
    forceDistrib  = *(ia->GetPointer(1));
  }

  double user_mag = -1.0, vec[3] = { -1.0,-1.0,-1.0 };
  if( da ) { 
    vec[0]   = *(da->GetPointer(0)); 
    vec[1]   = *(da->GetPointer(1)); 
    vec[2]   = *(da->GetPointer(2)); 
    user_mag = *(da->GetPointer(3));
  }

  //[2]. Define the glyphs we might attach (note cones found to be too obscuring)
  vtkSmartPointer<vtkGlyphSource2D> coneX = vtkSmartPointer<vtkGlyphSource2D>::New();
    coneX->SetGlyphTypeToTriangle();
    coneX->SetScale( this->GetGlyphScale() );
    coneX->FilledOff();
    coneX->Update();
  vtkSmartPointer<vtkGlyphSource2D> coneY = vtkSmartPointer<vtkGlyphSource2D>::New();
    coneY->SetGlyphTypeToTriangle();
    coneY->SetScale( this->GetGlyphScale() );
    coneY->FilledOff();
    coneY->Update();
  vtkSmartPointer<vtkGlyphSource2D> coneZ = vtkSmartPointer<vtkGlyphSource2D>::New();
    coneZ->SetGlyphTypeToTriangle();
    coneZ->SetScale( this->GetGlyphScale() );
    coneZ->FilledOff();
    coneZ->Update();

  //position the 'X' cone
  vtkSmartPointer<vtkTransform> transformX = vtkSmartPointer<vtkTransform>::New();
    transformX->Translate(-0.5 * this->GetGlyphScale(), 0.0, 0.0);
    transformX->RotateZ( 270.0 );
  vtkSmartPointer<vtkTransformPolyDataFilter> txPolyDataX = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    txPolyDataX->SetInputData( coneX->GetOutput() );
    txPolyDataX->SetTransform( transformX );
    
  //position the 'Y' cone
  vtkSmartPointer<vtkTransform> transformY = vtkSmartPointer<vtkTransform>::New();
    transformY->RotateY( 90.0 );
    transformY->Translate(0.0, -0.5 * this->GetGlyphScale(), 0.0);
  vtkSmartPointer<vtkTransformPolyDataFilter> txPolyDataY = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    txPolyDataY->SetInputData( coneY->GetOutput() );
    txPolyDataY->SetTransform( transformY );
    
  //position the 'Z' cone
  vtkSmartPointer<vtkTransform> transformZ = vtkSmartPointer<vtkTransform>::New();
    transformZ->Translate(0.0, 0.0, 0.5 * this->GetGlyphScale());
    transformZ->RotateX( 270.0 );
  vtkSmartPointer<vtkTransformPolyDataFilter> txPolyDataZ = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    txPolyDataZ->SetInputData( coneZ->GetOutput() );
    txPolyDataZ->SetTransform( transformZ );       
    
  vtkSmartPointer<vtkGlyphSource2D> arrow = vtkSmartPointer<vtkGlyphSource2D>::New();
    arrow->SetGlyphTypeToArrow();
    arrow->SetScale( 2.0 * this->GetGlyphScale() );  //doubled as arrow appears small to triangle
    arrow->FilledOff();
    arrow->Update();
  vtkSmartPointer<vtkTransform> transformArrow = vtkSmartPointer<vtkTransform>::New();
    transformArrow->Translate(0.5 * this->GetGlyphScale(), 0.0, 0.0);
  vtkSmartPointer<vtkTransformPolyDataFilter> txArrow = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    txArrow->SetInputData( arrow->GetOutput() );
    txArrow->SetTransform( transformArrow );

  // Glyph the gradient vector (with arrows)
  vtkSmartPointer<vtkGlyph3D> glyph = vtkSmartPointer<vtkGlyph3D>::New();
    glyph->SetInputData(inputBC);
   
  if( ia && da ) { 
    
    if( conditionType == VOXFE_BC_NODAL ) {
    
      cout << "voxfeGlyphAnnotationFilter found nodal BC \n"; 

      //add glyphs as per each restricted axis
      vtkSmartPointer<vtkAppendPolyData> append = vtkSmartPointer<vtkAppendPolyData>::New();
      if( vec[0] > 0.0 ) {
        cout << "Found X condition\n";
        append->AddInputConnection(txPolyDataX->GetOutputPort());
      }
      if( vec[1] > 0.0 ) {
        cout << "Found Y condition\n";
        append->AddInputConnection(txPolyDataY->GetOutputPort());
      }
      if( vec[2] > 0.0 ) {
        cout << "Found Z condition\n";
        append->AddInputConnection(txPolyDataZ->GetOutputPort());
      }

      glyph->SetSourceConnection( append->GetOutputPort() );
      glyph->ScalingOff();
      glyph->OrientOff();
      glyph->Update();
      
    }
    else    
    if( conditionType == VOXFE_BC_FORCE_PARALLEL ) {

      cout << "voxfeGlyphAnnotationFilter found force/parallel BC \n";
      
      glyph->SetSourceConnection(txArrow->GetOutputPort());
      glyph->ScalingOff();
      glyph->OrientOn();
      glyph->SetVectorModeToUseNormal();
      //glyph->SetColorModeToColorByVector();
      glyph->Update();
    }
    else
    if( conditionType == VOXFE_BC_FORCE_TO_POINT ) {

      cout << "voxfeGlyphAnnotationFilter found force/point BC \n";
      
      glyph->SetSourceConnection(txArrow->GetOutputPort());
      glyph->ScalingOff();
      glyph->OrientOn();
      glyph->SetVectorModeToUseNormal();
      //glyph->SetColorModeToColorByVector();
      glyph->Update();

    }
  } 

  //[4]. Output
  //we use vtkMultiBlockDataSet here to pass the geometry through
  vtkInformation *outputModelInfo = outputVector->GetInformationObject(0);
  vtkPolyData  *outputModel = vtkPolyData::SafeDownCast( outputModelInfo->Get(vtkDataObject::DATA_OBJECT()) );
  outputModel->ShallowCopy(glyph->GetOutput());

  //cout << "Added glyphs....\n";

  return 1;
}



////////// External Operators /////////////



