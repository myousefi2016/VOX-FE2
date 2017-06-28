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
#include "pqConstraintPropertyWidgetDecorator.h"

#include "pqCoreUtilities.h"
#include "pqPropertyWidget.h"
#include "vtkCommand.h"
#include "vtkSMProperty.h"
#include "vtkSMProxy.h"
#include "vtkSMUncheckedPropertyHelper.h"

#include "vtkPVXMLElement.h"

#include <iostream>

//-----------------------------------------------------------------------------
pqConstraintPropertyWidgetDecorator::pqConstraintPropertyWidgetDecorator(
    vtkPVXMLElement* config, pqPropertyWidget* parentObject)
  : Superclass(config, parentObject),ConstraintType("Unknown")
{
  vtkSMProxy* proxy = parentObject->proxy();
  
  vtkSMProperty* propC = proxy? proxy->GetProperty(VOXFE_PANEL_PARAMETER) : NULL;

  if (!propC)
    {
    qDebug("Could not locate property named " VOXFE_PANEL_PARAMETER ". "
      "pqConstraintPropertyWidgetDecorator will have no effect.");
    return;
    }

  this->ObservedObject = propC;
  this->ObserverId = pqCoreUtilities::connect(
    propC, vtkCommand::UncheckedPropertyModifiedEvent,
    this, SIGNAL(visibilityChanged()));

  ConstraintType = config->FindNestedElementByName("Property")->GetAttributeOrDefault("type", "NotFound" );

#ifdef VOXFE_DECORATOR_DEBUG
  std::cout << "\n=====================================================\n";
  std::cout << config->GetName() << "  " << config->GetId() << "\n"; //config->GetAttribute(???) << 
  config->PrintSelf( std::cout, vtkIndent(2) );
  std::cout << "ConstraintType: " << ConstraintType << "\n";
  std::cout << "\n=====================================================\n";
#endif
}

//-----------------------------------------------------------------------------
pqConstraintPropertyWidgetDecorator::~pqConstraintPropertyWidgetDecorator()
{
    if (this->ObservedObject && this->ObserverId)
    {
      this->ObservedObject->RemoveObserver(this->ObserverId);
    }
}

//-----------------------------------------------------------------------------
bool pqConstraintPropertyWidgetDecorator::canShowWidget(bool show_advanced) const
{
  pqPropertyWidget* parentObject = this->parentWidget();
  vtkSMProxy* proxy = parentObject->proxy();
  
  vtkSMProperty* propC = proxy? proxy->GetProperty(VOXFE_PANEL_PARAMETER) : NULL;
   
  if (propC ) {
  
    int valueC = vtkSMUncheckedPropertyHelper(propC).GetAsInt();
#ifdef VOXFE_DECORATOR_DEBUG
    std::cout << "@@ Got Constraint Mode: " << valueC << "  " << proxy->GetXMLName() << " @@ \n";
#endif
    
    if( ConstraintType.compare("axis") == 0 ) {
    
      if( valueC == VOXFE_BC_NODAL ) return true;
    }
    else if( ConstraintType.compare("force_vector") == 0 ) {
    
      if( valueC == VOXFE_BC_FORCE_PARALLEL ) return true;
    }
    else if( ConstraintType.compare("force_endpoint") == 0 ) {
    
      if( valueC == VOXFE_BC_FORCE_TO_POINT ) return true;
    }
    else if( ConstraintType.compare("force_magnitude") == 0 ) {
    
      if( (valueC == VOXFE_BC_FORCE_PARALLEL) || (valueC == VOXFE_BC_FORCE_TO_POINT) ) return true;
    }
    
  }
  
  return false;
}
