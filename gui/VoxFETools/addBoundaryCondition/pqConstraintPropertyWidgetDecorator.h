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

#ifndef __pqConstraintPropertyWidgetDecorator_h
#define __pqConstraintPropertyWidgetDecorator_h

#include "pqPropertyWidgetDecorator.h"
#include "vtkWeakPointer.h"
#include <string>

class vtkObject;


/** \class pqConstraintPropertyWidgetDecorator
 *  \brief Adjust the ParaView gui controls according to the Constraint mode selected
 *
 */
class pqConstraintPropertyWidgetDecorator : public pqPropertyWidgetDecorator
{
  Q_OBJECT
  typedef pqPropertyWidgetDecorator Superclass;
public:

  /** Constructor */
  pqConstraintPropertyWidgetDecorator( vtkPVXMLElement* config, pqPropertyWidget* parentObject );

  /** Destructor */
  virtual ~pqConstraintPropertyWidgetDecorator();

  virtual bool canShowWidget(bool show_advanced) const;   ///< Can we show the widget acc to supplied property (axis or force) and type?

private:
  Q_DISABLE_COPY(pqConstraintPropertyWidgetDecorator)

  vtkWeakPointer<vtkObject> ObservedObject;
  unsigned long ObserverId;
  std::string ConstraintType;
  
};

#endif
