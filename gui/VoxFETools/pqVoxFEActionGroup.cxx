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

#include "pqVoxFEActionGroup.h"
#include "pqVoxFEManager.h"

//=============================================================================
pqVoxFEActionGroup::pqVoxFEActionGroup(QObject *p) : QActionGroup(p)
{
  pqVoxFEManager *manager = pqVoxFEManager::instance();
  if (!manager)
    {
    qFatal("Cannot get VoxFE Tools manager.");
    return;
    }

  this->addAction(manager->actionDataLoadManager());
  this->addAction(manager->actionExtractBlock());
  this->addAction(manager->actionAddBoundaryCondition());
  this->addAction(manager->actionHighlightBoundaryConditions());
  this->addAction(manager->actionOutputSolverScript());
  this->addAction(manager->actionLoadStrainData());

  // Action groups are usually used to establish radio-button like
  // functionality.  We don't really want that, so turn it off.
  this->setExclusive(false);
}
