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

#ifndef __pqVoxFEDataLoadManager_h
#define __pqVoxFEDataLoadManager_h

#include <QDialog>

#include "ui_pqVoxFEMeshLoadDialog.h"
#include "ui_pqVoxFEStrainLoadDialog.h"

class pqServer;

/** \class pqVoxFEDataLoadManager
    \brief This dialog box provides an easy way to set up the readers in the pipeline
           and to ready them for the rest of the VoxFE tools.
*/
class pqVoxFEDataLoadManager : public QDialog
{
  Q_OBJECT;

public:

  /** This class is essentially a blank dialog widget; we can use use it to set up alternate
      entry methods by selecting from these
  */
  enum VoxFEDialogType {
    VOXFE_LOAD_MESH,
    VOXFE_LOAD_STRAIN
  };

  /** Constructor */
  pqVoxFEDataLoadManager(QWidget *p, VoxFEDialogType=VOXFE_LOAD_MESH, Qt::WindowFlags f = 0);

  /** Destructor */
  ~pqVoxFEDataLoadManager();

public slots:

  virtual void checkInputValid();   ///< Check eg valid file names
  virtual void setupPipeline();     ///< Set up the ParaView pipline
  virtual void colorMapStrainData();  ///< Begin import of displacement/strain data

signals:

  void createdPipeline();   ///< Signal to emit when pipeline initialized

protected:
  pqServer *Server;            ///< Convenience pointer
  VoxFEDialogType DialogType;  ///< Type: for loading mesh or displacement/strain data

private:
  pqVoxFEDataLoadManager(const pqVoxFEDataLoadManager &);  ///< Not implemented
  void operator=(const pqVoxFEDataLoadManager &);          ///< Not implemented

  /** Voxel data load dialog */
  Ui::pqVoxFEMeshLoadDialog* MeshImportUi;

  /** Displacement/strain data load dialog */
  Ui::pqVoxFEStrainLoadDialog* StrainImportUi;
};

#endif //__pqVoxFEDataLoadManager_h
