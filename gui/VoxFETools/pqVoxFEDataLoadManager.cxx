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


#include "pqVoxFEDataLoadManager.h"
#include "pqVoxFEManager.h"

#include "pqApplicationCore.h"
#include "pqDataRepresentation.h"
#include "pqDisplayPolicy.h"
#include "pqObjectBuilder.h"
#include "pqPipelineSource.h"
#include "pqSMAdaptor.h"
#include "pqUndoStack.h"
#include "vtkSMProperty.h"
#include "vtkSMSourceProxy.h"
#include "pqFiltersMenuReaction.h"
#include "pqActiveObjects.h"

#include <QPushButton>
#include <QtDebug>
#include <QMessageBox>

//=============================================================================
pqVoxFEDataLoadManager::pqVoxFEDataLoadManager(QWidget *p, VoxFEDialogType type,
                                             Qt::WindowFlags f/*=0*/)
  : QDialog(p, f), DialogType(type),
    MeshImportUi(0),
    StrainImportUi(0)
{
  pqVoxFEManager *manager = pqVoxFEManager::instance();
  this->Server = manager->getActiveServer();

  if( DialogType == VOXFE_LOAD_MESH ) {

    this->MeshImportUi = new Ui::pqVoxFEMeshLoadDialog;
    this->MeshImportUi->setupUi(this);

    this->MeshImportUi->meshFile->setServer(this->Server);
    this->MeshImportUi->meshFile->setForceSingleFile(true);
    this->MeshImportUi->meshFile->setExtension("VoxFE Mesh Files (*.script *.voxfe)");

    pqPipelineSource *meshReader = manager->getMeshReader();
    if (meshReader)
      {
      vtkSMProxy *meshReaderProxy = meshReader->getProxy();

      vtkSMProperty *meshFileName = meshReaderProxy->GetProperty("MeshFileName");
      this->MeshImportUi->meshFile->setFilenames(
                                  pqSMAdaptor::getFileListProperty(meshFileName));
      }

    QObject::connect(
                this->MeshImportUi->meshFile, SIGNAL(filenamesChanged(const QStringList &)),
                this, SLOT(checkInputValid()));

    //Initiate the pipeline if user OK's
    QObject::connect(this, SIGNAL(accepted()), this, SLOT(setupPipeline()));

  }
  else if( DialogType == VOXFE_LOAD_STRAIN ) {


    this->StrainImportUi = new Ui::pqVoxFEStrainLoadDialog;
    this->StrainImportUi->setupUi(this);

    this->StrainImportUi->strainFile->setServer(this->Server);
    this->StrainImportUi->strainFile->setForceSingleFile(true);
    this->StrainImportUi->strainFile->setExtension("VoxFE Strain Files (*.strain *.txt)");

    QObject::connect(this, SIGNAL(accepted()), this, SLOT(colorMapStrainData()));
  }

  this->checkInputValid();
}

pqVoxFEDataLoadManager::~pqVoxFEDataLoadManager()
{
  if( this->MeshImportUi ) delete this->MeshImportUi;
  if( this->StrainImportUi ) delete this->StrainImportUi;
}

//-----------------------------------------------------------------------------
void pqVoxFEDataLoadManager::checkInputValid()
{
  bool valid = true;

  if( DialogType == VOXFE_LOAD_MESH ) {

    if (this->MeshImportUi->meshFile->filenames().isEmpty()) valid = false;
    this->MeshImportUi->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(valid);
  }
  else if( DialogType == VOXFE_LOAD_STRAIN ) {

    //fixme:

  }

}

//-----------------------------------------------------------------------------
void pqVoxFEDataLoadManager::setupPipeline()
{
  pqApplicationCore *core = pqApplicationCore::instance();
  pqObjectBuilder *builder = core->getObjectBuilder();
  pqDisplayPolicy *displayPolicy = core->getDisplayPolicy();

  pqVoxFEManager *manager = pqVoxFEManager::instance();

  BEGIN_UNDO_SET("Load VoxFE Mesh Data");

  // Determine the views.  Do this before deleting existing pipeline objects.
  pqView *meshView = manager->getMeshView();

  // Delete existing pipeline objects.  We will replace them.
  manager->destroyPipelineSourceAndConsumers(manager->getMeshReader());
  manager->destroyPipelineSourceAndConsumers(manager->getStrainReader());

  QStringList meshFiles = this->MeshImportUi->meshFile->filenames();
  // This should never really be not empty.
  if (!meshFiles.isEmpty())
    {
    pqPipelineSource *meshReader
      = builder->createReader("sources", "VoxFEReader", meshFiles, this->Server);

    vtkSMSourceProxy* meshReaderProxy =
      vtkSMSourceProxy::SafeDownCast(meshReader->getProxy());

    //#ifdef PUSH_VTK_UPDATES NOTE:
    // These calls have been excluded as they were effectively forcing Auto-Apply on the pipeline
    // But we want to pause while we decide which slices to load
#ifdef PUSH_VTK_UPDATES
    // Push changes to server so that when the representation gets updated,
    // it uses the property values we set.
    meshReaderProxy->UpdateVTKObjects();

    // ensures that new timestep range, if any gets fetched from the server.
    meshReaderProxy->UpdatePipelineInformation();


    // Make representations.
    pqDataRepresentation *repr;
    repr = displayPolicy->createPreferredRepresentation(
                                 meshReader->getOutputPort(0), meshView, false);
    repr->setVisible(true);
    
    // We have already made the representations and pushed everything to the
    // server manager.  Thus, there is no state left to be modified.
    meshReader->setModifiedState(pqProxy::UNMODIFIED);
#endif
  }

  END_UNDO_SET();
  emit this->createdPipeline();

}

void pqVoxFEDataLoadManager::colorMapStrainData() {

  pqApplicationCore *core = pqApplicationCore::instance();
  pqObjectBuilder *builder = core->getObjectBuilder();

  pqVoxFEManager *manager = pqVoxFEManager::instance();
  pqPipelineSource *meshReader = manager->getMeshReader();
  if( !meshReader ) {

    QMessageBox::information(NULL, "pqVoxFEDataLoadManager::colorMapStrainData", "was invoked but meshReader not found\n");
    return;
  }

//See note below
#if 1

  //fixme: this could probably be done more directly via the builder?
  //createFilter accepts the mesh as the model input if it is active, so try to force it....
  pqActiveObjects::instance().setActiveSource(meshReader);
  pqPipelineSource* strainFilter = pqFiltersMenuReaction::createFilter( "filters", "AddStrainsFilter" );

  //this supplies the strain data file
  if(strainFilter) {

    BEGIN_UNDO_SET("Import VoxFE Strain Data");

    vtkSMProxy *strainReaderProxy = strainFilter->getProxy();
    QStringList dispFiles = this->StrainImportUi->strainFile->filenames();
    pqSMAdaptor::setFileListProperty(
                       strainReaderProxy->GetProperty("StrainFile"), dispFiles);

#ifdef PUSH_VTK_UPDATES
     strainReaderProxy->UpdateVTKObjects();
#endif

     END_UNDO_SET();
  }
#else

  // I put this in here as an option to make this into a reader - it works , but you have to set
  // the class as a source in the server xml  and change the input manually, so that this
  // filter also gets the mesh as an input.
  // (The point being that we might have to live without buttons on some machines eg. Archer)

  QStringList dispFiles = this->StrainImportUi->strainFile->filenames();
  if(! dispFiles.isEmpty() ) {

    BEGIN_UNDO_SET("Load Strain Data");

    pqPipelineSource *strainFilter
      = builder->createReader("sources", "AddStrainsFilter", dispFiles, this->Server);

    END_UNDO_SET();
  }
#endif

}

