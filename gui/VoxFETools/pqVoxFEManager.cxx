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


/*
  note on building with eclipse
  ccmake -G"Eclipse CDT4 - Unix Makefiles" -D CMAKE_BUILD_TYPE=Debug -i /home/richard/src/async/vox-fe/tools/VoxFETools
*/

#include "pqVoxFEManager.h"
#include "pqVoxFEDataLoadManager.h"

#include "vtkAlgorithm.h"
#include "vtkTable.h"

#include "vtkPVArrayInformation.h"
#include "vtkPVDataInformation.h"
#include "vtkPVDataSetAttributesInformation.h"
#include "vtkSMProxyManager.h"
#include "vtkSMSourceProxy.h"

#include "vtkMath.h"

//rh
#include "pqActiveObjects.h"
#include "pqFiltersMenuReaction.h"
#include "vtkSMSessionProxyManager.h"

//#include "pqActiveView.h"
#include "pqApplicationCore.h"
#include "pqDisplayPolicy.h"
#include "pqObjectBuilder.h"
#include "pqOutputPort.h"
#include "pqPipelineRepresentation.h"
#include "pqPipelineFilter.h"
#include "pqPipelineSource.h"
#include "pqRenderView.h"
#include "pqScalarsToColors.h"
#include "pqServer.h"
#include "pqServerManagerModel.h"
#include "pqSMAdaptor.h"
#include "pqUndoStack.h"
#include "pqXYChartView.h"
#include "pqSMAdaptor.h"

#include "pqSelectionManager.h"

#include "pqActiveObjects.h"

#include <QMainWindow>
#include <QPointer>
#include <QtDebug>
#include <QMessageBox>

//#include "ui_pqVoxFEActionHolder.h"

//=============================================================================
QPointer<pqVoxFEManager> pqVoxFEManagerInstance = NULL;

pqVoxFEManager *pqVoxFEManager::instance()
{
  if (pqVoxFEManagerInstance == NULL)
    {
    pqApplicationCore *core = pqApplicationCore::instance();
    if (!core)
      {
      qFatal("Cannot use the VoxFE Tools without an application core instance.");
      return NULL;
      }

    pqVoxFEManagerInstance = new pqVoxFEManager(core);
    }

  return pqVoxFEManagerInstance;
}

//-----------------------------------------------------------------------------
pqVoxFEManager::pqVoxFEManager(QObject *p) : QObject(p)
{
  // Create the ActionHolder structure (from designer)
  this->ActionHolderWidget = new QWidget(NULL);
  this->ActionHolder.setupUi(this->ActionHolderWidget);

  //Add the buttons in order
  QObject::connect(this->actionDataLoadManager(), SIGNAL(triggered(bool)),
                   this, SLOT(showDataLoadManager()));
  QObject::connect(this->actionExtractBlock(), SIGNAL(triggered(bool)),
                     this, SLOT(extractBlock()));
  QObject::connect(this->actionAddBoundaryCondition(), SIGNAL(triggered(bool)),
                   this, SLOT(addBoundaryCondition()));
  QObject::connect(this->actionOutputSolverScript(), SIGNAL(triggered(bool)),
                   this, SLOT(outputSolverScript()));
  QObject::connect(this->actionLoadStrainData(), SIGNAL(triggered(bool)),
                   this, SLOT(loadStrainData()));
  QObject::connect(this->actionHighlightBoundaryConditions(), SIGNAL(triggered(bool)),
                   this, SLOT(markPointsAsBCs()));

#ifdef CHECK_VOXFE_MESH_IS_LOADED
  this->checkActionEnabled();
#endif  

  //===============================================================
  // Switch off filter buttons when another filter is being applied
  // mirrors behaviour of ParaView as, eg (pqEditColorMapReaction.h/cxx)
  // the filters menu itself is greyed out until the user clicks apply.
  // NB. Unless we do this, it is possible to crash PV by selecting some
  // points, extracting, but before selecting 'Apply', clicking say
  // AddBoundaryCondition
  //===============================================================
  QObject::connect(&pqActiveObjects::instance(),
    SIGNAL(representationChanged(pqDataRepresentation*)),
    this, SLOT(updateEnableState()), Qt::QueuedConnection);  //fixme: QueuedConnection maps across threads, is this needed?
  this->updateEnableState();
}

pqVoxFEManager::~pqVoxFEManager()
{
  delete this->ActionHolderWidget;
}

//-----------------------------------------------------------------------------

void pqVoxFEManager::updateEnableState()
{
  pqPipelineRepresentation* repr = qobject_cast<pqPipelineRepresentation*>(
    pqActiveObjects::instance().activeRepresentation());

  if(repr != NULL) {
    this->actionAddBoundaryCondition()->setEnabled(true);
    this->actionExtractBlock()->setEnabled(true);
    this->actionAddBoundaryCondition()->setEnabled(true);
    this->actionOutputSolverScript()->setEnabled(true);
    this->actionLoadStrainData()->setEnabled(true);
    this->actionHighlightBoundaryConditions()->setEnabled(true);
  }
  else {
    this->actionAddBoundaryCondition()->setEnabled(false);
    this->actionExtractBlock()->setEnabled(false);
    this->actionAddBoundaryCondition()->setEnabled(false);
    this->actionOutputSolverScript()->setEnabled(false);
    this->actionLoadStrainData()->setEnabled(false);
    this->actionHighlightBoundaryConditions()->setEnabled(false);
  }
}


//-----------------------------------------------------------------------------
QAction *pqVoxFEManager::actionDataLoadManager()
{
  return this->ActionHolder.actionDataLoadManager;
}

QAction *pqVoxFEManager::actionAddBoundaryCondition()
{
  return this->ActionHolder.actionAddBoundaryCondition;
}

QAction *pqVoxFEManager::actionOutputSolverScript()
{
  return this->ActionHolder.actionOutputSolverScript;
}

QAction *pqVoxFEManager::actionLoadStrainData()
{
  return this->ActionHolder.actionLoadStrainData;
}

QAction *pqVoxFEManager::actionHighlightBoundaryConditions()
{
  return this->ActionHolder.actionHighlightBoundaryConditions;
}

QAction *pqVoxFEManager::actionExtractBlock()
{
  return this->ActionHolder.actionExtractBlock;
}


//-----------------------------------------------------------------------------
pqServer *pqVoxFEManager::getActiveServer()
{
  pqApplicationCore *app = pqApplicationCore::instance();
  pqServerManagerModel *smModel = app->getServerManagerModel();
  pqServer *server = smModel->getItemAtIndex<pqServer*>(0);
  return server;
}

//-----------------------------------------------------------------------------
QWidget *pqVoxFEManager::getMainWindow()
{
  foreach(QWidget *topWidget, QApplication::topLevelWidgets())
    {
    if (qobject_cast<QMainWindow*>(topWidget)) return topWidget;
    }
  return NULL;
}

//-----------------------------------------------------------------------------
pqView *pqVoxFEManager::findView(pqPipelineSource *source, int port,
                                const QString &viewType)
{
  // Step 1, try to find a view in which the source is already shown.
  if (source)
    {
    foreach (pqView *view, source->getViews())
      {
      pqDataRepresentation *repr = source->getRepresentation(port, view);
      if (repr && repr->isVisible()) return view;
      }
    }

  // Step 2, check to see if the active view is the right type.
  //pqView *view = pqActiveView::instance().current();  //removed in Qt5
  pqView *view = pqActiveObjects::instance().activeView();
  if (view->getViewType() == viewType) return view;

  // Step 3, check all the views and see if one is the right type and not
  // showing anything.
  pqApplicationCore *core = pqApplicationCore::instance();
  pqServerManagerModel *smModel = core->getServerManagerModel();
  foreach (view, smModel->findItems<pqView*>())
    {
    if (   view && (view->getViewType() == viewType)
        && (view->getNumberOfVisibleRepresentations() < 1) )
      {
      return view;
      }
    }

  // Give up.  A new view needs to be created.
  return NULL;
}

pqView *pqVoxFEManager::getMeshView()
{
  return this->findView(this->getMeshReader(), 0,
                        pqRenderView::renderViewType());
}

pqView *pqVoxFEManager::getViewOfSource( pqPipelineSource* src )
{
  return this->findView(src, 0, pqRenderView::renderViewType());
}


pqRenderView *pqVoxFEManager::getMeshRenderView()
{
  return reinterpret_cast<pqRenderView*>(this->getMeshView());
}

/*pqView *pqVoxFEManager::getPlotView()
{
  return this->findView(this->getPlotFilter(), 0,
                        pqXYChartView::XYChartViewType());
}*/

//-----------------------------------------------------------------------------
pqPipelineSource *pqVoxFEManager::findPipelineSource(const char *SMName)
{
  pqApplicationCore *core = pqApplicationCore::instance();
  pqServerManagerModel *smModel = core->getServerManagerModel();

  QList<pqPipelineSource*> sources
    = smModel->findItems<pqPipelineSource*>(this->getActiveServer());

#if 0
  //debug test -- searches through for "name=" tags in XML
  foreach(pqPipelineSource *s, sources)
    {
    cout << "Source: " << s->getProxy()->GetXMLName() << "\n";
    }
#endif

  foreach(pqPipelineSource *s, sources)
    {
    if (strcmp(s->getProxy()->GetXMLName(), SMName) == 0) return s;
    }

  return NULL;
}

QList<pqPipelineSource*> pqVoxFEManager::findAllPipelineSources(const char *SMName){

  pqApplicationCore *core = pqApplicationCore::instance();
  pqServerManagerModel *smModel = core->getServerManagerModel();

  QList<pqPipelineSource*> sources
    = smModel->findItems<pqPipelineSource*>(this->getActiveServer());

  QList<pqPipelineSource*> foundSources;

#if 0
  //debug test -- searches through for "name=" tags in XML
  foreach(pqPipelineSource *s, sources)
    {
    cout << "Source: " << s->getProxy()->GetXMLName() << "\n";
    }
#endif

  foreach(pqPipelineSource *s, sources) {
    if (strcmp(s->getProxy()->GetXMLName(), SMName) == 0) {
      foundSources.push_back( s );
    }
  }

  return foundSources;
}


pqPipelineSource *pqVoxFEManager::getMeshReader()
{
  return this->findPipelineSource("VoxFEReader");
}

pqPipelineSource *pqVoxFEManager::getStrainReader()
{
  return this->findPipelineSource("AddStrainsFilter"); //VoxFEStrainReader"); //fixme: ? xml 'name' here
}

QList<pqPipelineSource*> pqVoxFEManager::getBoundaryConditions()
{
  pqApplicationCore *core = pqApplicationCore::instance();
  pqServerManagerModel *smModel = core->getServerManagerModel();

  QList<pqPipelineSource*> BCs;

  QList<pqPipelineSource*> sources
    = smModel->findItems<pqPipelineSource*>(this->getActiveServer());

  foreach(pqPipelineSource *s, sources)
    {
    if (strcmp(s->getProxy()->GetXMLName(), "AddBoundaryCondition") == 0) BCs << s;
    }

  return BCs;
}


//-----------------------------------------------------------------------------
static void destroyPortConsumers(pqOutputPort *port)
{
  foreach (pqPipelineSource *consumer, port->getConsumers())
    {
    pqVoxFEManager::destroyPipelineSourceAndConsumers(consumer);
    }
}

void pqVoxFEManager::destroyPipelineSourceAndConsumers(pqPipelineSource *source)
{
  if (!source) return;

  foreach (pqOutputPort *port, source->getOutputPorts())
    {
    destroyPortConsumers(port);
    }

  pqApplicationCore *core = pqApplicationCore::instance();
  pqObjectBuilder *builder = core->getObjectBuilder();
  builder->destroy(source);
}

//-----------------------------------------------------------------------------
void pqVoxFEManager::showDataLoadManager()
{
  //this sets up callback to load data/build pipeline
  pqVoxFEDataLoadManager *dialog = new pqVoxFEDataLoadManager(this->getMainWindow());

  dialog->setAttribute(Qt::WA_DeleteOnClose, true);
  
#ifdef CHECK_VOXFE_MESH_IS_LOADED 
  //removed at present as loading from saved state does not trigger
  QObject::connect(dialog, SIGNAL(createdPipeline()), this, SLOT(checkActionEnabled()));
#endif  

  //QObject::connect(dialog, SIGNAL(createdPipeline()), this, SLOT(showEField()));
  //QObject::connect(dialog, SIGNAL(createdPipeline()), this, SLOT(showStandardViewpoint()));
  dialog->show();
}

//-----------------------------------------------------------------------------
void pqVoxFEManager::checkActionEnabled()
{
  //QMessageBox::information(NULL, "pqVoxFEManager", "checkActionEnabled was invoked\n");

  pqPipelineSource *meshReader = this->getMeshReader();
  pqPipelineSource *strainReader = this->getStrainReader();

  if (!meshReader)
    {
    this->actionAddBoundaryCondition()->setEnabled(false);
    this->actionOutputSolverScript()->setEnabled(false);
    this->actionLoadStrainData()->setEnabled(false);
    //this->actionGenerate3DMesh()->setEnabled(false);
    this->actionHighlightBoundaryConditions()->setEnabled(false);

    //cerr << "**************************************************\n";
    //cerr << "Disabling buttons\n";
    //cerr << "**************************************************\n";
    }
  else
    {
    this->actionAddBoundaryCondition()->setEnabled(true);
    this->actionOutputSolverScript()->setEnabled(true);
    this->actionLoadStrainData()->setEnabled(true);
    //this->actionGenerate3DMesh()->setEnabled(true);
    this->actionHighlightBoundaryConditions()->setEnabled(true);
  }

}

//-----------------------------------------------------------------------------

void pqVoxFEManager::addBoundaryCondition()
{
  //borrowed from PrismCore.cxx
  // Get the list of selected sources.
  pqApplicationCore* core = pqApplicationCore::instance();
  pqObjectBuilder* builder = core->getObjectBuilder();
  pqPipelineSource* source = 0;
  pqPipelineSource* filter_merge = 0;
  pqPipelineSource* filter_addBC = 0;
  pqServer* server = 0;
  QList<pqOutputPort*> inputs;

  source = this->getActiveSource();
  if(!source)
  {
      QMessageBox::warning(NULL, tr("No Object Selected"),
          tr("No pipeline object is selected.\n"
          "Please select a pointset object from the list on the left."),
          QMessageBox::Ok );

      return;
  }
  else {

	  const char* obj_selected = source->getProxy()->GetXMLName();

	  //Allow ExtractSelection and Clip objects only.... 
	  // (strstr is null if can't find substring)
	  if( !strstr( obj_selected, "ExtractSelection" ) && !strstr( obj_selected, "Clip" ) ) {  

		QString s("Selected object: ");
		s.append(obj_selected);

	    //check name of object
	    QMessageBox::warning(NULL, tr(s.toLocal8Bit().data()),
            tr("You must extract selected points (and click apply) \n"
			   "before attempting to add boundary conditions."),
            QMessageBox::Ok );

		return;
	  }
	  else {

 	   //fixme: Is there a way to tell if 'Apply' has been pressed??
		/*  ... this is protected anyway ...
		bool modified = source->getProxy()->ArePropertiesModified();
		if( modified ) 
			 QMessageBox::warning(NULL, tr("Selected object..."),
            tr("has been modified"),
            QMessageBox::Ok );
			
		//source->getProxy()->UpdateVTKObjects();
		source->getProxy()->UpdateSelfAndAllInputs();
		vtkSMSourceProxy* sourceProxy = vtkSMSourceProxy::SafeDownCast(source->getProxy());
		sourceProxy->UpdatePipeline(); //force apply

		return;
		*/
	  }


  }


  server = source->getServer();
  if(!server)
  {
    qDebug() << "No active server selected.";
  }

  inputs.push_back(source->getOutputPort(0));

  QMap<QString, QList<pqOutputPort*> > namedInputs;
  namedInputs["Input"] = inputs;

  //debug
  cout << "\n\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
       << "\nAddBC has received: " << namedInputs.size() << " inputs, with: "
       << source->getNumberOfOutputPorts() << " outputs\n";
  cout << "Port 0 has: " << source->getNumberOfConsumers(0) << " consumers.....\n";

  pqUndoStack *stack=core->getUndoStack();

  //-----------------------------------------------------------------------------------------
  if(stack)
  {
    stack->beginUndoSet("Add Boundary Condition");
  }

  filter_merge = builder->createFilter("filters", "MergeBlocks", namedInputs, server);
  //fixme: is there a way to hide this filter in the browser???

  inputs.clear();
  namedInputs.clear();
  unsigned int outport = 0;
  inputs.push_back( filter_merge->getOutputPort(0) );
  namedInputs["Input"] = inputs;

  //-----------------------------------------------------------------------------------------

  filter_addBC = builder->createFilter("filters", "AddBoundaryCondition", namedInputs, server);

  if(stack)
  {
    stack->endUndoSet();
  }

  if( filter_addBC ) {

    //fixme: This attempts to remove the selection, to prevent it getting passed on to another block
    // (though I think it may still exist in the system....)
    pqSelectionManager* slmanager = qobject_cast<pqSelectionManager*>(
      core->manager("SelectionManager"));
    if( slmanager )
      slmanager->clearSelection();
  }

  //==================================================

  //fixme: this works, but can be annoying if you're working with individual blocks
  // can we adjust accordingly??

  //Switch on the 'ExtractBlock' if found or the 'MeshReader' if not....
  pqPipelineSource *meshObj = this->findPipelineSource("ExtractBlock");
  if( !meshObj ) meshObj    = this->getMeshReader();
  if(meshObj) {

	pqView *view = this->getViewOfSource( meshObj );
	if(view) {

	  pqDataRepresentation *repr = meshObj->getRepresentation(0, view);
      if(repr) {
		if( !repr->isVisible() ) {
		  repr->setVisible(true);
		  view->render();
		}
	  }
	}
  }

  //===================================================
  //fixme: This is a bunch of test code where I was trying to update the gui with the voxel size
  //as a step size for adjusting sliders. I think, reading through Markmail, this is not
  //directly possible - we'd have to create a custom widget and allow access somehow to that
  // eg http://paraview.markmail.org/search/?q=update%20pipeline#query:update%20pipeline+page:2+mid:jeech4otraze2blx+state:results
  /*
  bool testqv;
  vtkSMSourceProxy* sourceProxy = vtkSMSourceProxy::SafeDownCast(this->getMeshReader()->getProxy());
  sourceProxy->UpdateProperty("MeshSize");
  sourceProxy->UpdatePropertyInformation( sourceProxy->GetProperty("MeshSize") );
  sourceProxy->UpdatePipelineInformation();
  sourceProxy->UpdatePipeline();
  QVariant qv = pqSMAdaptor::getElementProperty( sourceProxy->GetProperty("MeshSize") );
  double voxel_size = qv.toDouble(&testqv);
  if( !testqv ) {
    cerr << "Cannot obtain voxel size ... \n";
  }
  else {


    cout << "\n\n**Got MESH SIZE: " << voxel_size << " ** \n\n";


  }
  // --------
  meshObj = ;
  if( meshObj ) {

    vtkSMProxy *meshReaderProxy = meshReader->getProxy();
    QVariant qv;
	pqSMAdaptor::getElementProperty( meshReaderProxy->GetProperty("MeshSize"), qv);
	  
	  
	  
pqDataRepresentation *repr = filter_addBC->getRepresentation(0, view);
      if (!repr) {

        vtkSMProxy *reprProxy = repr->getProxy();
        pqSMAdaptor::setEnumerationProperty( reprProxy->GetProperty("Representation"), "Points");

        //view->render();
      }
  }
  */
}


void pqVoxFEManager::outputSolverScript()
{
  //borrowed from PrismCore.cxx
  // Get the list of selected sources.
  pqApplicationCore* core = pqApplicationCore::instance();
  pqObjectBuilder* builder = core->getObjectBuilder();
  pqPipelineSource* source = 0;
  pqPipelineSource* filter = 0;
  pqServer* server = 0;
  QList<pqOutputPort*> inputs;

  //fixme: most of the work here is being done by pqFiltersMenuReaction::createFilter below
  // -- can we make smarter use of this stuff or remove
  source = this->getActiveSource();
  if(!source /*|| (source->getNumberOfOutputPorts()<2)*/ )
  {
      QMessageBox::warning(NULL, tr("No Objects Selected"),
          tr("At least 2 pipeline objects must be selected.\n"
          "Please select the VoxFE mesh and constraint point data from the list on the left."),
          QMessageBox::Ok );

      return;
  }
  server = source->getServer();
  if(!server)
  {
    qDebug() << "No active server selected.";
  }

#if 0
  //FIXME: solver script has 2 or more inputs.....???? --- Change Input Dialog  ????
  //see pqChangeInputDialog  OR maybe void pqChangePipelineInputReaction::changeInput()
  //
  inputs = source->getOutputPorts();
  QMap<QString, QList<pqOutputPort*> > namedInputs;
  namedInputs["Input"] = inputs;


  //test
  std:: cout << "\n\n***************************\nOutputScript has ports:\n";
  for( int k=0; k<source->getNumberOfOutputPorts(); ++k ) {
    std::cout << "\t" << inputs[k]->getPortName().toStdString() << std::endl;
  }
  std:: cout << "\n****************************\n\n";
#endif


  //filter = builder->createFilter("filters", "OutputSolverScript", namedInputs, server);
  pqActiveObjects::instance().setActiveSource(this->getMeshReader());
  filter = pqFiltersMenuReaction::createFilter( "filters", "OutputSolverScript" );


}


void pqVoxFEManager::loadStrainData() {

  //this sets up callback to load data/build pipeline
  pqVoxFEDataLoadManager *dialog = new pqVoxFEDataLoadManager(this->getMainWindow(), pqVoxFEDataLoadManager::VOXFE_LOAD_STRAIN);

  dialog->setAttribute(Qt::WA_DeleteOnClose, true);

  //QObject::connect(dialog, SIGNAL(createdPipeline()), this, SLOT(checkActionEnabled()));

  dialog->show();

}


//Depracated  as superceded by GenerateSurfaceMesh option in panel
void pqVoxFEManager::generate3DMesh() {

  pqApplicationCore* core = pqApplicationCore::instance();
  pqObjectBuilder* builder = core->getObjectBuilder();
  pqPipelineSource* source = 0;
  pqPipelineSource* filter = 0;
  pqServer* server = 0;
  QList<pqOutputPort*> inputs;

  source = pqActiveObjects::instance().activeSource(); //this->getActiveSource();
  if(!source)
  {
      QMessageBox::warning(NULL, tr("No Object Selected"),
          tr("No pipeline object is selected.\n"
          "Please select a pointset object from the list on the left."),
          QMessageBox::Ok );

      return;
  }
  server = source->getServer();
  if(!server)
  {
    qDebug() << "No active server selected.";
  }

  inputs.push_back(source->getOutputPort(0));

  QMap<QString, QList<pqOutputPort*> > namedInputs;
  namedInputs["Input"] = inputs;

  pqUndoStack *stack=core->getUndoStack();
  if(stack)
  {
    stack->beginUndoSet("Remodel Point Data (Delaunay3D)");
  }

  filter = builder->createFilter("filters", "Delaunay3D", namedInputs, server);

  if(stack)
  {
    stack->endUndoSet();
  }

}

void pqVoxFEManager::extractBlock() {

  pqApplicationCore* core = pqApplicationCore::instance();
  pqObjectBuilder* builder = core->getObjectBuilder();
  pqPipelineSource* source = 0;
  pqPipelineSource* filter = 0;
  pqServer* server = 0;
  QList<pqOutputPort*> inputs;

  source = pqActiveObjects::instance().activeSource(); //this->getActiveSource();
  if(!source)
  {
      QMessageBox::warning(NULL, tr("No Object Selected"),
          tr("No pipeline object is selected.\n"
          "Please select a pointset object from the list on the left."),
          QMessageBox::Ok );

      return;
  }
  server = source->getServer();
  if(!server)
  {
    qDebug() << "No active server selected.";
  }

  inputs.push_back(source->getOutputPort(0));

  QMap<QString, QList<pqOutputPort*> > namedInputs;
  namedInputs["Input"] = inputs;

  pqUndoStack *stack=core->getUndoStack();
  if(stack)
  {
    stack->beginUndoSet("Extract Block");
  }

  filter = builder->createFilter("filters", "ExtractBlock", namedInputs, server);

  //we'd like to switch off the prune output to stop the scattering of data cf Utkarsh's
  // after my question:
  /* Try unchecking the "Purne Output" checkbox on the "Extract Block"
     filter's Property panel. The problem is that "Prune Output" changes
     the structure of the multiblock. Hence the selection defined in (3)
     resolves differently after (5). Utkarsh
   */
  vtkSMSourceProxy* proxy = vtkSMSourceProxy::SafeDownCast( filter->getProxy() );
  QVariant qv = 0;
  pqSMAdaptor::setElementProperty( proxy->GetProperty("PruneOutput", 0), qv );

  if(stack)
  {
    stack->endUndoSet();
  }

}


//-----------------------------------------------------------------------------

//Loop through BCs and mark with bigger point glyphs
// (modified from SLAC::showSolidMesh())
void pqVoxFEManager::markPointsAsBCs()
{
  //fixme: make members?
  QVariant ptsize = 5.0;
  double colour_step = 0.5;

  pqApplicationCore *core = pqApplicationCore::instance();
  pqUndoStack *stack = core->getUndoStack();

  QList<pqPipelineSource*> sourceBCs = this->getBoundaryConditions();

  if( stack && (sourceBCs.size() > 0) )
    stack->beginUndoSet("Mark points as BCs");
  else return;

  int k = 0;
  foreach(pqPipelineSource *src, sourceBCs)
  {
    if( !src ) {
      cerr << "Cannot find src object.... \n";
      continue;
    }

    pqView *view = this->findView(src, 0, pqRenderView::renderViewType());
    if (!view) {
      cerr << "Cannot find view... \n";
      continue;
    }

    pqDataRepresentation *repr = src->getRepresentation(0, view);
    if (!repr) {
      cerr << "Cannot find representation... \n";
      continue;
    }
    vtkSMProxy *reprProxy = repr->getProxy();

    //increase the point size (default=2)
    pqSMAdaptor::setEnumerationProperty( reprProxy->GetProperty("Representation"), "Points");
    pqSMAdaptor::setElementProperty( reprProxy->GetProperty("PointSize",0), ptsize );

    //find the BC type
    bool testqv;
    vtkSMSourceProxy* sourceProxy = vtkSMSourceProxy::SafeDownCast(src->getProxy());
    QVariant qv = pqSMAdaptor::getElementProperty( sourceProxy->GetProperty("BoundaryCondition") );
    int bctype = qv.toInt(&testqv);
    if( !testqv ) {
      cerr << "Cannot obtain BC type... \n";
      continue;
    }

	//try to randomize the colour a bit
	k++;
	double step = k*colour_step + vtkMath::Random(-0.25,0.25); 
	if( (step > 1.0) || (step < 0.0) ) {
	  k = 0;
	  step = vtkMath::Random(0.0,1.0);
	}
	QList< QVariant > nodalConstraintColor, forceParallelColor, forceToPointColor, noRemodelColor;
    nodalConstraintColor << (1.0-step) << (0.5+vtkMath::Random(-0.3,0.5))  << step;
    forceParallelColor   << (0.5+vtkMath::Random(-0.3,0.5))  << step << (1.0-step);
    forceToPointColor    << step << (1.0-step) << (0.5+vtkMath::Random(-0.3,0.5));
    noRemodelColor       << (1.0-step) << step << (0.5+vtkMath::Random(-0.3,0.5));

    //types hard-coded in SM XML
    if( bctype == 1 )
      pqSMAdaptor::setMultipleElementProperty( reprProxy->GetProperty("AmbientColor"), nodalConstraintColor );
    else if( bctype == 2 )
      pqSMAdaptor::setMultipleElementProperty( reprProxy->GetProperty("AmbientColor"), forceParallelColor );
    else if( bctype == 3 )
      pqSMAdaptor::setMultipleElementProperty( reprProxy->GetProperty("AmbientColor"), forceToPointColor );
    else if( bctype == 4 )
      pqSMAdaptor::setMultipleElementProperty( reprProxy->GetProperty("AmbientColor"),  noRemodelColor );


    reprProxy->UpdateVTKObjects();  //fixme: needed to force update -- is there a better way?
                                    // - and - maybe just do once ...? (left here in case mpi)

    view->render();
  }

  if (stack) stack->endUndoSet();


  //Switch off the 'MergeBlocks' to avoid this data muddying the highlighted data
  QList<pqPipelineSource*> mergeObj = this->findAllPipelineSources("MergeBlocks");

  foreach(pqPipelineSource *s, mergeObj) {

    //cout << "Source: " << s->getProxy()->GetXMLName() << "\n";
    pqView *view = this->getViewOfSource( s );
    if (!view) continue;

    pqDataRepresentation *repr = s->getRepresentation(0, view);
    if (!repr) continue;

    if( repr->isVisible() ) {
      repr->setVisible(false);
      view->render();                    //FIXME: are all these **separate** render calls needed....???!!
    }
  }

}


//----------------------------------------------------------------------------
pqPipelineSource* pqVoxFEManager::getActiveSource() const
{
  return pqActiveObjects::instance().activeSource();
}


//=====================================================================================================================
//=====================================================================================================================
//=====================================================================================================================

/*
void pqVoxFEManager::showWireframeSolidMesh()
{
  pqPipelineSource *reader = this->getMeshReader();
  if (!reader) return;

  pqView *view = this->getMeshView();
  if (!view) return;

  pqDataRepresentation *repr = reader->getRepresentation(0, view);
  if (!repr) return;
  vtkSMProxy *reprProxy = repr->getProxy();

  pqApplicationCore *core = pqApplicationCore::instance();
  pqUndoStack *stack = core->getUndoStack();

  if (stack) stack->beginUndoSet("Show Wireframe Mesh");

  pqSMAdaptor::setEnumerationProperty(
                reprProxy->GetProperty("Representation"), "Surface With Edges");
  pqSMAdaptor::setEnumerationProperty(
          reprProxy->GetProperty("BackfaceRepresentation"), "Follow Frontface");

  reprProxy->UpdateVTKObjects();

  if (stack) stack->endUndoSet();

  view->render();
}*/

/*
void pqVoxFEManager::showWireframeAndBackMesh()
{
  pqPipelineSource *reader = this->getMeshReader();
  if (!reader) return;

  pqView *view = this->getMeshView();
  if (!view) return;

  pqDataRepresentation *repr = reader->getRepresentation(0, view);
  if (!repr) return;
  vtkSMProxy *reprProxy = repr->getProxy();

  pqApplicationCore *core = pqApplicationCore::instance();
  pqUndoStack *stack = core->getUndoStack();

  if (stack) stack->beginUndoSet("Show Wireframe Front and Solid Back");

  pqSMAdaptor::setEnumerationProperty(
                         reprProxy->GetProperty("Representation"), "Wireframe");
  pqSMAdaptor::setEnumerationProperty(
                   reprProxy->GetProperty("BackfaceRepresentation"), "Surface");

  reprProxy->UpdateVTKObjects();

  if (stack) stack->endUndoSet();

  view->render();
}*/

//-----------------------------------------------------------------------------
/*
void pqVoxFEManager::createPlotOverZ()
{
  pqApplicationCore *core = pqApplicationCore::instance();
  pqObjectBuilder *builder = core->getObjectBuilder();
  pqUndoStack *stack = core->getUndoStack();
  pqDisplayPolicy *displayPolicy = core->getDisplayPolicy();

  pqPipelineSource *meshReader = this->getMeshReader();
  if (!meshReader) return;

  if (stack) stack->beginUndoSet("Plot Over Z");

  // Determine view.  Do this before deleting existing pipeline objects.
  pqView *plotView = this->getPlotView();

  // Delete existing plot objects.  We will replace them.
  this->destroyPipelineSourceAndConsumers(this->getPlotFilter());

  // Turn on reading the internal volume, which is necessary for plotting
  // through the volume.
  vtkSMProxy *meshReaderProxy = meshReader->getProxy();
  pqSMAdaptor::setElementProperty(
                      meshReaderProxy->GetProperty("ReadInternalVolume"), true);
  meshReaderProxy->UpdateVTKObjects();
  meshReader->updatePipeline();

  // Get the mesh data bounds (which we will use later to set up the plot).
  vtkPVDataInformation *dataInfo
    = meshReader->getOutputPort(1)->getDataInformation();
  double bounds[6];
  dataInfo->GetBounds(bounds);

  // Create the plot filter.
  QMap<QString, QList<pqOutputPort*> > namedInputs;
  QList<pqOutputPort *> inputs;
  inputs.push_back(meshReader->getOutputPort(1));
  namedInputs["Input"] = inputs;
  pqPipelineSource *plotFilter = builder->createFilter("filters", "ProbeLine",
                                                       namedInputs,
                                                       this->getActiveServer());

  // Set up the line for the plot.  The line is one of the inputs to the filter
  // which is implicitly set up through a proxy list domain.
  vtkSMProxy *plotProxy = plotFilter->getProxy();
  pqSMProxy lineProxy
    = pqSMAdaptor::getProxyProperty(plotProxy->GetProperty("Source"));
  if (!lineProxy)
    {
    qWarning() << "Could not retrieve plot line source.  "
               << "Plot not set up correctly.";
    }
  else
    {
    QList<QVariant> minPoint;
    minPoint << 0.0 << 0.0 << bounds[4];
    pqSMAdaptor::setMultipleElementProperty(lineProxy->GetProperty("Point1"),
                                            minPoint);
    QList<QVariant> maxPoint;
    maxPoint << 0.0 << 0.0 << bounds[5];
    pqSMAdaptor::setMultipleElementProperty(lineProxy->GetProperty("Point2"),
                                            maxPoint);
    pqSMAdaptor::setElementProperty(lineProxy->GetProperty("Resolution"), 1000);
    lineProxy->UpdateVTKObjects();
    }

  // Make representation
  pqDataRepresentation *repr;
  repr = displayPolicy->createPreferredRepresentation(
                                 plotFilter->getOutputPort(0), plotView, false);
  repr->setVisible(true);

  this->updatePlotField();

  // We have already made the representations and pushed everything to the
  // server manager.  Thus, there is no state left to be modified.
  meshReader->setModifiedState(pqProxy::UNMODIFIED);
  plotFilter->setModifiedState(pqProxy::UNMODIFIED);

  if (stack) stack->endUndoSet();
}*/


//-----------------------------------------------------------------------------
/*
void pqVoxFEManager::toggleBackgroundBW()
{
  pqRenderView *view = this->getMeshRenderView();
  if (!view) return;
  vtkSMProxy *viewProxy = view->getProxy();

  QList<QVariant> oldBackground;
  QList<QVariant> newBackground;

  oldBackground = pqSMAdaptor::getMultipleElementProperty(
                                          viewProxy->GetProperty("Background"));
  if (   (oldBackground[0].toDouble() == 0.0)
      && (oldBackground[1].toDouble() == 0.0)
      && (oldBackground[2].toDouble() == 0.0) )
    {
    newBackground << 1.0 << 1.0 << 1.0;
    }
  else if (   (oldBackground[0].toDouble() == 1.0)
           && (oldBackground[1].toDouble() == 1.0)
           && (oldBackground[2].toDouble() == 1.0) )
    {
    const int *defaultBackground = view->defaultBackgroundColor();
    newBackground << defaultBackground[0]/255.0
                  << defaultBackground[1]/255.0
                  << defaultBackground[2]/255.0;
    }
  else
    {
    newBackground << 0.0 << 0.0 << 0.0;
    }

  pqSMAdaptor::setMultipleElementProperty(viewProxy->GetProperty("Background"),
                                          newBackground);

  viewProxy->UpdateVTKObjects();
  view->render();
}*/

//-----------------------------------------------------------------------------
/*
void pqVoxFEManager::showStandardViewpoint()
{
  pqRenderView *view = qobject_cast<pqRenderView*>(this->getMeshView());
  if (view)
    {
    view->resetViewDirection(1, 0, 0,
                             0, 1, 0);
    }
  view->render();
}*/

//-----------------------------------------------------------------------------
/*
void pqVoxFEManager::resetRangeTemporal()
{
  this->ScaleFieldsByCurrentTimeStep = false;

  // Check to see if the ranges are already computed.
  if (this->getTemporalRanges())
    {
    this->showField(this->CurrentFieldName);
    return;
    }

  pqApplicationCore *core = pqApplicationCore::instance();
  pqObjectBuilder *builder = core->getObjectBuilder();
  pqUndoStack *stack = core->getUndoStack();

  pqPipelineSource *meshReader = this->getMeshReader();
  if (!meshReader) return;

  if (stack) stack->beginUndoSet("Compute Ranges Over Time");

  // Turn on reading the internal volume, which is necessary for plotting
  // through the volume.
  vtkSMProxy *meshReaderProxy = meshReader->getProxy();
  pqSMAdaptor::setElementProperty(
                      meshReaderProxy->GetProperty("ReadInternalVolume"), true);
  meshReaderProxy->UpdateVTKObjects();
  meshReader->updatePipeline();

  // Create the temporal ranges filter.
  pqPipelineSource *rangeFilter = builder->createFilter("filters",
                                                        "TemporalRanges",
                                                        meshReader, 1);

  this->showField(this->CurrentFieldName);

  // We have already pushed everything to the server manager, and I don't want
  // to bother making representations.  Thus, it is unnecessary to make any
  // further modifications.
  meshReader->setModifiedState(pqProxy::UNMODIFIED);
  rangeFilter->setModifiedState(pqProxy::UNMODIFIED);

  if (stack) stack->endUndoSet();
}*/

//-----------------------------------------------------------------------------
/*
void pqVoxFEManager::resetRangeCurrentTime()
{
  this->ScaleFieldsByCurrentTimeStep = true;
  this->showField(this->CurrentFieldName);
}*/

