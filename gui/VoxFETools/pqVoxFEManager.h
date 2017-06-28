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

#ifndef __pqVoxFEManager_h
#define __pqVoxFEManager_h

#include <QObject>

#include "ui_pqVoxFEActionHolder.h"

class QAction;
class pqPipelineSource;
class pqRenderView;
class pqServer;
class pqView;

/** \class pqVoxFEManager
    \brief This singleton class manages the callback actions and pipeline handling functions needed
           by VoxFE tools. Note the entry point for the plugin is really pqVoxFEActionGroup, which
           creates the manager and adds the actions held here to paraview.
*/
class pqVoxFEManager : public QObject
{
  Q_OBJECT;
public:
  static pqVoxFEManager *instance();

  ~pqVoxFEManager();

  /** Get the action (callback) for the respective operation. */
  QAction *actionDataLoadManager();

  QAction *actionAddBoundaryCondition();  ///< Add the callback
  QAction *actionOutputSolverScript();    ///< Add the callback
  QAction *actionLoadStrainData();        ///< Add the callback
  QAction *actionHighlightBoundaryConditions();    ///< Add the callback
  QAction *actionExtractBlock();                   ///< Add the callback

  /** Convenience function for getting the current server. */
  pqServer *getActiveServer();

  /** Convenience function for getting the main window. */
  QWidget *getMainWindow();

  /** Get the window used for viewing the mesh. */
  pqView *getMeshView();

  /** Get the window associated with this source */
  pqView *getViewOfSource( pqPipelineSource* src );

  /** Get the renderer used for viewing the mesh. */
  pqRenderView *getMeshRenderView();

  /** Get the reader objects.  Returns NULL if that reader was never created. */
  pqPipelineSource *getMeshReader();

  /** Get the reader objects.  Returns NULL if that reader was never created. */
  pqPipelineSource *getStrainReader();

  /// Get plotting object.  Returns NULL if that object was never created.
  QList<pqPipelineSource*> getBoundaryConditions();

  /** Convenience function for destroying a pipeline object and all of its consumers. */
  static void destroyPipelineSourceAndConsumers(pqPipelineSource *source);

public slots:

  void showDataLoadManager();  ///< Called back function
  void extractBlock();         ///< Called back function
  void checkActionEnabled();   ///< Called back function
  void addBoundaryCondition(); ///< Called back function
  void outputSolverScript();   ///< Called back function
  void loadStrainData();       ///< Called back function
  void generate3DMesh();       ///< Called back function
  void markPointsAsBCs();      ///< Called back function

public slots:

  /** Updates the enabled state. Applications need not explicitly call this.
      (see Qt/ApplicationComponents/pqEditColorMapReaction.h/cxx)
  */
  void updateEnableState();

protected:

  /** Finds a pipeline source with the given SM XML name.  If there is more than
      one, the first is returned.
  */
  virtual pqPipelineSource *findPipelineSource(const char *SMName);

  /** Find all pipeline sources with the given SM XML name. */
  virtual QList<pqPipelineSource*> findAllPipelineSources(const char *SMName);

  /** Finds a view appropriate for the data of the source and port given,
      constrained to those views with the given type.
  */
  virtual pqView *findView(pqPipelineSource *source, int port,
                           const QString &viewType);

  //borrowed from PrismCore
  pqPipelineSource* getActiveSource() const;

 private:

  /** Private constructor */
  pqVoxFEManager(QObject *p);

  /** This widget serves no real purpose other than initializing the ActionHolder
   *  structure created with designer
   */
  Ui::pqVoxFEActionHolder ActionHolder;

  /** This widget serves no real purpose other than initializing the ActionHolder
   *  structure created with designer
   */
  QWidget *ActionHolderWidget;


  pqVoxFEManager(const pqVoxFEManager &);        ///< Not implemented
  void operator=(const pqVoxFEManager &);        ///< Not implemented
};

#endif //__pqVoxFEManager_h
