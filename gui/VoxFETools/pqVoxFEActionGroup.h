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

#ifndef __pqVoxFEActionGroup_h
#define __pqVoxFEActionGroup_h

#include <QActionGroup>

/** \class pqVoxFEActionGroup
 *  \brief Adds actions that are helpful for setting up visualization of VoxFE simulation result files.
 */
class pqVoxFEActionGroup : public QActionGroup
{
  Q_OBJECT;
public:

  /** Constructor */
  pqVoxFEActionGroup(QObject *p);

private:
  pqVoxFEActionGroup(const pqVoxFEActionGroup &);    ///< Not implemented
  void operator=(const pqVoxFEActionGroup &);        ///< Not implemented
};

#endif //__pqVoxFEActionGroup_h
