# This plugin contains tools developed for the VoxFE voxel-based
# FE package.
#
# Its layout borrows from other plugins, notably SLACTools, but
# we also envisage releasing separate solver and remodelling tools.
# This CMake contains instructions for the main GUI component, which
# adds a number of buttons to the ParaView interface for loading
# models, adding boundary conditions, displaying strains etc.
#
# VoxFE Team, University of Hull, Hull, UK

cmake_minimum_required(VERSION 2.8)

SET(PLUGIN_NAME VoxFETools)
SET(PLUGIN_VERSION "0.9")

#find paraview build
IF(NOT ParaView_BINARY_DIR)
  FIND_PACKAGE(ParaView REQUIRED)
  INCLUDE(${PARAVIEW_USE_FILE})
ENDIF(NOT ParaView_BINARY_DIR)

OPTION(VOXFE_USE_ITK "Use ITK library to allow reading image data" OFF)

#find itk build, if wanted
IF(VOXFE_USE_ITK)
  FIND_PACKAGE(ITK REQUIRED)
  IF(ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})
  ENDIF(ITK_FOUND)
  
  ADD_DEFINITIONS(-DVOXFE_USE_ITK)
  
ENDIF(VOXFE_USE_ITK)  

INCLUDE_DIRECTORIES(
  ${VTK_INCLUDE_DIRS}
  ./itkReader
  ./addBoundaryCondition
  ./outputScript
  ./strainReader
  )

SET(SM_XML
  VoxFEToolsSM.xml
  )

SET(SM_SRC
  itkReader/voxfeITKReader.cxx
  addBoundaryCondition/voxfeAddConstraintFilter.cxx
  addBoundaryCondition/voxfeGlyphAnnotationFilter.cxx
  outputScript/voxfeOutputScriptFilter.cxx
  strainReader/voxfeStrainFilter.cxx
  )

SET(UI_FILES
  pqVoxFEActionHolder.ui
  pqVoxFEMeshLoadDialog.ui
  pqVoxFEStrainLoadDialog.ui
  )


SET(MOC_HEADERS
  pqVoxFEActionGroup.h
  pqVoxFEDataLoadManager.h
  pqVoxFEManager.h 
  )

SET(CLIENT_SRCS
  pqVoxFEActionGroup.cxx
  pqVoxFEDataLoadManager.cxx
  pqVoxFEManager.cxx
  )
  
SET(ADDITIONAL_SRCS  
  strainReader/voxfeStrainData.cxx
  strainReader/voxfeRecomputeStrainFilter.cxx
  itkReader/ReadVoxFEScript.cxx  
  )
  
# Slight hack... added to allow Windows to find pqFiltersMenuReaction::createFilter
IF(WIN32)
  FIND_LIBRARY(VOXFE_ADDITIONAL_LIBS 
    NAMES vtkpqApplicationComponents-pv${PARAVIEW_VERSION}
    # making an assumption here, but why else would you be building if not to debug
    #	PATHS ${ParaView_DIR}/lib/Debug               
    # this also chooses Debug by default, need to change manually in VS if not
    PATHS ${ParaView_DIR}/lib/${CMAKE_CFG_INTDIR}     
  )
ENDIF(WIN32)


LINK_LIBRARIES( ${PLUGIN_NAME}  ${ITK_LIBRARIES}  ${VOXFE_ADDITIONAL_LIBS} )

IF (PARAVIEW_BUILD_QT_GUI)

  IF (PARAVIEW_QT_VERSION VERSION_GREATER "4")
    SET (Qt5_FIND_COMPONENTS Widgets)
    INCLUDE (ParaViewQt5)
    QT5_WRAP_CPP(MOC_SRCS ${MOC_HEADERS} addBoundaryCondition/pqConstraintPropertyWidgetDecorator.h)
    QT5_WRAP_UI(UI_SRCS ${UI_FILES})
  ELSE ()
    INCLUDE(${QT_USE_FILE})
    QT4_WRAP_CPP(MOC_SRCS ${MOC_HEADERS} addBoundaryCondition/pqConstraintPropertyWidgetDecorator.h)
    QT4_WRAP_UI(UI_SRCS ${UI_FILES})
  ENDIF ()
  
  SET (outifaces2)
  SET (outsrcs2)
  ADD_PARAVIEW_PROPERTY_WIDGET_DECORATOR(outifaces2 outsrcs2
    TYPE "constraint_decorator"
    CLASS_NAME pqConstraintPropertyWidgetDecorator)
  
  ADD_PARAVIEW_ACTION_GROUP(VoxFE_TOOLBAR_IFACE VoxFE_TOOLBAR_IFACE_SRCS
    CLASS_NAME pqVoxFEActionGroup
    GROUP_NAME "ToolBar/VoxFE"
    )

  ADD_PARAVIEW_PLUGIN(${PLUGIN_NAME} ${PLUGIN_VERSION}
    DOCUMENTATION_DIR "${CMAKE_CURRENT_SOURCE_DIR}/plugin-doc"
    SERVER_MANAGER_XML ${SM_XML}
    SERVER_MANAGER_SOURCES ${SM_SRC}
    GUI_INTERFACES ${VoxFE_TOOLBAR_IFACE}  ${outifaces2}
    GUI_RESOURCES VoxFETools.qrc
    SOURCES ${UI_SRCS} ${VoxFE_TOOLBAR_IFACE_SRCS} ${CLIENT_SRCS} ${ADDITIONAL_SRCS}
    GUI_SOURCES ${outsrcs2} ${MOC_SRCS} addBoundaryCondition/pqConstraintPropertyWidgetDecorator.cxx
    )

ELSE (PARAVIEW_BUILD_QT_GUI)
  ADD_PARAVIEW_PLUGIN(${PLUGIN_NAME} ${PLUGIN_VERSION}
    SERVER_MANAGER_XML ${SM_XML}
    SERVER_MANAGER_SOURCES ${SM_SRC}
    )
ENDIF (PARAVIEW_BUILD_QT_GUI)
