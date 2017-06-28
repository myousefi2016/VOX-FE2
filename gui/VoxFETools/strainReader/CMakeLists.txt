cmake_minimum_required(VERSION 2.8)
 
#INCLUDE_DIRECTORIES(
#  .
#)
#/home/richard/build/paraview-4.1.0
 
 
FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})
 ADD_PARAVIEW_PLUGIN( strainReader "1.0"
   SERVER_MANAGER_XML strainReader.xml
   SERVER_MANAGER_SOURCES voxfeStrainReader.cxx
   )

