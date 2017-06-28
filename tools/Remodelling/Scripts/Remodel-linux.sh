#! /bin/bash

# This is a first attempt to perform remodelling by calling all the following routines 
# in sequence and saving the product of each loop in a subdirectory. 
# 
# We rely on 4 routines here:
#
# 1. Solver generates file "displacements.txt"
#
# 2. convertVoxFEScript creates a graph/list of the nodes in each element and a separate
#    element-map ('.elmap') which lists the global IDs of all the elements in the model
#
# 3. m2gmetis creates a graph/list of elements which is 1-based (we need to refer back to
#    the .elmap file to get global IDs
#
# 4. voxfeRemodelControl: reads the element graph and displacement, computes SED of surface
#    voxels and -- by comparison with lower and upper thresholds -- creates a new model
#    file.
#
# 5. Optional. convertVoxFEScript (again), this time to create a vtk file sequence for
#    reviewing in ParaView
#
# At present, we copy the results at each step into the subdirectory step. 
# The 'StartingDirectory' contains the initial files to be seed the loop; other
# directories will be created as needed

# A quick conversion script, while the new solver beds in
CreatePetScScript() {

  # clear out any old file 
  sed -n 's/VOXEL_SIZE/SET_VOXEL_SIZE/p'                $1  > $1.petsc    ;
  sed -n 's/LOAD_MATERIALS_FILE/LOAD_MATERIALS/p'       $1 >> $1.petsc    ;
  echo ""                                                  >> $1.petsc    ;
  sed -n 's/ALGORITHM_FEA\ JAPCG/SET_ALGORITHM_FEA\ KSPCG\ PCJACOBI/p' $1 >> $1.petsc    ;
  sed -n 's/MAX_ITER/SET_MAX_ITER/p'                    $1 >> $1.petsc    ;
  sed -n 's/TOLERANCE/SET_TOLERANCE/p'                  $1 >> $1.petsc    ;
  echo ""                                                  >> $1.petsc    ;
  sed -n 's/LOAD_MCTSCAN/LOAD_MODEL/p'                  $1 >> $1.petsc    ;
  sed -n 's/[ \t]*SELECT_NODE_FILE/LOAD_CONSTRAINTS/p'  $1 >> $1.petsc    ;
  echo ""                                                  >> $1.petsc    ; 
  sed -n 's/PRINT_X/PRINT_DISPLACEMENTS/p'              $1 >> $1.petsc    ;
}

# ==================================================================================

USAGE="./Remodel-linux.sh  StartDirectory  FinishDirectory  AdaptiveThreshold<Y/N> LowerThreshold UpperThreshold  PetScSolver<Y/N> BaseName"
EXAMPLE="./Remodel-linux.sh  1 25  y 0.05 0.95  y  NameWithoutScriptExtn"

if [ $# != 7 ] ; then
    echo $USAGE
    echo "eg. " $EXAMPLE
    echo ""
    exit 1;
else    
    # so we just might remember how this was called, once we've had several attempts
    echo $@ > Remodel-linux.args.txt
fi

#example calls....
# mpirun -n 8 /home/richard/src/async/vox-fe/solver/PARA_BMU z-test.script
# /home/richard/src/async/vox-fe/tools/VoxFETools/itkReader/convertVoxFEScript z-test.script M
# /home/richard/src/metis-5.1.0/build/Linux-x86_64/programs/m2gmetis -gtype=dual z-test.script.metis z-test.script.graph
# voxfeRemodelControl z-test.script z-test.script.elmap z-test.graph 0.1 1 displacement.txt

# 4 routines
#CCPforge
#PARABMU="mpirun -n 8 /home/richard/src/async/vox-fe/solver/PARA_BMU "
#PETSCSOLVER="mpirun -n 8 /home/richard/src/async/new_solver/src/VoxFESolver "
#CONVERTSCRIPT="/home/richard/src/async/vox-fe/tools/VoxFETools/itkReader/convertVoxFEScript "
#M2GMETIS="/home/richard/src/metis-5.1.0/build/Linux-x86_64/programs/m2gmetis -gtype=dual "
#REMODEL="/home/richard/src/async/vox-fe/tools/VoxFETools/remodelGraph/voxfeRemodelControl "

#Sourceforge 
PARABMU="mpirun -n 8 /home/richard/src/async/vox-fe/solver/PARA_BMU "
PETSCSOLVER="mpirun -n 8 /home/richard/src/async/new_solver/src/VoxFESolver "
CONVERTSCRIPT="/home/richard/src/async/vox-fe/sf/gui/VoxFETools/itkReader/convertVoxFEScript "
M2GMETIS="/home/richard/src/metis-5.1.0/build/Linux-x86_64/programs/m2gmetis -gtype=dual "
REMODEL="/home/richard/src/async/vox-fe/sf/tools/Remodelling/RemodelGraphControl/voxfeRemodelControl "


# Files/directories
VOXFESCRIPT=".script"
VOXFEMODEL=".model"
SCRIPTFILE=$7${VOXFESCRIPT}
MODELFILE=$7${VOXFEMODEL}

if [ ! -f $1/${SCRIPTFILE} ]
then
  echo "Cannot find opening script file: $1/${SCRIPTFILE}"
  exit
elif [ ! -f $1/${MODELFILE} ]
then
  echo "Cannot find opening model file: $1/${MODELFILE}"
  exit
fi
  
CONSTRAINTSFILE=`sed -n 's/[ \t]*SELECT_NODE_FILE\ //p' $1/${SCRIPTFILE}`
DISPFILE=`sed -n 's/PRINT_X\ //p' $1/${SCRIPTFILE}`
MATFILE=`sed -n 's/LOAD_MATERIALS_FILE\ //p' $1/${SCRIPTFILE}`
echo "Got files: " ${CONSTRAINTSFILE} ${DISPFILE} ${MATFILE}

# if we want to create a vtk sequence
VTK="vtk"
mkdir ${VTK}

######################### GET TO WORK ##########################

echo "Changing to dir: $1"
cd $1

# Start the loop
j=$1
j=$((j + 1))
for i in `seq $j $2`; do {

 echo "================================ SOLVE ==========================================="
 
	# Choose the solver
  if [ ! -f ${SCRIPTFILE} ]
  then
    echo "========================================================================================="
    echo "*** The file: ${SCRIPTFILE} cannot be found ***"
    echo "========================================================================================="
    exit 2
  fi

  case "$6" in
    [yY]) 
      # Create the alternate new solver script
      CreatePetScScript ${SCRIPTFILE} 
		  PERFORM_SOLVE="${PETSCSOLVER} ${SCRIPTFILE}.petsc"
		  PERFORM_REMODEL="${REMODEL} ${SCRIPTFILE} ${SCRIPTFILE}.elmap ${SCRIPTFILE}.graph $3 $4 $5 ${DISPFILE} Y 1"
		  ;;
    [nN]) 
      PERFORM_SOLVE="${PARABMU} ${SCRIPTFILE}"
 		  PERFORM_REMODEL="${REMODEL} ${SCRIPTFILE} ${SCRIPTFILE}.elmap ${SCRIPTFILE}.graph $3 $4 $5 ${DISPFILE} N 1"
      ;;
    * ) echo "wrong input" && exit;;
  esac

  # 1. do solve to get displacements (if not already done....)
  echo " choosing solver =>  ${PERFORM_SOLVE}"  
  if [[ ! -e FECompleted.txt ]]
  then
    echo "Starting solve...."
    ${PERFORM_SOLVE}
    echo ".....Solve done"    
  fi

 echo "================================ GRAPH ==========================================="
 
  #2. give the node connectivity of each element
  echo " ${CONVERTSCRIPT} ${SCRIPTFILE} M"
  ${CONVERTSCRIPT} ${SCRIPTFILE} M
  
  #3. generate the element graph
  echo "${M2GMETIS} ${SCRIPTFILE}.metis ${SCRIPTFILE}.graph"
  ${M2GMETIS} ${SCRIPTFILE}.metis ${SCRIPTFILE}.graph

  echo "================================ REMODEL ========================================"
    
  #4. Remodel
  echo "${PERFORM_REMODEL}"
  ${PERFORM_REMODEL}

  #5. make a vtk version of the existing model to help create a ParaView sequence/animation
  k=$((i - 1))
  echo "Creating VTK file: ../${VTK}/${SCRIPTFILE}.vtk.$k"
  if [ ! -f ../${VTK}/${SCRIPTFILE}.vtk.$k ]
  then 
    ${CONVERTSCRIPT} ${SCRIPTFILE} V
    mv ${SCRIPTFILE}.vtk ../${VTK}/${SCRIPTFILE}.vtk.$k
  fi

  echo "================================= COPY TO NEXT =================================="
  echo "                                    " $i
  echo "================================================================================="  

  # create the next directory in the chain, to store the new remodelled data
  # (we're still working in the previous directory)
  mkdir ../$i;

  # move the remodelled data (script/model) into the new directory (leaving log and 
  # displacements; materials and constraints remain the same and are just copied)
  cp -f ${SCRIPTFILE} ../$i/${SCRIPTFILE}
  mv -f ${MODELFILE}.001 ../$i/${MODELFILE}
  cp -f ${MATFILE} ${CONSTRAINTSFILE} ../$i

  # I can see why keep the models, but we don't need these again??
  rm -f ${SCRIPTFILE}.001 ${SCRIPTFILE}.elmap ${SCRIPTFILE}.graph ${SCRIPTFILE}.elmap.txt ${SCRIPTFILE}.metis
  
  #change things back in the working directory to start the next pass
  #cp ${MODELFILE} ${MODELFILE}.orig
  #cp ${MODELFILE}.001 ${MODELFILE}
  
  # change to the newly created data directory
  cd ../$i;

} done;

