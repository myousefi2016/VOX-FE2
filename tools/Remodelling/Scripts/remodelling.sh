#!/bin/bash

USAGE="./remodelling.sh  StartDirectory  FinishDirectory  LowerThreshold UpperThreshold BaseName"
EXAMPLE="./remodelling.sh  1 25  0.1 1.0  NameWithoutScriptExtn"

if [ $# != 5 ] ; then
    echo $USAGE
    echo "eg. " $EXAMPLE
    echo ""
    exit 1;
fi


module add metis

CONVERTSCRIPT="/work/ecse0115/ecse0115/richardh/convertVoxFEScript "
M2GMETIS="/work/y07/y07/cse/metis/5.1.0/bin/m2gmetis -gtype=dual "
REMODEL="/work/ecse0115/ecse0115/richardh/voxfeRemodelControl "
# Y==26, N==6 face connectivity
CONNECTIVITY="Y"  
BASENAME=$5
SCRIPTFILE=${BASENAME}.script
MODELFILE=${BASENAME}.model

CONSTRAINTSFILE=`sed -n 's/[ \t]*SELECT_NODE_FILE\ //p' $1/${SCRIPTFILE}`
MATFILE=`sed -n 's/LOAD_MATERIALS_FILE\ //p' $1/${SCRIPTFILE}`
DISPFILE=`sed -n 's/PRINT_X\ //p' $1/${SCRIPTFILE}`
echo "@@@ Got files: " ${CONSTRAINTSFILE} ${DISPFILE} ${MATFILE} " @@@"

# Create vtk files to allow ParaView to load the (re-)models as a sequence
VTK="vtk"
mkdir ${VTK}

# start remodelling
echo "@@@ Starting remodelling process @@@"

for (( i=$1; i<=$2; i++ ))
do {

  cd $i

  echo "Starting remodel at: " $i >> ../timing.txt
  date -I'ns'  >> ../timing.txt

  echo "@@@ Entering directory: $i, and waiting for solver to finish reading input files... @@@"
  while [[ ! -e InputCompleted.txt ]]
  do
    sleep 10
  done
  
  # now we can read the input files and build the initial graph structures
  # Get the node connectivity of each element  
  if [ ! -e ${SCRIPTFILE}.metis ]
  then
    echo " ${CONVERTSCRIPT} ${SCRIPTFILE} M"
    ${CONVERTSCRIPT} ${SCRIPTFILE} M
  fi
  
  # Generate the element graph
  if [ ! -e ${SCRIPTFILE}.graph ]
  then
    echo "${M2GMETIS} ${SCRIPTFILE}.metis ${SCRIPTFILE}.graph"
    ${M2GMETIS} ${SCRIPTFILE}.metis ${SCRIPTFILE}.graph
  fi

  # Make a vtk version of the existing model to help create a ParaView sequence/animation
  if [ ! -e  ../${VTK}/${BASENAME}.vtk.$i ]
  then
    echo "Creating VTK file: ../${VTK}/${SCRIPTFILE}.vtk.$i"
    ${CONVERTSCRIPT} ${SCRIPTFILE} V
    mv ${SCRIPTFILE}.vtk ../${VTK}/${BASENAME}.vtk.$i
  fi

  # Wait for the solver to complete
  while [[ ! -e FECompleted.txt ]]
  do
    sleep 10
  done
  
  # Check the solver finished ok
  output=`grep solver_done "FECompleted.txt"`
  if [ -z "$output" ] 
  then 
    echo "@@@ Remodel detected that solve failed at step: $i @@@"; 
    exit   
  else 
    echo "@@@ Step: $i solved ok @@@"; 
  fi

  # Do remodelling
  echo "${REMODEL} ${SCRIPTFILE} ${SCRIPTFILE}.elmap ${SCRIPTFILE}.graph $3 $4 ${DISPFILE} Y 1"
  ${REMODEL} ${SCRIPTFILE} ${SCRIPTFILE}.elmap ${SCRIPTFILE}.graph $3 $4 ${DISPFILE} Y 1

  # Move key files to the next directory
  j=$i
  j=$((j + 1))
  cp -f ${SCRIPTFILE} ../$j/${SCRIPTFILE}
  mv -f ${MODELFILE}.001 ../$j/${MODELFILE}
  cp -f ${MATFILE} ${CONSTRAINTSFILE} ../$j
  
  #remove the unwanted intermediates
  rm -f ${SCRIPTFILE}.001 ${SCRIPTFILE}.elmap ${SCRIPTFILE}.graph ${SCRIPTFILE}.elmap.txt ${SCRIPTFILE}.metis  
  
  #signal to the solver to start the next step
  echo "@@@ Remodelling completed iteration: $i @@@"
  touch Remodelling_done
  cd ..

} done

echo "Remodelling completed, exiting..."
