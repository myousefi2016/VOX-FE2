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

#ifndef CONVERT_TO_ANSYS_SCRIPT_HXX
#define CONVERT_TO_ANSYS_SCRIPT_HXX


#include "convertVoxFEScript.h"
#include "ReadVoxFEScript.hxx"


template<class T>
bool WriteAnsysSolid45File( const char* opfilename,
			    const map<T, T>& nodes,
          const Element<T>* element,
			    const T& n_elements,
			    const ModelMetaData<T>& model ) {

  std::ofstream fout( opfilename );
  if( !fout.is_open() ) return false;

  //output the header lines
  for( int k=0; k<ANSYSHEADER_LINES; k++ ) fout << ANSYSHEADER[k] << endl;
  fout <<  endl;

  //output materials as eg.
  /*
    "mp,EX,1,17e3	   	! defines Youngs modulus of material type 1",
    "mp,NUXY,1,0.3      ! defines Poissons ratio of material type 1"
  */
  map< int, MaterialType >::const_iterator mmapIt = model.materialMap.begin();
  for( ; mmapIt != model.materialMap.end(); ++mmapIt ) {
    fout << "mp,EX," << mmapIt->first << "," << mmapIt->second.YoungsModulus
	     << "\t! defines Young\'s modulus of this material type" << endl;
    fout << "mp,NUXY," << mmapIt->first << "," << mmapIt->second.PoissonsRatio
	     << "\t! defines Poisson\'s ratio of this material type" << endl;
  }
  fout << endl;

  //find min index to reduce to min vertex pos
  T min_index[3], max_index[3];
  rvsFindMinMax( element, n_elements, min_index, max_index, true );


#ifdef DEBUG
  cout << "Min indices:\n";
  for( int j=0; j<3; j++ ) cout << min_index[j] << "  ";
  cout << "\nVox size:\n";
  for( int j=0; j<3; j++ ) cout << model.voxSize[j] << "  ";
  cout << "\nOffsets:\n";
  cout << ANSYS_OFFSET_X << "  " << ANSYS_OFFSET_Y << "  " << ANSYS_OFFSET_Z << "\n";
  cout << "\nIndices:\n";
#endif

  T i[3], rem;
  typename map<T,T>::const_iterator mit = nodes.begin();
  for(  ;  mit != nodes.end(); mit++ ) {

    //get node indices from unique node index
    i[2] = (mit->first) / model.XY; // "z"
    rem  = (mit->first) % model.XY;
    i[1] = rem / model.X;     // "y"
    i[0] = rem % model.X;     // "x"

    //adjust to min index
    fout << "n," << mit->second << ","
         << (i[0] - min_index[0]) * model.voxSize[0] - ANSYS_OFFSET_X << ","
	 << (i[1] - min_index[1]) * model.voxSize[1] - ANSYS_OFFSET_Y << ","
	 << (i[2] - min_index[2]) * model.voxSize[2] - ANSYS_OFFSET_Z << endl;

#ifdef DEBUG
    cout << mit->second << "  " << i[0] << "  " << i[1] << "  " << i[2] << "     "
         << (i[0] - min_index[0]) * model.voxSize[0] - ANSYS_OFFSET_X << ","
	 << (i[1] - min_index[1]) * model.voxSize[1] - ANSYS_OFFSET_Y << ","
	 << (i[2] - min_index[2]) * model.voxSize[2] - ANSYS_OFFSET_Z << endl;
#endif
  }

  //elements
  T a[NODESPERELEMENT], seq[NODESPERELEMENT], last_material = 0;

  //Nodes: Anysis order different from VOXFE, use seq to reorder
  seq[0] = 0;
  seq[1] = 2;
  seq[2] = 3;
  seq[3] = 1;
  seq[4] = 4;
  seq[5] = 6;
  seq[6] = 7;
  seq[7] = 5;


  for( T k = 0; k < n_elements; k++ ) {

    //materials -- I'm not sure how ansys really handles multiple materials
    //this is my best guess for now.... assume the mat line is output first,
    //then the element node connections
    if( k == 0 ) {
      fout << endl << "mat," << element[k].elType << endl << endl;
      last_material = element[k].elType;
    }
    else if( element[k].elType != last_material ) {
      fout << endl << "mat," << element[k].elType << endl << endl;
      last_material = element[k].elType;
    }

    rvsGetVoxelNodes( model, element[k].index, a );

    fout << "en," << k + 1;
    for( T m = 0; m < NODESPERELEMENT; m++ ) {

      typename map<T,T>::const_iterator mit = nodes.find( a[ seq[m] ] );
      if( mit != nodes.end() )
        fout << "," << mit->second;
    }
    fout << endl;
  }

  fout.close();
  return true;
}

#undef ANSYS_OFFSET_X
#undef ANSYS_OFFSET_Y
#undef ANSYS_OFFSET_Z


template<class T>
bool WriteFeltMeshFile( const char* opfilename,
                        const map<T, T>& nodes,
                        const Element<T>* element,
                        const T& n_elements,
			                  const ModelMetaData<T>& model ) {

  std::ofstream fout( opfilename );
  if( !fout.is_open() ) return false;

  //======================================================================================
  //output the header lines
  T current_line = 0, count = 0;
  fout << FELTHEADER[current_line++];
  fout << FELTHEADER[current_line++];

  //nodes/elements
  fout << nodes.size() << FELTHEADER[current_line++];
  fout << n_elements   << FELTHEADER[current_line++];

  //======================================================================================
  //for FELT, constraints are given with the nodes (in a very joined up way)
  //we need to process the constraints to find which global node they're
  //related to
  T i[3], rem;
  map<T,T> gnodeToFixedConstraintMap;
  map<T,T> gnodeToForceConstraintMap;
  for( T m = 0; m < model.constraint.size() ; m++ ) {

    //get the global (unique) node index
    i[0] =  (model.constraint[m].nodeIndex[2] * model.XY) +
            (model.constraint[m].nodeIndex[1] * model.X) +
             model.constraint[m].nodeIndex[0];

    if( model.constraint[m].type == 'N' ) {

      gnodeToFixedConstraintMap.insert( pair<T,T>( i[0] , m ) );
    }
    else if( model.constraint[m].type == 'F' ) {

      gnodeToForceConstraintMap.insert( pair<T,T>( i[0] , m ) );
    }

#ifdef CONSTRAINT_DEBUG
    //debug=========================================================
    for(int p=0; p<3; p++)
      cout << model.constraint[m].nodeIndex[p] + 1 << "  ";
    cout << "     ";
    for(int p=0; p<3; p++)
      cout << model.constraint[m].nodeIndex[p] * model.voxSize[0] << "  ";
    cout << "         " << model.constraint[m].type << "     ";

    if( gnodeToFixedConstraintMap.find( i[0] ) != gnodeToFixedConstraintMap.end() )
      cout << "IN-FIXED\n";
    else
    if( gnodeToForceConstraintMap.find( i[0] ) != gnodeToForceConstraintMap.end() )
      cout << "IN-FORCE\n";
    else
      cout << "CONSTRAINT-ERROR\n";
    //debug=========================================================
#endif

  }

  //======================================================================================
  //output nodal data
  fout << FELTHEADER[current_line++];

  //work out padding for force 'naming'
  unsigned int numForces = gnodeToForceConstraintMap.size();
  unsigned int numDigits = 1 + (T)log10( (double)numForces );
  char* spec = new char[ numDigits + 20 ];
  sprintf( spec,"force_%%0%uu", (numDigits+1) );

  cout << "Spec: " << spec << endl;

  char* buff = new char[ numDigits + 20 ];
  count = 1;

  typename map<T,T>::const_iterator mit = nodes.begin();
  typename map<T,T>::const_iterator forceConIt, nodeConIt;
  for(  ;  mit != nodes.end(); mit++ ) {

    //get node indices from unique node index
    i[2] = (mit->first) / model.XY; // "z"
    rem  = (mit->first) % model.XY;
    i[1] = rem / model.X;     // "y"
    i[0] = rem % model.X;     // "x"

    fout << mit->second
         << "    x=" << i[0] * model.voxSize[0]
         << "    y=" << i[1] * model.voxSize[1]
         << "    z=" << i[2] * model.voxSize[2];

    //add node constraint
    nodeConIt  = gnodeToFixedConstraintMap.find( mit->first );
    forceConIt = gnodeToForceConstraintMap.find( mit->first );

    if( nodeConIt != gnodeToFixedConstraintMap.end()   ) {

      //must be some nodal constraint, so work out which
      const T& m = nodeConIt->second;
      int sum = 0;
      for( int axis = 0; axis<3; axis++ ) sum += model.constraint[m].axisFixed[axis];
      if( sum == 1 ) {
	if(      model.constraint[m].axisFixed[0] == 1 )
	  fout << "   constraint=roll_yz";
	else if( model.constraint[m].axisFixed[1] == 1 )
	  fout << "   constraint=roll_zx";
	else if( model.constraint[m].axisFixed[2] == 1 )
	  fout << "   constraint=roll_xy";
      }
      else if( sum == 2 ) {
	if(      model.constraint[m].axisFixed[0] == 0 )
	  fout << "   constraint=roll_x";
	else if( model.constraint[m].axisFixed[1] == 0 )
	  fout << "   constraint=roll_y";
	else if( model.constraint[m].axisFixed[2] == 0 )
	  fout << "   constraint=roll_z";
      }
      else {  //sum == 3
        fout << "   constraint=fixed";
      }
    }
    else fout << "   constraint=free";

    //add any force constraint, or add \n
    //FIXME: could make this a lot neater by maintaining a set of forces, for now list all of them
    forceConIt = gnodeToForceConstraintMap.find( mit->first );
    if( forceConIt != gnodeToForceConstraintMap.end() ) {

      //compile the force name
      sprintf( buff, spec, count );
      fout << "     force=" << buff << endl;
      count++;
    }
    else fout << endl;

  }

  //======================================================================================
  //output elements -- Nodes: FELT order different from VOXFE, use seq to reorder
  fout << FELTHEADER[current_line++];
  T seq[NODESPERELEMENT], a[NODESPERELEMENT];
  seq[0] = 0;
  seq[1] = 1;
  seq[2] = 3;
  seq[3] = 2;
  seq[4] = 4;
  seq[5] = 5;
  seq[6] = 7;
  seq[7] = 6;

  char material_num[8];
  for( T k = 0; k < n_elements; k++ ) {

    rvsGetVoxelNodes( model, element[k].index, a );
    fout << k+1 << " nodes=[";
    for( T m = 0; m < NODESPERELEMENT; m++ ) {

      typename map<T,T>::const_iterator mit = nodes.find( a[ seq[m] ] );
      if( mit != nodes.end() )
        fout << mit->second;
      if( m < (NODESPERELEMENT - 1) ) fout << ",";
    }
    fout << "]   material=material_";  //=bone; was assumed
    
    sprintf( material_num, "%04u", element[k].elType );
    fout << material_num << endl;   //FIXME: Name of materials???
  }

  //======================================================================================
  //material props
  fout << FELTHEADER[current_line++];
  map< int, MaterialType >::const_iterator matprop; // = model.materialMap.begin();
  //fout << "bone nu=" << matprop->second.PoissonsRatio
  //     << " E="      << matprop->second.YoungsModulus << endl;
  for(int p=1; p<VOXFE_MAX_MATERIALS; p++ ) {
  
    matprop = model.materialMap.find( p );
    if( matprop != model.materialMap.end() ) {
    
      sprintf( material_num, "%04d", p );
      fout << "material_" << material_num 
           << " nu=" << matprop->second.PoissonsRatio
           << " E="  << matprop->second.YoungsModulus << endl;
    }
  }

  //======================================================================================
  //forces
  count = 1;
  fout << FELTHEADER[current_line++];
  forceConIt = gnodeToForceConstraintMap.begin();
  for( ; forceConIt != gnodeToForceConstraintMap.end(); ++forceConIt ) {

    //compile the force name
    sprintf( buff, spec, count );
    fout << buff;
    fout << std::fixed;  //std::scientific;  //  //
    fout << "  Fx=" << setprecision(10)  << model.constraint[forceConIt->second].force[0];
    fout << "  Fy=" << setprecision(10)  << model.constraint[forceConIt->second].force[1];
    fout << "  Fz=" << setprecision(10)  << model.constraint[forceConIt->second].force[2] << endl;
    count++;
  }


  //======================================================================================
  //output the last block of constraint names
  for( ; current_line < FELTHEADER_LINES; current_line++ ) {

     fout << FELTHEADER[current_line] << endl;
  }

  delete[] buff;
  delete[] spec;
  fout.close();

  return true;
}


template<class T>
bool WriteMetisMeshFile( const char* opfilename,
			const map<T, T>& nodes,
                         const Element<T>* element,
			const T& n_elements,
			const ModelMetaData<T>& model ) {

  std::ofstream fout( opfilename );
  if( !fout.is_open() ) return false;

  //output number of elements [optionally number of weights, here zero, so ignored]
  fout << n_elements << endl;

  //Metis doesn't care about node locations, only the nodes in each element (order is not important either)
  T a[NODESPERELEMENT];
  for( T k = 0; k < n_elements; k++ ) {

    rvsGetVoxelNodes( model, element[k].index, a );

    for( T m = 0; m < NODESPERELEMENT; m++ ) {

      typename map<T,T>::const_iterator mit = nodes.find( a[m] );
      if( mit != nodes.end() )
        fout << mit->second;
      if( m == (NODESPERELEMENT-1) )
        fout << endl;
      else fout << " ";
    }
    //fout << endl;
  }

  fout.close();
  return true;
}


template<class T>
bool WriteVTKMeshFile( const char* opfilename,
                       const map<T, T>& nodes,
                       const Element<T>* element,
                       const T& n_elements,
                       const ModelMetaData<T>& model ) {

  std::ofstream fout( opfilename );
  if( !fout.is_open() ) return false;

  //======================================================================================
  //output the header lines
  T current_line = 0;
  for( int k=0; k<5; k++ ) fout << VTKHEADER[current_line++];

  //nodes
  fout << nodes.size() << VTKHEADER[current_line++];

  //======================================================================================
  //output nodal data
  T i[3], rem;
  typename map<T,T>::const_iterator mit = nodes.begin();
  for(  ;  mit != nodes.end(); mit++ ) {

    //get node indices from unique node index
    i[2] = (mit->first) / model.XY; // "z"
    rem  = (mit->first) % model.XY;
    i[1] = rem / model.X;     // "y"
    i[0] = rem % model.X;     // "x"

    fout << i[0] * model.voxSize[0]  <<  "  "
         << i[1] * model.voxSize[1]  <<  "  "
         << i[2] * model.voxSize[2]  <<  "\n";

  }

  //======================================================================================
  //output elements -- Nodes: VTK same order as VOXFE, no need to reorder
  fout << "\n" << VTKHEADER[current_line++];
  fout << n_elements << "  " << n_elements*(NODESPERELEMENT+1) << "\n";  //+1 for additional size term, given first
  T a[NODESPERELEMENT];

  for( T k = 0; k < n_elements; k++ ) {

    rvsGetVoxelNodes( model, element[k].index, a );
    fout << NODESPERELEMENT << " ";
    for( T m = 0; m < NODESPERELEMENT; m++ ) {

      typename map<T,T>::const_iterator mit = nodes.find( a[m] );
      if( mit != nodes.end() )
        fout << mit->second - 1 << " ";  //0-based indices in VTK!!!
    }
    fout << "\n";
  }
  fout << "\n";

  //======================================================================================
  //cell types
  fout << VTKHEADER[current_line++];
  fout << " " << n_elements << "\n";
  for( T k = 0; k < n_elements; k++ ) {
    fout << VTK_VOXEL_CELLTYPE << " ";
    if( !((k+1)%10) )
      fout << "\n";
  }
  fout << "\n\n";

  //======================================================================================
  //mimic different material types with attached cell data
  fout << VTKHEADER[current_line++];
  fout << " " << n_elements << "\n";
  fout << VTKHEADER[current_line++];
  fout << "  MaterialType  "  << VTKHEADER[current_line++];
  fout << VTKHEADER[current_line++];

  for( T k = 0; k < n_elements; k++ ) {

    fout << element[ k ].elType << " ";
    if( !((k+1)%10) ) fout << "\n";
  }
  fout << "\n";

  //======================================================================================
  fout.close();

  return true;

}



template<class T>
void CollectTopAndBottomNodes( const T  minmaxZ[2], const map<T, T>& nodes,
                               const Element<T>* element, const T& n_elements,
                               const ModelMetaData<T>& model,
                               set<T>& topNodes, set<T>& topElements,
                               set<T>& botNodes, set<T>& botElements ) {

  topNodes.clear();
  botNodes.clear();
  topElements.clear();
  botElements.clear();

  T minNodeZ = (minmaxZ[0]+1) * model.XY; //min unique node id (we want nodes below this)
  T maxNodeZ = minmaxZ[1] * model.XY;     //max unique node id (we want nodes above this)

  T zNode;
  typename map<T,T>::const_iterator mit = nodes.begin();
  for(  ;  mit != nodes.end(); mit++ ) {

    //get node indices from unique node index
    zNode = (mit->first) / model.XY; // "z"

    if( mit->first >= maxNodeZ ) topNodes.insert( mit->second );
    if( mit->first <  minNodeZ ) botNodes.insert( mit->second );
  }

  for( T k = 0; k < n_elements; k++ ) {

    if( element[k].index[2] == minmaxZ[1] ) topElements.insert( k+1 );  //originally 1-based, so keep for now
    if( element[k].index[2] == minmaxZ[0] ) botElements.insert( k+1 );
  }
}

template<class T>
void PrintNodeSet( ofstream& fout, const set<T>& elements ) {

  typename set<T>::const_iterator sit = elements.begin();
  fout << elements.size() << endl;
  int pp = 1;
  for( ; sit != elements.end() ; sit++ ) {
    fout << "  " << *sit;
    if( (pp%5) == 0 ) fout << endl; //try to keep the output display pretty
    pp++;
  }
  fout << endl;

}

#endif

