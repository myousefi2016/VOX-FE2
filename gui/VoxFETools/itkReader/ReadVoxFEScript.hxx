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

#ifndef READ_VOXFE_SCRIPT_HXX
#define READ_VOXFE_SCRIPT_HXX


#include "ReadVoxFEScript.h"
using namespace std;

template<class T>
bool rvsReadVoxelFile( const char* file,
                    Element<T>*& element, T& n_elements,
		            const ModelMetaData<T>& model,
                    T  minmaxZ[2],
                    map<T, T>& nodes,
                    const int _border[3],
                    T startValue ) {


  ifstream fin( file );
  if( !fin.is_open() ) return false;

  minmaxZ[1] = 0; //max
  minmaxZ[0] = std::numeric_limits<T>::max(); //min

  fin >> skipws >> n_elements;

  if( element ) { delete[] element; element = 0;  }

  if( n_elements > 0 ) {
      element = new (std::nothrow) Element<T>[ n_elements ];

    if( !element ) {
      fin.close(); return false;
    }
  }

  //collect set of unique nodes
  nodes.clear();

#ifndef VOXFE_READ_DEBUG
  T elNum;
#endif

  T zero = 0, a[NODESPERELEMENT]; //voxel file contains index of voxels, we want nodes

  for( T k = 0; k < n_elements; k++ ) {

    fin >>
#ifdef VOXFE_READ_DEBUG
           element[ k ].elNum
#else
           elNum
#endif
        >> element[ k ].elType
        >> element[ k ].index[0] >> element[ k ].index[1] >> element[ k ].index[2];

    //voxfe seems to use 1-based indices, store as zero-base for now but add in the border
    for( T m = 0; m < 3; m++ ) {
      element[ k ].index[m] += (_border[m] - voxfeVoxelOffset);
    }

    //find min/max indices
    if( element[ k ].index[2] < minmaxZ[0] ) minmaxZ[0] = element[ k ].index[2];
    if( element[ k ].index[2] > minmaxZ[1] ) minmaxZ[1] = element[ k ].index[2];

    if( fin.fail() ) {
      fin.close(); return false;
    }

    //collect set of unique nodes, 8 for each voxel
    rvsGetVoxelNodes( model, element[k].index, a );
    for( T m = 0; m < NODESPERELEMENT; m++ )
      nodes.insert( pair<T, T>( a[m], zero ) );  //zero is a dummy here
  }

  //create the reverse mapping from node index (map key) to position in map (map value)
  rvsMapNodes( nodes, startValue );

  fin.close();
  return true;
}


template<class T>
bool rvsWriteVoxelFile( const char* file,
                        const Element<T>*& element,
                        const T& n_elements,
                        const int _border[3] ) {


  ofstream fout( file );
  if( !fout.is_open() ) return false;

  fout << "\n" << n_elements << "\n";

  for( T k = 0; k < n_elements; k++ ) {

    //voxfe seems to use 1-based indices, store as zero-base for now but add in the border
    fout <<
#ifdef VOXFE_READ_DEBUG
           element[ k ].elNum
#else
           k+1 << " "
#endif
        << element[ k ].elType   << " "
        << element[ k ].index[0] + _border[0] + voxfeVoxelOffset << " "
        << element[ k ].index[1] + _border[1] + voxfeVoxelOffset << " "
        << element[ k ].index[2] + _border[2] + voxfeVoxelOffset << "\n";
  }

  fout.close();
  return true;
}




template<class T>
bool rvsCreateNodeMap( const char* file,
                    T& n_elements,
		                const ModelMetaData<T>& model,
                    T  minmaxZ[2],
                    map<T, T>& nodes,
                    T startValue,
                    T* nodeMapToFile,
                    const char* nodeMapFile ) {


  ifstream fin( file );
  if( !fin.is_open() ) return false;

  minmaxZ[1] = 0; //max
  minmaxZ[0] = std::numeric_limits<T>::max(); //min

  fin >> skipws >> n_elements;

  //elements are not retained here, create one for reading in
  Element<T> element;
  if( n_elements <= 0 ) {
    fin.close();
    return false;
  }

  //collect set of unique nodes
  nodes.clear();
  T zero = 0, a[NODESPERELEMENT]; //voxel file contains index of voxels, we want nodes

  if( !nodeMapToFile ) { //store map data in 'nodes'

    for( T k = 0; k < n_elements; k++ ) {

      //voxfe seems to use 1-based indices, but store as zero-base for now
#ifdef VOXFE_READ_DEBUG
  fin >> element.elNum
      >> element.elType
      >> element.index[0] >> element.index[1] >> element.index[2];
#else
  T elNum;
  fin >> elNum
      >> element.elType
      >> element.index[0] >> element.index[1] >> element.index[2];
#endif

      for( T m = 0; m < 3; m++ ) {
        //element.index[m] -= 1;
        element[ k ].index[m] += (model.borderOffset[m] - 1);
      }

      //find min/max indices
      if( element.index[2] < minmaxZ[0] ) minmaxZ[0] = element.index[2];
      if( element.index[2] > minmaxZ[1] ) minmaxZ[1] = element.index[2];

      if( fin.fail() ) {
        fin.close();
        return false;
      }

      //collect set of unique nodes, 8 for each voxel
      rvsGetVoxelNodes( model, element.index, a );

      for( T m = 0; m < NODESPERELEMENT; m++ )
        nodes.insert( pair<T, T>( a[m], zero ) );  //zero is a dummy here
    }

    //create the reverse mapping from node index (map key) to position in map (map value)
    rvsMapNodes( nodes, startValue );

  }
  else {

    fstream* mapFile = NULL;  //pointer used here for scope
    set<T> *topSet, *botSet;  //pointers here as convenient to swap
    T currentLayer = 0, nodeCount = 0;

    if( !nodeMapFile ) {
      cerr << "\nCannot open the arg (Null) nodeMapFile**\n";
      return false;
    }
    mapFile = new (std::nothrow) fstream(nodeMapFile, ios::out | ios::binary);

    if( !mapFile || !mapFile->is_open() ){
      cerr << "\nCannot open the arg (" << nodeMapFile << ")**\n";
      return false;
    }

    //populate the map file with -1s
    T dummy = -1; //fixme: what if an unsigned type is used ??
    for( T k = 0; k < model.XYZ; k++ ) {
      mapFile->write((char*)&dummy, sizeof(T) );
    }
    mapFile->close(); //need to close, so we can, open for reading+writing
    mapFile->open(nodeMapFile, ios::in | ios::out | ios::binary);

    topSet = new (std::nothrow) set<T>;
    botSet = new (std::nothrow) set<T>;

    //Loop through elements: here we have to take more care as nodes are written out
    //as we complete each layer of voxels. So
    //  (1) We read a voxel
    //  (2) if it's the first voxel, set the first layer as current.
    //  (3) If the next layer has been reached output the bottom layer
    //      of nodes and swap the top and bottom pointers to merge the node
    //      layers appropriately.
    //  (4) At the last voxel we repeat (3) if at the next layer and the output the last top
    //      layer of nodes.

    //Process elements
    for( T k = 0; k < n_elements; k++ ) {

      //voxfe seems to use 1-based indices, but store as zero-base for now
      if( !rvsCreateNodeMap_ReadElement( fin, element, model, minmaxZ, a ) ) {
        cerr << "Failed to read element\n";
        return false;
      }
#ifdef VOXFE_READ_DEBUG
      cerr << "k: " << k << "\n";
#endif
      if( k == 0 ) currentLayer = element.index[2];  //fixme: maybe we could 'peek' instead of doing this
                                                     //a billion times

      if( element.index[2] > currentLayer ) {

        rvsCreateNodeMap_SwapLayers( mapFile,
                                     topSet, botSet,
                                     currentLayer, nodeCount );
      }

      //check if the last element---------------------------------------------------
      T pos, entry;
      if( k == (n_elements-1) ) {

        //"Add-In"
        for( int m = 0; m < HALFNODESPERELEMENT; m++ ) {
          botSet->insert( a[m] );
          topSet->insert( a[m+HALFNODESPERELEMENT] );
        }

        //need to output bottom layer after above add-in
        //cerr << "last bot set size: " << botSet->size() << "\n";
        typename set<T>::const_iterator sit = botSet->begin();
        for( ; sit != botSet->end(); ++sit ) {

          pos = (*sit) * sizeof(T);
          if( rvsCreateNodeMap_WriteMapEntry( mapFile, pos, nodeCount ) )
            nodeCount++;
        }

        //output data for the last top layer
        //cerr << "last top set size: " << topSet->size() << "\n";
        sit = topSet->begin();
        for( ; sit != topSet->end(); ++sit ) {

          T pos = (*sit) * sizeof(T);
          if( rvsCreateNodeMap_WriteMapEntry( mapFile, pos, nodeCount ) )
            nodeCount++;
        }
      }
      else {    //if not last element

        //"Add-In"
        for( int m = 0; m < HALFNODESPERELEMENT; m++ ) {
          botSet->insert( a[m] );
          topSet->insert( a[m+HALFNODESPERELEMENT] );
        }
      }

    }  //for all elements

    //return the node count
    dummy = rvsCreateNodeMap_RemapNodes( mapFile, model.XYZ );
    if( nodeCount != dummy ) {
      cerr << "**** Recount failed ****\n";
    }
    *nodeMapToFile = nodeCount;

    //wrap up
    if( mapFile ) { mapFile->close(); delete mapFile; }
    if( topSet )  { topSet->clear();  delete topSet; }
    if( botSet )  { botSet->clear();  delete botSet; }

  }  // map to file

  fin.close();

  return true;
}


template<class T>
bool rvsCreateNodeMap_ReadElement( std::ifstream& fin,
                                   Element<T>& element,
                                   const ModelMetaData<T>& model,
                                   T minmaxZ[2],
                                   T gnode[NODESPERELEMENT] ) {
#ifdef VOXFE_READ_DEBUG
  fin >> element.elNum
      >> element.elType
      >> element.index[0] >> element.index[1] >> element.index[2];
#else
  T elNum;
  fin >> elNum
      >> element.elType
      >> element.index[0] >> element.index[1] >> element.index[2];
#endif

  for( T m = 0; m < 3; m++ ) {
    //element.index[m] -= 1;
    element.index[m] += (model.borderOffset[m] - 1);
  }

  //find min/max indices
  if( element.index[2] < minmaxZ[0] ) minmaxZ[0] = element.index[2];
  if( element.index[2] > minmaxZ[1] ) minmaxZ[1] = element.index[2];

  if( fin.fail() ) {
    fin.close();
    return false;
  }

  rvsGetVoxelNodes( model, element.index, gnode );     //collect set of unique nodes, 8 for each voxel

  return true;
}


template<class T>
void rvsCreateNodeMap_SwapLayers( std::fstream* mapFile,
                                  set<T>*& topSet,
                                  set<T>*& botSet,
                                  T& currentLayer,
                                  T& nodeCount ) {

   //if we've moved to the next layer, output the bottom set of nodes
  //and swap the pointers

  //update the map with the bottom set nodes
  typename set<T>::const_iterator sit = botSet->begin();
  for( ; sit != botSet->end(); ++sit ) {

    T pos = (*sit) * sizeof(T);

    if( rvsCreateNodeMap_WriteMapEntry( mapFile, pos, nodeCount ) )
      nodeCount++;
  }

  //the bottom set can now be cleared
  botSet->clear();

  set<T>* tmpSet = topSet;
  topSet = botSet;  //the new top set will be empty
  botSet = tmpSet;  //the old top set will now be merged with bottom nodes from the next layer
  currentLayer++;

  //cerr << "after swap bot set size: " << botSet->size() << "  top set size: " << topSet->size() << "\n";
}


template<class T>
bool rvsCreateNodeMap_WriteMapEntry( std::fstream* bin_out,
                                     const T& location,
                                     const T& new_value ) {

  T entry;
  bin_out->seekg( location );
  bin_out->read( (char*)&entry, sizeof(T) );

  if( (entry < 0) || (new_value < entry) ) {
    bin_out->seekp( location );
    bin_out->write( (char*)&new_value, sizeof(T) );

#ifdef VOXFE_READ_DEBUG
    cerr << "writing node to mapFile: " << location/sizeof(T) << "  (count: " << new_value << ")\n";
#endif
    return true;
  }

  return false;
}


template<class T>
T rvsCreateNodeMap_RemapNodes( std::fstream* bin_out, const T& size ) {

  T entry, recount=0;
  bin_out->seekg( 0 );
  bin_out->seekp( 0 );
  for( T k = 0; k < size; k++ ) {

    bin_out->read( (char*)&entry, sizeof(T) );

    if( entry < 0 ) continue;

    bin_out->seekp( k*sizeof(T) );
    bin_out->write( (char*)&recount, sizeof(T) );
    recount++;
  }

  return recount;
}




template<class T>
T rvsReadVoxelFileInChunks( ifstream& fin, Element<T>*& element,
                            T& n_elements, const unsigned int& chunk ) {

#ifndef VOXFE_READ_DEBUG
  T elNum; //elNum is dummy
#endif

  T k = 0;

  if( !fin.is_open() ) return k;

  fin >> skipws >> n_elements;

  if( element ) { delete[] element; element = 0;  }

  if( n_elements > 0 )
      element = new (std::nothrow) Element<T>[ n_elements ];

  if( !element ) {
      fin.close(); return k;
  }

  for( k = 0; k < chunk; k++ ) {

    fin >>
#ifdef VOXFE_READ_DEBUG
             element[ k ].elNum
#else
             elNum
#endif
        >> element[ k ].elType
        >> element[ k ].index[0] >> element[ k ].index[1] >> element[ k ].index[2];

    //voxfe seems to use 1-based indices, but store as zero-base for now
    for( T m = 0; m < 3; m++ ) element[ k ].index[m] -= 1;

    if( !fin.good() ) break;
  }

  return k;
}


template<class T>
T rvsReadVoxelFileInChunks( ifstream& fin, Element<T>*& element, const unsigned int& chunk ) {

#ifndef VOXFE_READ_DEBUG
  T elNum; //elNum is dummy
#endif

  T k = 0;
  for( k = 0; k < chunk; k++ ) {

    fin >>
#ifdef VOXFE_READ_DEBUG
             element[ k ].elNum
#else
             elNum
#endif
        >> element[ k ].elType
        >> element[ k ].index[0] >> element[ k ].index[1] >> element[ k ].index[2];

    //voxfe seems to use 1-based indices, but store as zero-base for now
    for( T m = 0; m < 3; m++ ) element[ k ].index[m] -= 1;

    if( !fin.good() ) break;
  }

  return k;
}


template<class T>
bool rvsReadModelData( const char* file, const int _border[3], ModelMetaData<T>& model ) {

  enum ConstrainFileStatus {
    CONSTRAINT_FILE_TO_BE_CREATED,  //Is expected to be created, eg using ParaView
    CONSTRAINT_FILE_EXISTS,         //Has been created
    NO_SEPARATE_CONSTRAINT_FILE     //Constraints are listed in the script file (old)
  };


  string dummy;
  MaterialType material;
  int materialNum[2];

  ifstream fin( file );
  if( !fin.is_open() ) {
    return false;
  }

  //read material data, at least first 2 lines
  fin >> dummy;
  cout << dummy << endl;  //test

  while( dummy.compare("VOXEL_SIZE") ) {

    string mat_file;
    if( !dummy.compare("LOAD_MATERIALS_FILE") ) {
      fin >> mat_file;
      cout << "found material file: " << mat_file << endl;

      model.ReadMaterialsFile( mat_file.c_str() );
    }
    else {
      if( !dummy.compare("YOUNG_MODULUS") ) fin >> materialNum[0] >> material.YoungsModulus;
      cout << dummy << endl;  //test

      fin >> dummy;
      if( !dummy.compare("POISSON_RATIO") ) fin >> materialNum[1] >> material.PoissonsRatio;
      cout << dummy << endl;  //test

      if( materialNum[0] == materialNum[1] ) {

        model.materialMap.insert(  std::pair<int, MaterialType>(materialNum[1],material) );

      } else {
        cout << "prob at material stage" << endl;  //test
        fin.close(); return false;
      }
    }

    fin >> dummy;
  }

  if( fin.fail() ) {
    fin.close(); return false;
  }

  //at this point we're at the VOXEL_SIZE line
  fin >> model.voxSize[0] >> model.voxSize[1] >> model.voxSize[2] >> model.voxSize[3];

  fin >> skipws >> dummy;
  cout << dummy << endl;  //test
  if( !dummy.compare("ALGORITHM_FEA") ) {

    fin >> model.algorithm;

    //for PetSc code, we want to be able to add a sub-type
    if( !model.algorithm.compare("KSPCG") ) {
      fin >> dummy;
      model.algorithm += dummy;
    }
  }
  else {
    fin.close(); return false;
  }

  fin >> skipws >> dummy;
  cout << dummy << endl;  //test
  if( !dummy.compare("MAX_ITER") ) fin >> model.maxIterations;
  else {
    fin.close(); return false;
  }

  fin >> skipws >> dummy;
  cout << dummy << endl;  //test
  if( !dummy.compare("TOLERANCE") ) fin >> model.tolerance;
  else {
    fin.close(); return false;
  }

  fin >> skipws >> dummy;
  cout << dummy << endl;  //test
  if( !dummy.compare("LOAD_MCTSCAN") ) {
    model.scanType = dummy;
    fin >> model.scanDim[0] >> model.scanDim[1] >> model.scanDim[2];
    fin >> model.voxelFile;
  }
  else {
    fin.close(); return false;
  }

  //===============================================================================
  //Read rest of constraint file NB - no error checking here, returns true whatever
  //===============================================================================

  model.scanDim[0] += (2 * _border[0]);
  model.scanDim[1] += (2 * _border[1]);
  model.scanDim[2] += (2 * _border[2]);

  model.X   = model.scanDim[0] + 1;
  model.XY  = model.X  * (model.scanDim[1] + 1);
  model.XYZ = model.XY * (model.scanDim[2] + 1);

  //If border was read, need to set up the next entry
  fin >> skipws >> dummy;
  cout << dummy << endl;


  //===============================================================================
  //===============================================================================

  if( !dummy.compare("SELECTION_OF_NODES") ) {

    fin >> skipws >> dummy;

    ifstream fstr;

    //determine the status of constraints -- are they in a separate file,
    //if so determine if the file is available (it may not be)
    ConstrainFileStatus constraint_file_status;
    if( !dummy.compare("SELECT_NODE_3D") || !dummy.compare("PRESERVE_NODE_3D") ) {
      constraint_file_status = NO_SEPARATE_CONSTRAINT_FILE;
    }
    else if( !dummy.compare("SELECT_NODE_FILE") ) {
      fin >> skipws >> dummy;
      fstr.open( dummy.c_str() );
      if( !fstr.is_open() ) {
        constraint_file_status = CONSTRAINT_FILE_TO_BE_CREATED;
        cerr << "\n ** Cannot open constraints file ** \n\n";
        fstr.close();
      }
      else constraint_file_status = CONSTRAINT_FILE_EXISTS;

      fstr >> skipws >> dummy;
    }

    if( constraint_file_status == NO_SEPARATE_CONSTRAINT_FILE ) {

      while( !fin.eof() && dummy.compare("COMPUTE_SED") ) {

        if( !dummy.compare("SELECT_NODE_3D") ) {
          model.AddConstraint( fin );
        }
        else if( !dummy.compare("PRESERVE_NODE_3D") ) {
          model.AddConstraint( fin, true );
        }
        fin >> skipws >> dummy;
      }

      if( !dummy.compare("COMPUTE_SED") )
        cout << "Finished reading constraints\n";

    }
    else if( constraint_file_status == CONSTRAINT_FILE_EXISTS ) {

      while( !fstr.eof() ) {

        if( !dummy.compare("SELECT_NODE_3D") ) {
          model.AddConstraint( fstr );
        }
        else if( !dummy.compare("PRESERVE_NODE_3D") ) {
          model.AddConstraint( fstr, true );
        }
        fstr >> skipws >> dummy;
      }
      fstr.close();
    }

    fin >> skipws >> dummy;  //should be the COMPUTE_SED line
    if( !dummy.compare("COMPUTE_SED") )
      cout << "Finished reading constraints\n";
  }

  fin.close();
  return true;
}


template<class T>
void rvsGetVoxelNodes( const ModelMetaData<T>& model, const T index[3],
                    T a[NODESPERELEMENT] )  {

  //as first defined in Global assembly
  a[0] = (index[2] * model.XY) + (index[1] * model.X) + index[0]; //(el_k * XY) + (el_j * NodesX) + el_i;
  a[1] = a[0] + 1;
  a[2] = a[0] + model.X;
  a[3] = a[2] + 1;
  a[4] = a[0] + model.XY;
  a[5] = a[1] + model.XY;
  a[6] = a[2] + model.XY;
  a[7] = a[3] + model.XY;
}


template<class T>
void rvsMapNodes( map<T, T>& nodes, T startValue ) {

  T count = startValue;  //eg. Ansys uses 1-based indexing

  typename map<T,T>::iterator mit = nodes.begin();
  for(  ;  mit != nodes.end(); mit++ ) {

    //record the count as value in a reverse mapping
    mit->second = count;
    count++;
  }
}


template<class T>
void rvsPrintNodes( const ModelMetaData<T>& model, const map<T, T>& nodes,
                 ofstream& fout ) {

  fout << "# Nodes" << endl;
  T i[3], rem;
  typename map<T,T>::const_iterator mit = nodes.begin();
  for(  ;  mit != nodes.end(); mit++ ) {

    model.GetIndexXYZ( mit->first, i );

    fout << mit->first << "  :  "
         << i[0] << "  "
	       << i[1] << "  "
	       << i[2] << endl;
  }
}


template<class T>
bool rvsReadLSMs( const char* lsmFile, const int& lsmSize, ModelMetaData<T>& model ) {

  std::ifstream fin( lsmFile );
  if( !fin.is_open() )
    return false;

  string str;
  char cbuff[LSM_BUFF_LENGTH];
  map< int, MaterialType >::iterator it = model.materialMap.begin();
  for( ; it != model.materialMap.end(); ++it ) {

    cout << "looking for material: " << it->first << "\n";

    //loop through the lsm file looking for it->first
    long int material, msize, p;
    bool found = false;
    while ( !found ) {
      std::getline ( fin, str );
      if( !fin.good() ) break;
      //cout << "got line: " << str << endl;

      if( str[0] = '#' ) { //look for the first space

        const char* cstr = str.c_str();
        cstr++;
        material = strtol( cstr, 0, 10 );

        if( material == it->first ) {
          cout << "found material: " << material << " from: " << cstr << "\n";
          found = true;

          //read in rest of LSM
          msize = lsmSize * lsmSize;
          it->second.LSM = new (std::nothrow) double[ msize ];
          if( !it->second.LSM ) return false;
          for( p=0; p<msize; ++p ) fin >> it->second.LSM[p];
        }
        else it->second.LSM = 0;  //should not occur
      }
      cbuff[0] = '\0';  //appends otherwise
      str.clear();      //ditto
    }
  }

  fin.close();
  return true;

}

template<class T>
void rvsFindMinMax( const Element<T>* element, const T& n_elements,
                    T min_index[3], T max_index[3], bool initialise ) {

  if( initialise ) {

    //if not done, assume the supplied args are treated as running mins/maxs
    for( int p = 0; p < 3; p++ ) {
      min_index[p] = std::numeric_limits<T>::max();
      max_index[p] = std::numeric_limits<T>::min();
    }
  }

  for( T k = 0; k < n_elements; k++ ) {

    for( int p = 0; p < 3; p++ ) {

      if( element[ k ].index[p] < min_index[p] )
         min_index[p] = element[ k ].index[p];

      if( element[ k ].index[p] > max_index[p] )
         max_index[p] = element[ k ].index[p];
    }
  }

}


ImageGroupData* rvsReadImageGroupFile( const char* filepath ) {

  ifstream fin( filepath );
  if( !fin.is_open() ) return 0;

  ImageGroupData* idata = new ImageGroupData;
  string dummy;

  //Allow any number of comment lines starting with a single '#'
  char dummy_line[VOXFE_COMMENT_LINE_LENGTH];
  while ( fin.peek() == '#' ) {
    fin.getline( dummy_line, VOXFE_COMMENT_LINE_LENGTH );
  }


  fin >> dummy;
  if( !dummy.compare("LOAD_IMAGE") )  fin >> idata->ImageFile;

  fin >> dummy;

  //some options here... allow any order
  while( dummy.compare("GROUPS") ) {  //ie not GROUPS

    //check if voxel size is being supplied
    if( !dummy.compare("VOXEL_SIZE") ) {
      fin >> idata->VoxelSize;
    }

    //BORDER_OFFSET -- made optional here
    if( !dummy.compare("BORDER_OFFSET") ) {
      fin >> idata->BorderOffset[0] >> idata->BorderOffset[1] >> idata->BorderOffset[2];
    }

    fin >> skipws >> dummy;
    //cout << dummy << endl;
  }

  int nGroups = 0;
  if( !dummy.compare("GROUPS") ) fin >> nGroups;

  short startGroup, endGroup, remodelFlag;
  for( int k=0; k<nGroups; k++ ) {

    fin >> startGroup >> endGroup >> dummy >> remodelFlag;
    idata->Group.push_back( ImageGroupData::GroupDef( dummy, startGroup, endGroup, remodelFlag ) );
  }

  fin.close();

  return idata;
}

template<class T>
void rvsGetMaterialSet( const Element<T>* element, const T& n_elements, set<short>& label_set ) {

  label_set.clear();
  for( T k = 0; k < n_elements; k++ ) {

    short sh = element[ k ].elType;
    label_set.insert( sh );
  }
}

template<class T>
bool rvsWriteNodeMapKeys( const char* bin_file,
                      const map<T, T>& nodes ) {

  fstream bin_out( bin_file, ios::out | ios::binary );
  if( !bin_out.is_open() ){
    cerr << "\nCannot open the file (" << bin_file << ")**\n";
    return false;
  }

  T numNodes = nodes.size();
  bin_out.write( (char*)&numNodes, sizeof(T) );

  typename map<T,T>::const_iterator mit = nodes.begin();
  for(  ;  mit != nodes.end(); mit++ ) {

    bin_out.write( (char*)&(mit->first), sizeof(T) );
  }

  bin_out.close();
  return true;
}

template<class T>
bool rvsWriteNodeMapKeysASCII( const char* asc_file,
                      const map<T, T>& nodes ) {

  ofstream asc_out( asc_file );
  if( !asc_out.is_open() ){
    cerr << "\nCannot open the file (" << asc_file << ")**\n";
    return false;
  }

  T numNodes = nodes.size();
  asc_out << numNodes << "\n";

  typename map<T,T>::const_iterator mit = nodes.begin();
  for(  ;  mit != nodes.end(); mit++ ) {

    asc_out << mit->first << "\n";
  }

  asc_out.close();
  return true;
}


template<class T>
T* rvsReadNodeMapKeys( const char* bin_file, T& numNodes ) {

  numNodes = 0;
  fstream bin_in( bin_file, ios::in | ios::binary );
  if( !bin_in.is_open() ){
    cerr << "\nCannot open the file (" << bin_file << ")**\n";
    return 0;
  }

  bin_in.read( (char*)&numNodes, sizeof(T) );

  T* array = 0;
  if( numNodes > 0 ) array = new (std::nothrow) T[ numNodes ];
  if( !array ) {
    cerr << "\nCannot allocate nodemap...\n";
    return 0;
  }

  for( T k=0; k<numNodes; k++ ) {

    bin_in.read( (char*)&(array[k]), sizeof(T) );
  }

  bin_in.close();
  return array;
}

//===========================================

template<class T>
bool rvsWriteElementIDs( const char* bin_file, const Element<T>* element, const T& n_elements, const ModelMetaData<T>& model ) {

  fstream bin_out( bin_file, ios::out | ios::binary );
  if( !bin_out.is_open() ){
    cerr << "\nCannot open the file (" << bin_file << ")**\n";
    return false;
  }

  T numEls = n_elements;
  bin_out.write( (char*)&numEls, sizeof(T) );

  //create a set to order the indices, ... may have been read in in any order
  map<T,unsigned char> indices;
  for( T k = 0; k < n_elements; k++ ) {

    T index = element[k].index[2]*model.XY + element[k].index[1]*model.X + element[k].index[0];
    indices.insert(pair<T,unsigned char>(index, static_cast<unsigned char>(element[k].elType) ) );
  }

  typename map<T,unsigned char>::iterator mit = indices.begin();
  for( ; mit != indices.end(); ++mit ) {

    bin_out.write( (char*)&(mit->first),  sizeof(T) );
    bin_out.write( (char*)&(mit->second), sizeof(unsigned char) );
  }

  bin_out.close();
  return true;
}


template<class T>
bool rvsWriteElementIDsASCII( const char* asc_file, const Element<T>* element, const T& n_elements, const ModelMetaData<T>& model ) {

  ofstream asc_out( asc_file );
  if( !asc_out.is_open() ){
    cerr << "\nCannot open the file (" << asc_file << ")**\n";
    return false;
  }

  asc_out << n_elements << "\n";

  //create a set to order the indices, ... may have been read in in any order
  map<T,unsigned short> indices;
  for( T k = 0; k < n_elements; k++ ) {

    T index = element[k].index[2]*model.XY + element[k].index[1]*model.X + element[k].index[0];
    indices.insert(pair<T,unsigned short>(index, element[k].elType));
  }  //for all elements

  typename map<T,unsigned short>::iterator mit = indices.begin();
  for( ; mit != indices.end(); ++mit ) {

    asc_out << mit->first << " " << mit->second << "\n";
  }

  asc_out.close();
  return true;
}

template<class T>
void rvsReadElementIDs( const char* bin_file, T& numNodes, T*& IDs, unsigned char*& materials ) {

  numNodes = 0;
  fstream bin_in( bin_file, ios::in | ios::binary );
  if( !bin_in.is_open() ){
    cerr << "\nCannot open the file (" << bin_file << ")**\n";
    return;
  }

  bin_in.read( (char*)&numNodes, sizeof(T) );

  IDs = 0; materials = 0;
  if( numNodes > 0 ) {
    IDs       = new (std::nothrow) T[ numNodes ];
    materials = new (std::nothrow) unsigned char[ numNodes ];
  }

  if( !IDs ) {
    cerr << "\nCannot allocate nodemap array...\n";
    return;
  }

  if( !materials ) {
    cerr << "\nCannot allocate materials array...\n";
    delete[] IDs; //free this anyway..
    IDs = 0;
    return;
  }


  for( T k=0; k<numNodes; k++ ) {

    bin_in.read( (char*)&(IDs[k]),       sizeof(T) );
    bin_in.read( (char*)&(materials[k]), sizeof(unsigned char) );
  }

  bin_in.close();
}



//create specializations for some functions, so they will exist ie can be called in later
//included files
template rvsIntegerType* rvsReadNodeMapKeys( const char* bin_file, rvsIntegerType& numNodes );


#endif

