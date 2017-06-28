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

#include <stdlib.h>
#include <algorithm>

#include "voxfeRemodelControl.h"
#include "itkReader/ReadVoxFEScript.hxx"
#include "strainReader/voxfeStrainData.h"



#define VOXFE_REMODEL_PRINT_SPEC ".%03d" //".%06d"
#define VOXFE_REMODEL_DEFLT_SED  -99.0
#define VOXFE_REMODEL_INTEGER_BUFF 128
#define VOXFE_SCRIPT_EXTN          ".script"


voxfeRemodelControl::VertexData::VertexData() {
  connectivity = 0;
  status   = 0;
  material = 1;
  remodel  = 1;
  SED      = VOXFE_REMODEL_DEFLT_SED;
  component = 0;  //not allocated
}




voxfeRemodelControl::
voxfeRemodelControl(const char* scriptFile,
                    const char* elementMapFile,
                    const char* metisGraphFile ) :
                      Model(0),
                      NeighbourData(0),
                      Displacement(0),
                      UpdateStep(1),
                      NumVoxelsRemoved(0),
                      NumVoxelsAdded(0),
                      ScriptFile(scriptFile),
                      AdaptiveThresholding(false) {

  this->Graph.clear();
  this->Threshold[0] = this->Threshold[1] = 0.0;

  int border[3]; border[0] = border[1] = border[2] = 0;  //dummy here, no wish to add anything

  Model = new ModelMetaData<voxfeRemodelIntType>();  //don't need to read in the whole model, just the metadata

  this->CreateNeighbourDataBits();

  voxfeRemodelIntType* elementIDList = 0;        // Global ID element list
  unsigned char*       materialList  = 0;

  if( Model && rvsReadModelData( scriptFile, border, *Model ) ) {

    cout << "Read model data....\n";

    //get the global IDs of elements with their materials
    voxfeRemodelIntType numMappedElements=0 ;
    rvsReadElementIDs( elementMapFile, numMappedElements, elementIDList, materialList );

    if( !elementIDList || !materialList ) {
      cerr << "Remodeller cannot allocate Element data:\n";
      if( elementIDList ) delete[] elementIDList;
      if( materialList )  delete[] materialList;
      return;
    }

    //construct the control graph with elemnt data and remodel flags
    if( !CreateGraphFromMetis( metisGraphFile, elementIDList, materialList ) ) {
      cerr << "Remodeller cannot build graph\n";
      if( elementIDList ) delete[] elementIDList;
      if( materialList )  delete[] materialList;
      return;
    }

    //add in extra remodel restrictions from constraints -- elements with
    //boundary conditions or loads upon their nodes cannot remodel
    this->CreateConstraintRestrictions();

    //Zero Connnected components data
    this->ResetComponents();

    cout << "\n Graph control built... \n";
  }
  else {
    cout << "Could not read script file: " << scriptFile << "\n";
  }

  if( elementIDList ) delete[] elementIDList;
  if( materialList  ) delete[] materialList;
}


voxfeRemodelControl::~voxfeRemodelControl() {

  //cout << "\nDeleting the control....\n";

  //free the associated graph data
  GraphTypeIterator git = Graph.begin();
  for( ; git != Graph.end(); ++git ) {
    delete git->second;
  }
  this->Graph.clear();

  this->NodeToDisplacementMap.clear();  //pointers refer to data in Displacement array

  //cout << "\nDeleting data\n";

  if( this->NeighbourData ) delete[] this->NeighbourData;
  if( this->Model )         delete   this->Model;
  if( this->Displacement )  delete[] this->Displacement;
}


bool voxfeRemodelControl::
IsOnSurface( const IDType& vertexID ) {

  GraphTypeIterator it = this->Graph.find( vertexID );
  if( it != this->Graph.end() ) {

    return ( (it->second->connectivity) < NeighbourData[ND_ARRAY_LAST] );
  }
  else {
    cerr << "\n** Element: " << vertexID << " is not in the graph **\n\n";
    return false;
  }
}


bool voxfeRemodelControl::
CreateGraphFromMetis( const char* metisGraphFile,
                 const voxfeRemodelIntType* elementIDList,
                 const unsigned char* materialList ) {

  ifstream fin( metisGraphFile );
  if( !fin.is_open() ) {
    cerr << "Cannot open metis file: " << metisGraphFile << "\n";
    return false;
  }

  IDType nElements, nConnections, elNum = 1, neighID;
  fin >> nElements >> nConnections >> skipws;

  vector<IDType> neighbours;
  string s;

  std::getline( fin, s ); //fixme: for some reason, we have to force the reading
                          //of the rest of the (empty) line....

#ifdef VOXFE_REMODEL_DEBUG
  cout << "String: " << s << " gives: " << nElements << " elements and: " << nConnections << " connections\n";
#endif

  for( IDType k=0; k<nElements; k++ ) {

#ifdef VOXFE_REMODEL_DEBUG
    cout << "Elem: " << elNum << " has: ";
#endif

    //read in the list of neighbours (1-based, so convert to zero-base)
    s.clear();
    neighbours.clear();

    //each row in the metis output stores neighbours, so we have to read the whole row and process
    std::getline( fin, s );
    SplitLine(s, ' ', neighbours);

#ifdef VOXFE_REMODEL_DEBUG
    cout << neighbours.size() <<  " neighbours: ";
    for(int p=0; p<neighbours.size(); p++)
      cout << neighbours[p] << " ";
    cout << "\n";
#endif

    //process the neighbours
    this->CreateConnections( elNum, neighbours, elementIDList, materialList );

    elNum++;
  }

  return true;
}

void voxfeRemodelControl::
SplitLine(const string& s, char delim, vector<IDType>& elems) {

  elems.clear();

  stringstream ss(s);
  string item;
  IDType id;
  while (std::getline(ss, item, delim)) {
    ss >> id;
    elems.push_back(id);
  }

#if VOXFE_REMODEL_DEBUG
  for(int k=0; k<elems.size(); k++) {
    cout << elems[k] << "  ";
  }
  cout << "\n";
#endif
}


void voxfeRemodelControl::
CreateConstraintRestrictions() {

  //for each node/force/etc constraint
  IDType index;
  for( int p=0; p < Model->constraint.size(); ++p ) {

    //create refs to make life easier
    const voxfeRemodelIntType& x = Model->constraint[p].nodeIndex[0];
    const voxfeRemodelIntType& y = Model->constraint[p].nodeIndex[1];
    const voxfeRemodelIntType& z = Model->constraint[p].nodeIndex[2];

    //fixme: We could use the graph here and GetNodeNeighbours, but we need to do
    //lookups anyway, so we may as well process all the same (whether nodes
    //are associated with voxels or not)
    for(int p=-1; p<=0; p++ ) {
      for(int q=-1; q<=0; q++ ) {
        for(int r=-1; r<=0; r++ ) {

          index = (z+r)*this->Model->XY + (y+q)*this->Model->X  + (x+p);
          GraphTypeIterator git = this->Graph.find( index );
          if( git != this->Graph.end() ) {

            //tag the voxel as being out-of-bounds for remodelling
            git->second->remodel = 0;
          }
        }
      }
    }

  }
}


bool voxfeRemodelControl::
SetDisplacements( const char* dispFile, bool withLocation,
                  bool removeDisconnectedElements, bool use26NeighConnectivity ) {

  ifstream fin( dispFile );
  if( !fin.is_open() ) return false;

  if( this->Displacement ) {
    delete[] this->Displacement;
    this->Displacement = 0;
  }

  voxfeStrainData* strainData = new (std::nothrow) voxfeStrainData;
  if( strainData ) {

    if( !withLocation ) {  //eg. ParaBMU output

      //this really needs to be redone every time, as the graph has probably changed
      this->InializeNodeToDisplacementMap();

      unsigned long nodeCount = this->NodeToDisplacementMap.size();
      if( !nodeCount ) {
        cerr << "Cannot determine number of nodes\n";
        return false;
      }

      this->Displacement = strainData->ReadBMUSolverDisplacementData( fin, (unsigned long)nodeCount );
      if( !this->Displacement ) {
        cerr << "Cannot read displacement data\n";
        return false;
      }

      cout << "Read displacement data for: " << nodeCount << " nodes\n";
      this->UpdateNodeToDisplacementMap();
    }
    else {

      this->Displacement = strainData->ReadPETSCSolverDisplacementData( fin, this->NodeToDisplacementMap );
      if( !this->Displacement ) {
        cerr << "Cannot read displacement data\n";
        return false;
      }

      cout << "Read displacement data for: " << this->NodeToDisplacementMap.size() << " nodes\n";
    }

  }
  else return false;

  //we can free this now
  delete strainData;
  strainData = 0;

  //compute SEDs and threshold to get voxels to add/remove
  this->SelectBySED( true );

  //Add and remove (status < 0 ) elements and update the model connectivity
  this->Remodel();

  if( removeDisconnectedElements ) {

    this->SetConnectedComponents( use26NeighConnectivity, true );
  }

  this->WriteModel(removeDisconnectedElements);

  this->UpdateStep += 1;
  return true;
}


void voxfeRemodelControl::
Remodel() {

  //shouldn't get an overlap, but we ought to check
  IDSet intersectionSet, diffSet;
  std::set_intersection( voxels_to_add_onto.begin(), voxels_to_add_onto.end(),
                         voxels_to_remove.begin(), voxels_to_remove.end(),
                         std::inserter( intersectionSet, intersectionSet.begin() ) );

  if( !intersectionSet.empty() ) {

    cerr << "\n\n** Intersection set of voxels is not empty! **\n"
                "     -- Check your thresholds --            **\n\n";

    //the only sane thing would seem to be to remove overlapping voxels from both sets
    std::set_difference(voxels_to_add_onto.begin(), voxels_to_add_onto.end(),
                        intersectionSet.begin(), intersectionSet.end(),
                        std::inserter( diffSet, diffSet.begin() ) );
    voxels_to_add_onto = diffSet; //copy back

    diffSet.clear();
    std::set_difference(voxels_to_remove.begin(), voxels_to_remove.end(),
                        intersectionSet.begin(), intersectionSet.end(),
                        std::inserter( diffSet, diffSet.begin() ) );
    voxels_to_remove = diffSet; //copy back

  }
  intersectionSet.clear();
  diffSet.clear();

  //======================= REMOVE VOXELS =============================================
  //remove voxels: here we set status to negative, rather than deleting from the graph
  // -- the advantage being that we will not be left with disconnected regions (although
  //the solution might be badly conditioned, it should not be singular). All connections
  //to and from removed voxels are severed however so that elements at the end of each
  //cut connection will regard themselves as surface voxels
  IDSet& neighbours     = intersectionSet;  //re-use
  IDSet& non_neighbours = diffSet;
  voxfeRemodelIntType vertexXYZ[3], neighXYZ[3];
  int offset[3];
  GraphTypeIterator gitV, gitN;

  IDSetIterator vit = voxels_to_remove.begin();
  for( ; vit != voxels_to_remove.end(); ++vit ) {

    gitV = this->Graph.find( *vit );
    if( gitV != this->Graph.end() && ( !gitV->second->remodel ) ) continue;

    //find all the neighbours of this voxel, so connections can be cut
    if( GetNeighbours( *vit, neighbours, non_neighbours ) ) {

      this->Model->GetIndexXYZ( voxfeRemodelIntType(*vit), vertexXYZ );

      IDSetIterator nit = neighbours.begin();
      for( ; nit != neighbours.end(); ++nit ) {

        //if the neighbour is also to be removed, skip this step
        if( voxels_to_remove.find(*nit) != voxels_to_remove.end() ) {
          continue;
        }

        this->Model->GetIndexXYZ( voxfeRemodelIntType(*nit), neighXYZ );

        offset[0] = vertexXYZ[0] - neighXYZ[0];
        offset[1] = vertexXYZ[1] - neighXYZ[1];
        offset[2] = vertexXYZ[2] - neighXYZ[2];

        //generate the 'bit'
        uint index = VOXFE_GET_N26_INDEX(offset);

        //update the neighbour's connections
        gitN = this->Graph.find( *nit );
        if( gitN != this->Graph.end() ) {

          gitN->second->connectivity ^= (1<<index);  //XOR will remove the generated bit
        }
      }

      gitV->second->connectivity = 0;
      gitV->second->status = -(this->UpdateStep);
      this->NumVoxelsRemoved++;

    }  //GetNeighbs
  }  //voxels to remove


  //========================= ADD  VOXELS =============================================
  //Here we fill the complement of neighbours to selected voxels. There are a couple of
  //catches. 1) Some neighbours may have been removed in the above, so these we simply
  //want to revert. 2) Some neighbours may be being added back twice, so we simply adopt
  //a first-come-first-serve approach. 3) New neighbours of new neighbours might also
  //occur, so all new neighbours have to be added (with defaults) first before
  //connections can be set
  neighbours.clear(); non_neighbours.clear();
  IDSet voxels_added;

  vit = voxels_to_add_onto.begin();
  for( ; vit != voxels_to_add_onto.end(); ++vit ) {

    gitV = this->Graph.find( *vit );  //fixme: ?? GetNeighbours could return the iterator
    if( gitV != this->Graph.end() && ( !gitV->second->remodel ) ) continue;

    //find all the potential neighbours of this voxel to fill
    if( GetNeighbours( *vit, neighbours, non_neighbours ) ) {

      IDSetIterator nit = non_neighbours.begin();
      for( ; nit != non_neighbours.end(); ++nit ) {

        //if the prospective neighbour already exists, we can add it back
        //in if status is < 0
        gitN = this->Graph.find( *nit );
        if( gitN != this->Graph.end() ) {  //if it exists, with <0 status
          if( gitN->second->status < 0 ) {
            gitN->second->status = +(this->UpdateStep);
            this->NumVoxelsRemoved--;  //It's as if the voxel is being add back in
            this->NumVoxelsAdded++;
          }
        }
        else {  //if it doesn't exist create it...

          VertexData* vd = new (std::nothrow) VertexData();
          if( vd ) {
            vd->material = gitV->second->material;
            pair<GraphTypeIterator,bool> ret;
            ret = this->Graph.insert( pair< IDType, VertexData*>( *nit, vd ) );

            if( !ret.second ) cerr << "Attempted to add non-neighbour: " << *nit << "but failed\n";
          }
          else
            cerr << "Alloc error... could not add non-neighbour: " << *nit << "\n";
        }

        voxels_added.insert( *nit );  //accumulate all the voxels that have been added
      }

      //all neighbours should now be filled
      gitV->second->connectivity = NeighbourData[ND_ARRAY_LAST];

    }  //GetNon-Neighbs
  }  //voxels to add onto

  //process the added voxels to set their connections
  voxfeRemodelIntType index;
  vit = voxels_added.begin();
  for( ; vit != voxels_added.end(); ++vit ) {

    gitV = this->Graph.find( *vit );
    if( gitV == this->Graph.end() ) {
      cerr << "Voxel: " << *vit << " has not been added as it should....\n";
      continue;
    }
    this->Model->GetIndexXYZ( voxfeRemodelIntType(*vit), vertexXYZ );

    //voxels_added.
    for(int p=-1; p<=1; p++ ) {
      for(int q=-1; q<=1; q++ ) {
        for(int r=-1; r<=1; r++ ) {

          if( !p && !q && !r ) continue; //current voxel is not of interest!

          index = (vertexXYZ[2]+r)*Model->XY + (vertexXYZ[1]+q)*Model->X  + (vertexXYZ[0]+p);
          gitN = this->Graph.find( index );
          if( gitN != this->Graph.end() && ( gitN->second->status >= 0 ) ) {

            this->Model->GetIndexXYZ( voxfeRemodelIntType(index), neighXYZ );

            //compute the offset of the neighbour relative to the vertex
            offset[0] = neighXYZ[0] - vertexXYZ[0];
            offset[1] = neighXYZ[1] - vertexXYZ[1];
            offset[2] = neighXYZ[2] - vertexXYZ[2];

            //generate the 'bit'
            uint index = VOXFE_GET_N26_INDEX(offset);

            gitV->second->connectivity |= (1<<index);  //OR in the new bit

            gitV->second->status = +(this->UpdateStep);
          }
        }  //r
      }  //q
    }  //p


    this->NumVoxelsAdded++;
  }


}


bool voxfeRemodelControl::
GetNeighbours( const IDType& vertexID, IDSet& neighbours, IDSet& non_neighbours ) {

  GraphTypeIterator git = this->Graph.find( vertexID );
  if( git == this->Graph.end() ) {

    cerr << "** Could not access vertex data for: " << vertexID << " **\n";
    return false;
  }

  neighbours.clear();
  non_neighbours.clear();
  set<int> indices_filled, indices_unfilled;

  for( int p=0; p<ND_BLOCK_SIZE; p++ ) {  //loop through 2^0 through to 2^26
    if( git->second->connectivity & NeighbourData[p] )
      indices_filled.insert(p);
    else
      indices_unfilled.insert(p);
  }

  //convert the index values into global ID values to return
  int offset[3];
  voxfeRemodelIntType vertexXYZ[3], neighXYZ[3];
  this->Model->GetIndexXYZ( vertexID, vertexXYZ );

  set<int>::iterator sit = indices_filled.begin();
  for( ; sit!=indices_filled.end(); ++sit ) {

    VOXFE_GET_N26_XYZ( *sit, offset );

    //the index of this neighbour will then be...
    neighXYZ[0] = vertexXYZ[0] + offset[0];
    neighXYZ[1] = vertexXYZ[1] + offset[1];
    neighXYZ[2] = vertexXYZ[2] + offset[2];

    //... and the equivalent global ID
    IDType index = neighXYZ[2]*Model->XY + neighXYZ[1]*Model->X + neighXYZ[0];
    neighbours.insert( index );
  }

  sit = indices_unfilled.begin();
  for( ; sit!=indices_unfilled.end(); ++sit ) {

    VOXFE_GET_N26_XYZ( *sit, offset );

    //the index of this neighbour will then be...
    neighXYZ[0] = vertexXYZ[0] + offset[0];
    neighXYZ[1] = vertexXYZ[1] + offset[1];
    neighXYZ[2] = vertexXYZ[2] + offset[2];

    //... and the equivalent global ID
    IDType index = neighXYZ[2]*(Model->XY) + neighXYZ[1]*(Model->X) + neighXYZ[0];
    non_neighbours.insert( index );
  }

  return true;
}


bool voxfeRemodelControl::
GetFaceNeighbours( const IDType& vertexID, IDSet& neighbours ) {

  GraphTypeIterator git = this->Graph.find( vertexID );
  if( git == this->Graph.end() ) {

    cerr << "** Could not access vertex data for: " << vertexID << " **\n";
    return false;
  }

  neighbours.clear();
  set<int> indices_filled;

  //search only the face-connected neighbours (6-connected) corresponding to these bits
  //Note: hard-wiring here to bits as commented in NeighbourDataConstants
  int face_connected[ND_NUM_FACES];
  face_connected[0] = 4;  face_connected[1] = 10;
  face_connected[2] = 12; face_connected[3] = 14;
  face_connected[4] = 16; face_connected[5] = 22;

  for( int p=0; p<ND_NUM_FACES; p++ ) {  //loop through 2^0 through to 2^26
    if( git->second->connectivity & NeighbourData[ face_connected[p] ] )
      indices_filled.insert( face_connected[p] );
  }

  //convert the index values into global ID values to return
  int offset[3];
  voxfeRemodelIntType vertexXYZ[3], neighXYZ[3];
  this->Model->GetIndexXYZ( vertexID, vertexXYZ );

  set<int>::iterator sit = indices_filled.begin();
  for( ; sit!=indices_filled.end(); ++sit ) {

    VOXFE_GET_N26_XYZ( *sit, offset );

    //the index of this neighbour will then be...
    neighXYZ[0] = vertexXYZ[0] + offset[0];
    neighXYZ[1] = vertexXYZ[1] + offset[1];
    neighXYZ[2] = vertexXYZ[2] + offset[2];

    //... and the equivalent global ID
    IDType index = neighXYZ[2]*Model->XY + neighXYZ[1]*Model->X + neighXYZ[0];
    neighbours.insert( index );
  }

  return true;
}


bool voxfeRemodelControl::
GetNodeNeighbours( const IDType& vertexID, IDSet& neighbours ) {

  GraphTypeIterator git = this->Graph.find( vertexID );
  if( git == this->Graph.end() ) {

    cerr << "** Could not access vertex data for: " << vertexID << " **\n";
    return false;
  }

  neighbours.clear();
  set<int> indices_filled;

  for( int p=0; p<ND_NODE_NEIGH; p++ ) {  //loop through 2^0 through to 2^12
    if( git->second->connectivity & NeighbourData[ND_ARRAY_PRESERVE] )
      indices_filled.insert(p);
  }

  //convert the index values into global ID values to return
  int offset[3];
  voxfeRemodelIntType vertexXYZ[3], neighXYZ[3];
  this->Model->GetIndexXYZ( vertexID, vertexXYZ );

  /*typename*/ set<int>::iterator sit = indices_filled.begin();
  for( ; sit!=indices_filled.end(); ++sit ) {

    VOXFE_GET_N26_XYZ( *sit, offset );

    //the index of this neighbour will then be...
    neighXYZ[0] = vertexXYZ[0] + offset[0];
    neighXYZ[1] = vertexXYZ[1] + offset[1];
    neighXYZ[2] = vertexXYZ[2] + offset[2];

    //... and the equivalent global ID
    IDType index = neighXYZ[2]*Model->XY + neighXYZ[1]*Model->X + neighXYZ[0];
    neighbours.insert( index );
  }

  return true;
}

void voxfeRemodelControl::
WriteModel(bool removeDisconnectedElements) {

  char* buff = new char[ strlen(this->ScriptFile) + VOXFE_REMODEL_INTEGER_BUFF ];
  string outScript( this->ScriptFile );

  //test
  cout << "Starting to write: " << outScript << "\n\n";

  //We assume we want to number the output file
  sprintf(buff, VOXFE_REMODEL_PRINT_SPEC, this->UpdateStep);
  string extn;

#if 0
  //put the identifier (eg ".000001") before ".script"
  extn = VOXFE_SCRIPT_EXTN;
  outScript.replace( outScript.find(extn), extn.length(), buff);
  outScript.append(VOXFE_SCRIPT_EXTN);
#else
  //put the identifier (eg ".000001") after ".script"
  outScript.append(buff);
#endif

  string outModel( this->Model->voxelFile );
  extn = "./";
  outModel.replace( outModel.find(extn), extn.length(), "");  //remove the opening "./"
  outModel.append( buff );

  rvsWriteVoxFEScriptStub( outScript.c_str(), outModel.c_str(),
                           this->Model->voxSize, this->Model->scanDim,
                           "KSPCG PCJACOBI",
                           true, buff );
  delete[] buff;

  //now write the model
  ofstream modelFile( outModel.c_str() );
  if( !modelFile.is_open() ) {
    cerr << "\n\n ** Cannot open model file ** \n\n";
    return;
  }

  voxfeRemodelIntType count = 1, vertexXYZ[3],
    count_disconnected=0, count_neg_status=0, count_pos_status=0;
  GraphTypeIterator git;

//Go through counts to make sure the output model count (at the *head* of the file is good
//and to help the user via a log file, with counts of voxels removed, added etc
#define VOXFE_COUNT_FIRST
#ifdef VOXFE_COUNT_FIRST
  voxfeRemodelIntType count_valid=0;
  git = this->Graph.begin();

  //determine the valid voxel count first
  for( ; git != this->Graph.end(); ++git ) {

    //test
    if( !git->second->component ) continue;
    //if( !git->second->connectivity ) continue;
    if(  git->second->status < 0 )   continue;

    count_valid++;
  }
  modelFile << "\n" << count_valid << "\n";


  git = this->Graph.begin();
  if( removeDisconnectedElements ) {

    cout << "\n\n *** Writing model file, removing disconnected *** \n\n";

    for( ; git != this->Graph.end(); ++git ) {

      //log data
      if( !git->second->connectivity ) count_disconnected++;
      if(  git->second->status < 0 )   count_neg_status++;
      else if(  git->second->status > 0 )   count_pos_status++;

      if(  (!git->second->component) ||
           /*(!git->second->connectivity) || */   //fixme: does not seem to be needed ie solver copes
           ( git->second->status < 0 ) ) continue;

      this->Model->GetIndexXYZ( git->first, vertexXYZ );
      modelFile << count << " "
                << short(git->second->material) << " "
                //<< git->first << " "
                //<< git->second->component << " "
               << (vertexXYZ[0] + voxfeVoxelOffset) << " "
               << (vertexXYZ[1] + voxfeVoxelOffset) << " "
               << (vertexXYZ[2] + voxfeVoxelOffset) << "\n";

      count++;
    }
  }
  else {

    cout << "\n\n *** Writing model file, without removing disconnected *** \n\n";

    for( ; git != this->Graph.end(); ++git ) {

      //test
      if( !git->second->connectivity ) count_disconnected++;
      if(  git->second->status < 0 )   count_neg_status++;

      if( (!git->second->connectivity) || ( git->second->status < 0 ) ) continue;

      this->Model->GetIndexXYZ( git->first, vertexXYZ );
      modelFile << count << " " << short(git->second->material)  << " "
               << (vertexXYZ[0] + voxfeVoxelOffset) << " "
               << (vertexXYZ[1] + voxfeVoxelOffset) << " "
               << (vertexXYZ[2] + voxfeVoxelOffset) << "\n";

      count++;
    }
  }


#else

  modelFile << "\n" << (this->Graph.size() - this->NumVoxelsRemoved) << "\n";
  git = this->Graph.begin();
  for( ; git != this->Graph.end(); ++git ) {

    //test
    if( !git->second->connectivity ) count_disconnected++;
    if(  git->second->status < 0 )   count_neg_status++;

    if( /*(!git->second->connectivity) ||*/ ( git->second->status < 0 ) ) continue;

    this->Model->GetIndexXYZ( git->first, vertexXYZ );
    modelFile << count << " " << short(git->second->material)  << " "
             << (vertexXYZ[0] + voxfeVoxelOffset) << " "
             << (vertexXYZ[1] + voxfeVoxelOffset) << " "
             << (vertexXYZ[2] + voxfeVoxelOffset) << "\n";

    count++;
  }
#endif

  //output log data for user
  string log = outScript + ".log";
  ofstream fout( log.c_str(), std::ios::app );

  fout << "Writing out: " << outScript << "\n\n";
  fout << "\n\nGraph.size(): " << this->Graph.size()
       << "\n\tNum voxels removed: "         << this->NumVoxelsRemoved
       << "\n\tNum voxels added:   "         << this->NumVoxelsAdded
       << "\n\tNum disconnected voxels: "    << count_disconnected
       << "\n\tNum negative status voxels: " << count_neg_status
       << "\n\tNum positive status voxels: " << count_pos_status
       << "\n\tOutput voxels: " << (this->Graph.size() - this->NumVoxelsRemoved)
       << "   Final count: " << (count-1)
       << "\n\n&&&&&&&&&&&&&&&&&&&&&&&&&\n\n" ;

  fout.close();

  modelFile.close();
}


bool voxfeRemodelControl::
CreateConnections( const IDType& metisElementID,
                const vector<IDType>& metisNeighIDs,
                const voxfeRemodelIntType* elementIDList,
                const unsigned char* materialList ) {

  if( !Model || !elementIDList || !materialList ) return false;

  IDType vertexID = elementIDList[metisElementID-1], neighID;
  unsigned char vertexMaterial = materialList[metisElementID-1];

#ifdef VOXFE_REMODEL_DEBUG
  cout << "CreateConnections for metisID: " <<  metisElementID << " elementID: " << elementIDList[metisElementID-1]
       << " material: " <<  (short)vertexMaterial << "\n";
#endif

  pair<GraphTypeIterator,bool> ret;

  //check if vertex exists in the graph, or else create an entry
  GraphTypeIterator git = this->Graph.find( vertexID );
  if( git == this->Graph.end() ) {

    VertexData* vd = new (std::nothrow) VertexData();
    if( !vd ) {
      cerr << "** Could not create the Vertex Data **\n";
      return false;
    }

    //set material and remodel flag
    vd->material = vertexMaterial;
    vd->remodel  = (Model->materialMap[int(vertexMaterial)]).RemodelFlag;  //fixme: yuk -- may be adding to Model->materialMap

    ret = this->Graph.insert( pair< IDType, VertexData*>( vertexID, vd ) );
    if( ret.second ) git = ret.first;
    else {
      // ? this should have created the entry where it did not exist
      cerr << "** Could not create the graph entry **\n";
      return false;
    }
  }

  voxfeRemodelIntType vertexXYZ[3], neighXYZ[3];
  int offset[3];
  this->Model->GetIndexXYZ( voxfeRemodelIntType(vertexID), vertexXYZ );

#ifdef VOXFE_REMODEL_DEBUG
  cout << "Adding in connec data for vertex: " << vertexID
       << " (" << vertexXYZ[0] << "," << vertexXYZ[1] << "," << vertexXYZ[2]
       << ")\n";
#endif


  //add a bit for every connecting element
  git->second->connectivity = 0;
  for( int p=0; p<metisNeighIDs.size(); p++ ) {

    neighID = elementIDList[metisNeighIDs[p]-1];

    this->Model->GetIndexXYZ( voxfeRemodelIntType(neighID), neighXYZ );

    offset[0] = neighXYZ[0] - vertexXYZ[0];
    offset[1] = neighXYZ[1] - vertexXYZ[1];
    offset[2] = neighXYZ[2] - vertexXYZ[2];

    //generate the 'bit'
    uint index = VOXFE_GET_N26_INDEX(offset);

    //add it in to the graph entry
    git->second->connectivity |= (1 << index);

#ifdef VOXFE_REMODEL_DEBUG
    cout << " -- neighbour: " << neighID
         << " gives index: "  << index
         << " (bit: "    << std::hex << (1 << index)
         << ") at loc: " << std::dec
         << neighXYZ[0] << "," << neighXYZ[1] << "," << neighXYZ[2]
         << " (at offset: "
         << offset[0] << "," << offset[1] << "," << offset[2]
         << ")\n";
#endif
  }


#ifdef VOXFE_REMODEL_DEBUG
  cout << "Setting all connec data: " << std::hex << git->second->connectivity
       <<  std::dec << "\n\n";
#endif

  return true;
}

void voxfeRemodelControl::
CreateNeighbourDataBits() {

  if( !NeighbourData ) NeighbourData = new uint[ ND_ARRAY_SIZE ];

  NeighbourData[ND_ARRAY_LAST] = 0;

#ifdef VOXFE_REMODEL_DEBUG
  cout <<"Setting neighbour bits:\n";
#endif
  for( int p=0; p<ND_BLOCK_SIZE; p++ ) {  //store 27 bit-values values for 2^0 through to 2^26

    NeighbourData[p] = ( 1 << p );  //NOTE: (1 << 0) ie 1, represents the first bit with bitshift 0
    NeighbourData[ND_ARRAY_LAST] |= NeighbourData[p];
#ifdef VOXFE_REMODEL_DEBUG
    cout << p << " : " << NeighbourData[p] <<"\n";
#endif
  }
#ifdef VOXFE_REMODEL_DEBUG
  cout << "\n\n";
#endif

  //create some 'special' bit patterns
  NeighbourData[ND_ARRAY_LAST] -= (1 << 13); // all neighbours except the current voxel

  //all the neighbour voxels surrounding a given node
  NeighbourData[ND_ARRAY_PRESERVE] = NeighbourData[0]  |
                                     NeighbourData[1]  |
                                     NeighbourData[3]  |
                                     NeighbourData[4]  |
                                     NeighbourData[9]  |
                                     NeighbourData[10] |
                                     NeighbourData[12] |
                                     NeighbourData[13];

}


void voxfeRemodelControl::
InializeNodeToDisplacementMap() {

  voxfeRemodelIntType vertexXYZ[3], a[NODESPERELEMENT];
  double* d = 0;

  NodeToDisplacementMap.clear();

  GraphTypeIterator git = this->Graph.begin();
  for( ; git != this->Graph.end(); ++git ) {

    this->Model->GetIndexXYZ( voxfeRemodelIntType(git->first), vertexXYZ );
    rvsGetVoxelNodes( *(this->Model), vertexXYZ, a );

    for( int p=0; p<NODESPERELEMENT; p++ )
      NodeToDisplacementMap.insert( pair< IDType, double* >( a[p], d ) );
  }

}

void voxfeRemodelControl::
UpdateNodeToDisplacementMap() {

  //fixme: do check on num nodes
  voxfeRemodelIntType count;
  if( this->Displacement ) {

    NodeToDisplacementMapIterator nit = this->NodeToDisplacementMap.begin();
    for( count=0 ; nit != this->NodeToDisplacementMap.end(); ++nit, count+=3 ) {
      nit->second = (this->Displacement + count);
    }
  }
}


void voxfeRemodelControl::
SelectBySED(bool surface_only) {

  this->voxels_to_add_onto.clear();
  this->voxels_to_remove.clear();
  if( !ComputeSED(surface_only) ) {
    cerr << "Could not compute SED\n";
    return;
  }

  double localThreshold[2];
  localThreshold[0] = this->Threshold[0];
  localThreshold[1] = this->Threshold[1];

  //===============================================================
  // If we're using adaptive remodelling, we set up a histogram, to reset the thresholds
  // Experience suggests most of the data lie in a spike in the first few bins of the histogram
  // (which is nearly "half-normal" but not exactly, hence we need the histogram),
  // so use logs if min > 0. This saves having to have many millions of bins
  if( this->AdaptiveThresholding ) {

    long nbins = VOXFE_NUM_HISTOGRAM_BINS;
    double hmin = this->SED_stats.Min(), hmax = this->SED_stats.Max();

    bool useLogs = false;
    if( hmin > 0 ) {
      hmin = log10( hmin );
      hmax = log10( hmax );
      useLogs = true;
    }

    RunningHistogram histogram( nbins, hmin, hmax );

    //run through the graph adding data to the histogram
    double sum_to_min = 0.; double sum_to_max = 0.;
    GraphTypeIterator git = this->Graph.begin();
    if( useLogs ) {
      for( ; git != this->Graph.end(); ++git ) {

        if( ( git->second->status < 0 ) || ( !git->second->remodel ) ) continue;
        if( surface_only && (!this->IsOnSurface( git->first )) ) continue;

        histogram.AddValue( log10( git->second->SED ) );
      }

      localThreshold[0] = pow( 10.0, histogram.GetBinValueAtCumulativeDensity( this->Threshold[0], sum_to_min ) );
      localThreshold[1] = pow( 10.0, histogram.GetBinValueAtCumulativeDensity( this->Threshold[1], sum_to_max ) );
    }
    else {
      for( ; git != this->Graph.end(); ++git ) {

        if( ( git->second->status < 0 ) || ( !git->second->remodel ) ) continue;
        if( surface_only && (!this->IsOnSurface( git->first )) ) continue;

        histogram.AddValue( git->second->SED );
      }

      localThreshold[0] = histogram.GetBinValueAtCumulativeDensity( this->Threshold[0], sum_to_min );
      localThreshold[1] = histogram.GetBinValueAtCumulativeDensity( this->Threshold[1], sum_to_max );

    }

    //debug
    cout << "\nSetting adaptive thresholds: \n";
    cout << "Lower threshold: " << localThreshold[0] << " (sum: " << sum_to_min << ")\n";
    cout << "Upper threshold: " << localThreshold[1] << " (sum: " << sum_to_max << ")\n\n";

    ofstream fout("DumpHist.txt");
    histogram.Print( fout );
    fout.close();
  }


  GraphTypeIterator git = this->Graph.begin();
  for( ; git != this->Graph.end(); ++git ) {

	  //if 'removed' voxel or not for remodelling, pass
	  if( ( git->second->status < 0 ) || ( !git->second->remodel ) ) continue;

	  //Are we only checking surface elements?
	  if( surface_only && (!this->IsOnSurface( git->first )) )
	    continue;


   //work through voxels to remove
    if( git->second->SED <= localThreshold[0] ) {

      this->voxels_to_remove.insert( git->first );
    }
    else if( git->second->SED >= localThreshold[1] ) {

      this->voxels_to_add_onto.insert( git->first );
    }

	}

}


bool voxfeRemodelControl::
ComputeSED(bool surface_only) {

  if( NodeToDisplacementMap.empty() ) {
    cerr << "Need to construct node--displacement mapping\n";
    return false;
  }

  this->SED_stats.Clear();

  voxfeRemodelIntType vertexXYZ[3], neighXYZ[3], a[NODESPERELEMENT];
  int offset[3];

	//... needed to compute the strains
  voxfeRecomputeStrainFilter::BMatrix B;
	voxfeRecomputeStrainFilter::uVector u;
	voxfeRecomputeStrainFilter::eVector e;
  voxfeRecomputeStrainFilter recomp;
	recomp.computeVoxelGradientMatrix( this->Model->voxSize[0], B );
  double MOR, LAM ;

  //debug
  ofstream fout("DumpSED.txt");


  GraphTypeIterator git = this->Graph.begin();
	for( ; git != this->Graph.end(); ++git ) {

	  //if 'removed' voxel or not for remodelling, pass
	  if( ( git->second->status < 0 ) || ( !git->second->remodel ) ) continue;

	  //Are we only checking surface elements?
	  if( surface_only && (!this->IsOnSurface( git->first )) ) {

	    git->second->SED = VOXFE_REMODEL_DEFLT_SED; //set a default non-value (legal, but unlikely)
	    continue;
	  }

	  //get material props
	  int material = git->second->material;
    map< int, MaterialType >::iterator it = this->Model->materialMap.find( material );
    if( it != this->Model->materialMap.end() ) {
      MOR = recomp.getModulusOfRigidity( it->second.YoungsModulus, it->second.PoissonsRatio );
      LAM = recomp.getLambda( it->second.YoungsModulus, it->second.PoissonsRatio );
    }
    else {
      cerr << "Cannot determine material properties for: " << material << "\n";
      continue;
    }

    this->SetVoxelDisplacements( git, u );
    recomp.getStrainTensor( B, u, e ); //compute the basic strains
	  recomp.computeStrainEnergyDensity2( e.data(), MOR, LAM, git->second->SED );

    //track the stats for SED
    this->SED_stats.Push( git->second->SED );

    //debug
    fout << git->second->SED << "\n";

   //work through voxels to remove
   /*
    if( git->second->SED <= this->Threshold[0] ) {

      this->voxels_to_remove.insert( git->first );
    }
    else if( git->second->SED >= this->Threshold[1] ) {

      this->voxels_to_add_onto.insert( git->first );
    }
    */

	}

	//debug -- print stats:
  cout << "Min: " << this->SED_stats.Min() << "  "
       << "Max: " << this->SED_stats.Max() << "  "
       << "Mean: " << this->SED_stats.Mean() << "  "
       << "Var: "  << this->SED_stats.Variance() << "  "
       << "(Count: " << this->SED_stats.NumDataValues() << ")\n";

  //debug
  fout.close();

  return true;

}



void voxfeRemodelControl::
SetVoxelDisplacements( GraphTypeIterator& git, voxfeRecomputeStrainFilter::uVector& u ) {

  voxfeRemodelIntType vertexXYZ[3], a[NODESPERELEMENT];

  this->Model->GetIndexXYZ( voxfeRemodelIntType(git->first), vertexXYZ );
  rvsGetVoxelNodes( *(this->Model), vertexXYZ, a );

  //set the displacements around the voxel
  for( int m=0; m<NODESPERELEMENT; m++ ) {

    NodeToDisplacementMapIterator nit = NodeToDisplacementMap.find( a[m] );
    if( nit != NodeToDisplacementMap.end() ) {

      int m3 = 3*m;
      u(m3  ) = nit->second[0];
      u(m3+1) = nit->second[1];
      u(m3+2) = nit->second[2];
    }
  }
}


void voxfeRemodelControl::
SetThresholds( const double& lowerSEDThreshold, const double& upperSEDThreshold, bool adaptive ) {

  //fixme: add sanity checks?
  this->Threshold[0] = lowerSEDThreshold;
  this->Threshold[1] = upperSEDThreshold;
  AdaptiveThresholding = adaptive;

}


void voxfeRemodelControl::
SetConnectedComponents( bool use26NeighConnectivity, bool selectLargest ) {

  //reset the equivalence table and component data
  ResetComponents();

  //1st pass: assign labels for everything in the graph...
  AssignComponentLabels( use26NeighConnectivity );

#ifdef VOXFE_REMODEL_DEBUG
  PrintEquivTable( "EquivTable-Full.txt" );
#endif

  //Find parent component labels
  FlattenEquivalenceTable();

#ifdef VOXFE_REMODEL_DEBUG
  PrintEquivTable( "EquivTable-Flat.txt" );
#endif

  //Create a map to count elements in the largest group
  map<LabelType,long> labelCount;

  //2nd pass. Iterate through the graph, assigning parent labels and keeping count
  LabelType parent;
  LabelSetIterator sit;
  GraphTypeIterator git = this->Graph.begin();
  for( ; git != this->Graph.end(); ++git ) {

    parent = 0;
    sit = this->EquivTable[ git->second->component ].begin();
    if( sit != this->EquivTable[ git->second->component ].end() )
      parent = *sit;
    else {
      cerr << "\n\n** Cannot find a parent for element: " << git->first << " -- IGNORING! **\n\n";
      continue;
    }

    git->second->component = parent;

    //keep count
    if( labelCount.find( parent ) == labelCount.end() )
      labelCount.insert( pair<LabelType,long> (parent, 0) );

    labelCount[ parent ]++;
   }


  //if select==true, zero the labels of everything other than the largest object
  //find the most frequent label
  map<LabelType,long>::iterator mit = labelCount.begin();
  long max_count = 0;
  LabelType largest_obj = 0;
  for( ; mit != labelCount.end(); ++mit ) {
    if( mit->second > max_count ) {
      max_count   = mit->second;
      largest_obj = mit->first;
    }
  }

#ifdef VOXFE_REMODEL_DEBUG
  ofstream fcount( "LabelCounts.txt");
  mit = labelCount.begin();
  for( ; mit != labelCount.end(); ++mit ) {
    fcount << mit->first << " : " << mit->second << "\n";
  }
  fcount.close();
#endif


  if( selectLargest ) {

    //Final pass, if desired to zero components of all but the largest object in the graph
    GraphTypeIterator git = this->Graph.begin();
    for( ; git != this->Graph.end(); ++git ) {

      if( git->second->component != largest_obj )
        git->second->component = 0;
    }
  }

}



void voxfeRemodelControl::
AssignComponentLabels( bool use26NeighConnectivity ) {

  //1st pass of CC algorithm, loop through graph, find neighbours of each element
  //and assign the local parent to all labels in the group
  IDSet neighbours, non_neighbours;

  //cout << "Assigning components:\n";
  GraphTypeIterator git = this->Graph.begin();

  if( use26NeighConnectivity ) {

    for( ; git != this->Graph.end(); ++git ) {

      //if the graph has been remodelled, ignore these
      if(  git->second->status < 0 )   continue;

      if( !this->GetNeighbours( git->first, neighbours, non_neighbours ) ) continue;

      //add the current voxel too
      neighbours.insert( git->first );

      //cout << "Element: " << git->first << " got: " << neighbours.size()-1 << " neighbours...\n";
      SetToLocalParentLabel( neighbours );
    }

  } else {

    for( ; git != this->Graph.end(); ++git ) {

      //if the graph has been remodelled, ignore these
      if(  git->second->status < 0 )   continue;

      if( !this->GetFaceNeighbours( git->first, neighbours ) ) continue;

      //add the current voxel too
      neighbours.insert( git->first );

      //cout << "Element: " << git->first << " got: " << neighbours.size()-1 << " neighbours...\n";
      SetToLocalParentLabel( neighbours );
    }

  }
}


void voxfeRemodelControl::
FlattenEquivalenceTable() {

  //Work backwords through the table looking for the best high to low
  //(child to parent) connecting label
  LabelSet parents;
  LabelType parent, i = this->EquivTable.size() - 1;
  for(  ; i > 0; i-- ) {   //NB with note use of unsigned type, >=0 fails here !!!

    //a small optimisation here, if i has already been looked up, it may have been
    //updated and lie in 'parents', so skip it if it does
    //if( parents.find(i) != parents.end() ) continue;

    parents.clear();

    //Find all the labels associated with this upwards (towards the parent)
    FindAllParentLabels( i, parents );

    //set all labels in 'parents' to parent
    LabelSetIterator sit = parents.begin();
    if( sit != parents.end() )
      parent = *sit;
    else {
      cout << "** Cannot find a parent for: " << i << " -- IGNORING! **\n";
      continue;
    }

    // 'flatten' the table by replacing with the parent value
    ++sit;
    while( sit != parents.end() ) {

      this->EquivTable[*sit].clear();
      this->EquivTable[*sit].insert(parent);
      ++sit;
    }

    //fixme: we could retain the parents set at this point to reduce the search through EquivTable
    //but, in the interest of keeping memory needs down, and assuming the later searches for the
    //same parent should be much reduced, we leave things to run....
  }
}


void voxfeRemodelControl::
FindAllParentLabels( const LabelType& J, LabelSet& parents ) {

  parents.insert(J);  //include self

#ifdef VOXFE_REMODEL_DEBUG
  if( J >= this->EquivTable.size() ) {

    cerr << "** Label J: " << J << " exceeds table size! (" << this->EquivTable.size() <<") **\n";
    return;
  }
#endif

  //Return when EquivTable[J] set is empty or when only contains J
  //(Assume all neighbour connecting labels are stored from high (child)
  //to low (parent)
  LabelSetIterator sit = this->EquivTable[J].begin();
  for( ; sit != this->EquivTable[J].end(); ++sit ) {

    if( J == *sit ) continue;
    else {
      this->FindAllParentLabels( *sit, parents );
    }
  }
}


bool voxfeRemodelControl::
SetToLocalParentLabel( IDSet& elems ) {

  if( elems.empty() ) return false;

  vector< GraphTypeIterator > pElems;
  LabelSet labels;
  LabelSetIterator sit;
  IDSetIterator it;

  //gather the labels and iterators for all elems
  GraphTypeIterator git;
  for( it=elems.begin(); it!=elems.end(); ++it ) {

    git = this->Graph.find( *it );
    if( git != this->Graph.end() ) {

      //gather all the labels possibly associated with this group of elements
      if( git->second->component )
        labels.insert( git->second->component );

      pElems.push_back( git );
    }
  }

  //select a parent (lowest) label
  LabelType parent = 0;
  if( labels.empty() ) {
    GetNextLabel( parent );
    //test    //cout << "Using next label as parent: " << parent << "\n";
  }
  else {
    sit = labels.begin();
    parent = *sit;
    //test   //cout << "Found parent: " << parent << "  among: " << labels.size() << " labels: ";
    while( sit != labels.end() ) {

      this->EquivTable[*sit].insert(parent);
      //test //cout << *sit << " ";
      ++sit;
    }
  }

  //could not find a parent value, something went wrong -- what??
  if( !parent ) {
    cerr << "\n\n** Could not find the parent label... **\n\n";
    return false;
  }

  //assign all graph elements to this component
  for( int k=0; k < pElems.size(); k++ )
    pElems[k]->second->component = parent;

  return true;
}


void voxfeRemodelControl::
ResetComponents() {

  //reset all component values
  GraphTypeIterator git = this->Graph.begin();
  for( ; git != this->Graph.end(); ++git )
    git->second->component = 0;

  this->EquivTable.clear();

  //initiate EquivTable(we don't want to consider EquivTable[0])
  LabelType l;
  this->GetNextLabel( l );
}


void voxfeRemodelControl::GetNextLabel( LabelType& label ) {

  label = EquivTable.size();
  LabelSet s;

  s.insert( label );  //fixme: this seems unnecessary, but causes much more checking to be needed
  EquivTable.push_back( s );
}


void voxfeRemodelControl::
PrintEquivTable( const char* tableFile ) {

  ofstream fout( tableFile );
  if( !fout.is_open() ) return;

  for( LabelType i = 0; i < this->EquivTable.size(); i++ ) {

    fout << i << " : ";
    LabelSetIterator sit = this->EquivTable[i].begin();
    if( sit != this->EquivTable[i].end() )
      fout << *sit << " ";

    fout << "\n";
  }

  fout.close();
}


void voxfeRemodelControl::
PrintGraph( std::ostream& out, bool withNeighbours ) {

  voxfeRemodelIntType vertexXYZ[3], neighXYZ[3];
  IDSet neighbours, non_neighbours;
  IDSetIterator sit;
  GraphTypeIterator git = this->Graph.begin();
  for( ; git != this->Graph.end(); ++git ) {

    this->Model->GetIndexXYZ( voxfeRemodelIntType(git->first), vertexXYZ );
    out << "Element: " << git->first
        << " (" << vertexXYZ[0] << "," << vertexXYZ[1] << "," << vertexXYZ[2]
        << ") Connectivity: " << std::hex << git->second->connectivity
        << std::dec
        << " Material: " << (short)git->second->material
        << " Status: "   << git->second->status
        << " Remodel: "  << (short)git->second->remodel
        << " SED: "  << git->second->SED;

    if( !withNeighbours )
      out << "\n";
    else {
      out << "\n -- Neighbs: ";

			if( this->GetNeighbours( git->first, neighbours, non_neighbours ) ) {

			  sit = neighbours.begin();
			  for( ; sit != neighbours.end(); ++sit ) {

			    this->Model->GetIndexXYZ( voxfeRemodelIntType(*sit), neighXYZ );
			    out << *sit << " ("
			        << neighXYZ[0] << "," << neighXYZ[1] << "," << neighXYZ[2]
			        << ")  ";
			  }

			  out << "\n -- Non_neighbs: ";

			  sit = non_neighbours.begin();
			  for( ; sit != non_neighbours.end(); ++sit ) {

			    this->Model->GetIndexXYZ( voxfeRemodelIntType(*sit), neighXYZ );
			    out << *sit << " ("
			        << neighXYZ[0] << "," << neighXYZ[1] << "," << neighXYZ[2]
			        << ")  ";
			  }
			  out << "\n";
			}
    }

  }

  out << "\n\nNum Voxels Removed: " << this->NumVoxelsRemoved << "\n\n";

#ifdef VOXFE_REMODEL_DEBUG
  out << "\nDisplacements:\n";
  NodeToDisplacementMapIterator nit = NodeToDisplacementMap.begin();
  for( ; nit != NodeToDisplacementMap.end(); ++nit ) {

    out << nit->first << "  ";
    if( nit->second )
      out << nit->second[0] << " " << nit->second[1] << " " << nit->second[2] << "\n";
    else
      out << "\n";
  }
#endif

}




