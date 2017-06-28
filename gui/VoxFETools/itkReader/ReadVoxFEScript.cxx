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

#include "ReadVoxFEScript.h"

const char* VOXFESCRIPTHEADER[] = {
  "LOAD_MATERIALS_FILE ./materials.txt\n",
  "VOXEL_SIZE ",
  " 1\n\n",

  "ALGORITHM_FEA ",  //JAPCG
  "MAX_ITER 2000\n",
  "TOLERANCE 1e-09\n",
  "LOAD_MCTSCAN ",  //eg. 128 128 16  ./hole-cylinder.model

  "\nSELECTION_OF_NODES\n",
  "    SELECT_NODE_FILE ./constraints.txt\n",
  "COMPUTE_SED\n\n",

  "PRINT_X ", "./displacement.txt",
  "PRINT_FULL_STRAIN ", "./fullstrains.txt",
  "REACTION ", "./reactionforces.txt"
};

/*
const char* PETSCSCRIPTHEADER[] = {
    "SET_VOXEL_SIZE ",  // 0.1 0.1 0.1  1
    " 1\n\n",
    "LOAD_MATERIALS ./materials.txt\n\n",

    "SET_ALGORITHM_FEA KSPCG PCJACOBI\n",
    "SET_MAX_ITER 2000\n",
    "SET_TOLERANCE 1e-09\n\n",

    "LOAD_MODEL ", //520 340 90  ./bar.model
    "LOAD_CONSTRAINTS ./constraints.txt\n\n",

    "PRINT_DISPLACEMENTS ./displacement.txt"
};*/



//Write out VOXFESCRIPTHEADER as scriptfile, nominating the model
void rvsWriteVoxFEScriptStub( const char* scriptfile, const char* modelfile,
                           double spacing[3], int size[3],
                           const char* algorithm,
                           bool displacement_only, const char* counter ) {

  ofstream sout( scriptfile );
  int hcount = 0;
  for( int q=0; q<2; q++ ) sout << VOXFESCRIPTHEADER[ hcount++ ];

  sout << spacing[0] << " " << spacing[1] << " " << spacing[2] << " ";   //VOXEL_SIZE
  sout << VOXFESCRIPTHEADER[ hcount++ ];                                 //USF

  sout << VOXFESCRIPTHEADER[ hcount++ ];                                 //ALGORITHM
  sout << algorithm << "\n";

  for( int q=0; q<3; q++ ) sout << VOXFESCRIPTHEADER[ hcount++ ];        //MAX_ITER -- LOAD_MCTSCAN

  //the file name given needs to be local here
  sout << size[0] << " " << size[1] << " " << size[2] << "  ./" << modelfile;  //MODEL
  sout << "\n";  //assume any necessary count identifier is already in the name

  cout << "sizeof(VOXFESCRIPTHEADER): " << sizeof(VOXFESCRIPTHEADER) << "  ";
  cout << (sizeof(VOXFESCRIPTHEADER)/sizeof(const char*)) << "\n";

  for( int q=0; q<3; q++ ) sout << VOXFESCRIPTHEADER[ hcount++ ];        //SELECTION -- COMPUTE_SED

  //PRINT_X displacements
  sout << VOXFESCRIPTHEADER[ hcount++ ];
  if( counter )
    sout << VOXFESCRIPTHEADER[ hcount++ ] << counter << "\n";
  else
    sout << VOXFESCRIPTHEADER[ hcount++ ] << "\n";

  if( !displacement_only ) {
    //dump the last n lines
    for( ; hcount < (sizeof(VOXFESCRIPTHEADER)/sizeof(const char*)); hcount+=2 ) {
      sout << VOXFESCRIPTHEADER[ hcount ];
      if( counter )
        sout << VOXFESCRIPTHEADER[ hcount+1 ] << counter << "\n";
      else
        sout << VOXFESCRIPTHEADER[ hcount+1 ] << "\n";
    }
  }
  sout.close();
}


// ====================================================================================================
// ============================== ImageGroupData ======================================================
// ====================================================================================================

ImageGroupData::ImageGroupData() {

  for( int i=0; i<MAX_MATERIALS; i++ ) Tally[i] = 0; //initialize the state array
  Group.clear();
  VoxelSize = -1;  //not needed for specified types such as 'mhd'
                   //but must be supplied by user for png, bmp etc
  BorderOffset[0] = BorderOffset[1] = BorderOffset[2] = 0;
  TallyDone = false;
}

ImageGroupData::~ImageGroupData() {
  //delete[] Group;
}

  /** Maintain a list of labels present (and add single groups if any) */
void ImageGroupData::TallyGroups( const set<short>& labels_present ) {

  if( TallyDone ) return;  //why would you want to repeat ??

  //initialize the state array
  //for( int i=0; i<MAX_MATERIALS; i++ ) Tally[i] = 0;

  set<short>::const_iterator it = labels_present.begin();
  for( ; it != labels_present.end(); ++it )
    Tally[*it] = LABEL_PRESENT;

  set<short> singles;
  GetNonGroups( singles );

  //add in single labels (not already in a group) as additional groups
  char buff[32]; //fixme: sstream?
  short remodel = 1;
  for( it=singles.begin() ; it!=singles.end(); ++it ) {

    sprintf( buff, "%04d", *it );
    Group.push_back( ImageGroupData::GroupDef( string(buff), *it, *it, remodel ) );
  }

  TallyDone = true;
}

void ImageGroupData::WriteMaterialsFile( const char* MatFileName, const RefMaterialMap& RefMaterials ) {

  if( !TallyDone ) {
    cerr << "Group data has not been tallied with labels which actually exist\n";
    return; //fixme: output all possible instead?
  }

  //expand the tally into a full list of names
  vector<string> nameTally;
  for( int m=0; m<MAX_MATERIALS; ++m ) {
    nameTally.push_back("00");  //a 'null' material
  }

  //put real materials in the name tally
  for( short p=0; p<Group.size(); p++ ) {
    string s = Group[p].MaterialName;
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);  //set nameTally all to lower case
    for( short q=Group[p].Start; q<=Group[p].Stop; q++ ) {
      nameTally[q] = s;
    }
  }

  ofstream mat_out( MatFileName  );
  for( int m=0; m<MAX_MATERIALS; ++m ) {

    if( !Tally[m] ) continue;

    double ym = VOXFE_YOUNGSMODULUS;
    double pr = VOXFE_POISSONSRATIO;

    RefMaterialMapConstIterator mit = RefMaterials.find( nameTally[m].c_str() );
    if( mit != RefMaterials.end() ) {
      ym = mit->second.YoungsModulus;
      pr = mit->second.PoissonsRatio;
    }

    mat_out << m  << "  "
	    << ym << "  "
	    << pr << "  "
	    << (((Tally[m] & LABEL_REMOD)>0)?1:0)
	    << std::endl;
  }
  mat_out.close();
}


void ImageGroupData::GetNonGroups( set<short>& non_groups ) {

  non_groups.clear();

  //if the tally shows a label is in a group, promote the flag to show
  //whether it can be remodelled or not //fixme: is this still needed?
  for( short p=0; p< Group.size(); p++ ) {
    for( short q=Group[p].Start; q<=Group[p].Stop; q++ ) {

      if( Tally[q] == LABEL_PRESENT ) {
	if( Group[p].RemodelFlag > 0 )
	  Tally[q] = LABEL_PRESENT_IN_GROUP_REMOD;  //LABEL_PRESENT|LABEL_IN_GROUP|LABEL_REMOD
	else
	  Tally[q] = LABEL_PRESENT_IN_GROUP_NO_REMOD;
      }
    }
  }

  //find labels that occur outside known groups -- strays?
  for( short i=0; i<MAX_MATERIALS; i++ ) {

    //everything else here must be outside the user-defined groups
    if( Tally[i] == LABEL_PRESENT ) {
      Tally[i] |= LABEL_REMOD;
      non_groups.insert(i);
    }
    //else if( Tally[i] > 1 ) Tally[i] = 1; //reset to '1' for next call
  }
}





