
#include "voxfeRemodelControl.h"



#define VOXFE_REMODEL_ORIG_MAIN
#ifdef VOXFE_REMODEL_ORIG_MAIN
int main( int argc, char* argv[] ) {

  //this 'main' is intended to be used once before quitting creating only
  //a single remodelled file

  if( argc < 9 ) {
    cout << "Usage: " << argv[0]
        << "VoxFEFile.script "
           "ElementMapFile.elmap "
           "ElementGraphFile.graph "
           "AdaptiveThresholds<Y/N>"
           "LowerThreshold "
           "UpperThreshold "
           "DisplacementFile.txt "
           "PetScSolver<Y/N> "
           "[Use-26-Neighbour-Connectivity=0]";
    cout << "\n\n";
    exit(-1);
  }

  bool useAdaptiveThresholds = false;
  if( (argv[4][0] == 'Y') || (argv[4][0] == 'y') )
    useAdaptiveThresholds = true;


  bool usePetScSolver = false;
  if( (argv[8][0] == 'Y') || (argv[8][0] == 'y') )
    usePetScSolver = true;


  bool use26NeighConnectivity = 0;
  if( argc == 10 )
    use26NeighConnectivity = atoi( argv[9] );

  voxfeRemodelControl vrc( argv[1], argv[2], argv[3] );

  vrc.SetThresholds( atof( argv[5] ), atof( argv[6] ), useAdaptiveThresholds );

  //ofstream pg1( "PrintGraph-Initial.txt");
  //vrc.PrintGraph( pg1, false );
  //pg1.close();

#ifndef VOXFE_GET_CONNECTIVITY_MODEL
  if( !vrc.SetDisplacements( argv[7], usePetScSolver, true, use26NeighConnectivity ) ) {
    cerr << "** \n\n Setting displacements failed for: " << argv[7] << " **\n\n";
    return -1;
  }
#else
  vrc.SetConnectedComponents(use26NeighConnectivity, true);
  vrc.WriteModel(true);
#endif

  //ofstream pg2( "PrintGraph-Final.txt");
  //vrc.PrintGraph( pg2, false );
  //pg2.close();

  cout << "\n\n** Remodelling done.... have a nice day ;-) **\n\n";

  return 0;
}

// ============================================================================================================
// ============================================================================================================
#else
int main( int argc, char* argv[] ) {


  //this 'main' is intended to be used repeatedly, reusing the graph
  //to create many remodelled files -- currently broken

  if( argc != 8 ) {
    cout << "Usage: " << argv[0]
        << "VoxFEFile.script "
           "ElementMapFile.elmap "
           "ElementGraphFile.graph "
           "LowerThreshold "
           "UpperThreshold "
           "DisplacementFile.txt "
           "Number-Of-Remodelling-Steps";
    cout << "\n\n";
    exit(-1);
  }

  char buff[32];
  string ParaBMU("mpirun -n 8 /home/richard/src/async/vox-fe/solver/PARA_BMU ");
  //string ParaBMU("/home/richard/src/async/vox-fe/solver/PARA_BMU ");
  string ParaBMUExec, dfile;

  ParaBMUExec = ParaBMU + argv[1];

  //ugly but it lets run and pauses until complete...
  int i = system( ParaBMUExec.c_str() );
  cout << "\n\n ** Running solver: " << ParaBMUExec
       << "\n\n";

  int step = 0;  //We consider the starting script to be "0"
  voxfeRemodelControl vrc( argv[1], argv[2], argv[3] );
  vrc.SetThresholds( atof( argv[4] ), atof( argv[5] ) );
  vrc.SetDisplacements( argv[6] );

  cout << "\n\n ** Running remodeller, step: " << step << "\tDisplacements: " << argv[6]
         << "\n\n";


  int nsteps = atoi( argv[7] );
  for( step=1; step<=nsteps; step++ ) {

    sprintf(buff, VOXFE_REMODEL_PRINT_SPEC, step);
    ParaBMUExec = ParaBMU + argv[1] + buff;

    int i = system( ParaBMUExec.c_str() );
    cout << "\n\n ** Running solver: " << ParaBMUExec
         << "\n\n";

    dfile = string(argv[6]) + string(buff);
    cout << "\n\n ** Running remodeller, step: " << step << "\tDisplacements: " << dfile
         << "\n\n";
    vrc.SetDisplacements( dfile.c_str() );

    //vrc.PrintGraph( cout, false );

  }

  return 0;
}
#endif



