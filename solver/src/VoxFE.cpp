//
//  VoxFE.cpp
//
//
//  Created by N. Banglawala
//
//

#include "VoxFE.h"


//========== CONSTRUCTOR / DESTRUCTOR

VoxFE::VoxFE(const MPI_Comm mpicomm, const int rank) : comm(mpicomm), MPIrank(rank), solver(mpicomm, rank)
{}

VoxFE::~VoxFE()
{}

// script command
void VoxFE::create_command_map()
{
    
    // **** Must match enums in Commands.h ****
    // MUST have all options set before LOAD CONSTRAINTS called
    
    CommandM["SET_VOXEL_SIZE"] = CMD_SET_VOXEL_SIZE;
    CommandM["VOXEL_SIZE"] = CMD_SET_VOXEL_SIZE;
    
    CommandM["LOAD_MATERIALS"] = CMD_LOAD_MATERIALS;
    CommandM["LOAD_MATERIALS_FILE"] = CMD_LOAD_MATERIALS;
    
    CommandM["LOAD_MODEL"] = CMD_LOAD_MODEL;
    CommandM["LOAD_MCTSCAN"] = CMD_LOAD_MODEL;
    
    CommandM["SET_TOLERANCE"] = CMD_SET_TOLERANCE;
    CommandM["TOLERANCE"] = CMD_SET_TOLERANCE;
    
    CommandM["SET_MAX_ITER"] = CMD_SET_MAX_ITER;
    CommandM["MAX_ITER"] = CMD_SET_MAX_ITER;
    
    CommandM["SET_ALGORITHM_FEA"] = CMD_SET_ALGORITHM_FEA;
    CommandM["ALGORITHM_FEA"] = CMD_SET_ALGORITHM_FEA;
  
    // last command before solve is prompted
    CommandM["SELECTION_OF_NODES"] = CMD_SELECTION_OF_NODES; // DOES NOTHING !!
    CommandM["LOAD_CONSTRAINTS"] = CMD_LOAD_CONSTRAINTS;
    CommandM["SELECT_NODE_FILE"] = CMD_LOAD_CONSTRAINTS;
    
    CommandM["SELECT_NODE_3D"] = CMD_SELECT_NODE_3D;
    CommandM["PRESERVE_NODE_3D"] = CMD_PRESERVE_NODE_3D;
    
    CommandM["COMPUTE_SED"] = CMD_COMPUTE_SED; // DOES NOTHING !!
    
    CommandM["SOLVE"] = CMD_SOLVE; // not specified in script
    
    CommandM["PRINT_DISPLACEMENTS"] = CMD_PRINT_DISPLACEMENTS;
    CommandM["PRINT_X"] = CMD_PRINT_DISPLACEMENTS;
    
    CommandM["FINISH"] = CMD_FINISH; // not specified in script
}// create_command_map()


// FIX ME : add some more checking...
// break input string into command and data
void VoxFE::parse_input_line(std::string line)
{
    int i(0);
    int index(0);
    
    size_t cstart = line.find_first_not_of(" \t\n\r");

    
    // only add line if non-white space characters found
    if( cstart != std::string::npos )
    {
        size_t cend = line.find_first_of(" \t\n\r", cstart);
        std::string commandstr = line.substr(cstart, cend-cstart);
       
        // search for command in command map
        if( (CommandM.find(commandstr) != CommandM.end()) )
            index = CommandM[ commandstr ];
            else
        {
            printf("%d: ERROR: command %s does not exist, please check script file\n", MPIrank, commandstr.c_str());
            return;
        }
        
        if(cend < line.length())
            InputM[ index ] = line.substr(cend, line.length());
        else
            InputM[ index ] = "";
        
        printf("%d: COMMAND %d:  %s,  DATA: %s \n", MPIrank, index, commandstr.c_str(), (InputM[index]).c_str() );
    }// if non-empty line
    
}// parse_input_line()


// MASTER ONLY
void VoxFE::output_end_flag(bool allOK)
{
    char filename[] = "FECompleted.txt";
    
    FILE * finalFile = fopen(filename,"w");
    
    if(allOK)
        fprintf(finalFile, "solver_done\n");
    else
        fprintf(finalFile, "solver_failed\n");
    
    fclose(finalFile);
}

// MASTER ONLY
void VoxFE::output_input_flag()
{
    char filename[] = "InputCompleted.txt";
    
    FILE * doneFile = fopen(filename,"w");
    
    fprintf(doneFile, "solver_files_input\n");
    
    fclose(doneFile);
}

//========== execute commands

void VoxFE::read_script_file(const char *filename)
{
    MPI_Barrier(comm);
    startT = MPI_Wtime();
    
    bool          OK(true);
    FILE *        scriptFile;
    
    MPI_Comm_size(comm, &MPIsize);
    
    // SET SCRIPTFILE
    strcpy(solver.SCRIPTFILE, filename);
    
    if( MPIrank == MASTER )
        printf("%d: MPIsize = %zd \n", MPIrank, MPIsize);
    
    int command(0);
    int totalscript(0);
    InputM_it itr;
    
    // create command map
    create_command_map();
    
    // read commands from script into Input map
    if(MPIrank == MASTER)
    {
        // open script file
        std::ifstream scriptF(filename);
        if(scriptF)
        {
            try
            {
                printf("%d: about to read from the script file : %s \n", MPIrank,filename);
                int i(0);
                
                // for each line in script
                for( std::string line; getline(scriptF, line); )
                    parse_input_line(line);

                // add finish command
                InputM[ CMD_FINISH ] = " ";
                
                totalscript = InputM.size();
                itr = InputM.begin();
            }// try this
            catch( std::exception &e )
            {
                printf("%d: ERROR: could not read script file \n", MPIrank);
                std::cerr << "exception caught: " << e.what() << std::endl ;
                OK = false;
                
                output_end_flag(OK);
                
                return;
            }// catch this
            scriptF.close();
        }// if file exists
        else
        {
            printf("%d: ERROR: could not find script file \n", MPIrank);
            return;
        }
        
        itr = InputM.begin();
        
    }// if MASTER
    
    MPI_Bcast( &totalscript, 1, MPI_INT, MASTER, comm );
    
    if( MPIrank == MASTER )
        printf("%d: script read, totalinput = %d \n", MPIrank, totalscript);
    
    int in(0);
    int count(0);
    
    
    // execute commands
    while( OK && (command != CMD_FINISH) )
    {// while OK and not finished
        
        // check if setup done, else read from script input
        if(MPIrank == MASTER)
        {// get command
            
            if( solver.IsSetupDone() && !solver.SOLVE_DONE )
            {
                printf("%d: set up is done, about to solve \n", MPIrank);
                command = CMD_SOLVE;
                // no data string for SOLVE;
                
                //output input file flag
                output_input_flag();
            }
            else
            {
                command = (itr->first);
                strcpy( data, (itr->second).c_str() );
                
                ++itr;
            }
        }// if MASTER
        
        MPI_Bcast(&command, 1, MPI_INT, MASTER, comm);
        
        OK = execute_command(command);
        if( MPIrank == MASTER )
           printf("%d: command : %d executed and Ok = %d\n", MPIrank, command, OK);
        
    }// for each input line
    
    if( !OK ) // if something went wrong with setup
        printf("%d: ERROR: could not continue reading script, something went wrong \n", MPIrank);
    else
    {
        if( MPIrank == MASTER )
            printf("%d: all OK \n", MPIrank);
    }
    
    endT = MPI_Wtime();
    
    if( MPIrank == MASTER )
    {
        output_end_flag(OK);
        printf("%d: total time = %le \n", MPIrank, endT-startT);
    }
    
}// read_script_file();


// execute commands here
bool VoxFE::execute_command(int command)
{
    //printf("%d: about to execute with command = %d \n", MPIrank, command);
    bool OK(false);
    switch ( command )
    {
        case CMD_SET_VOXEL_SIZE : // SET VOXEL SIZE, SCALE FACTOR
        {
            double vsize[4]; // {a, b, c, sf};
            if( MPIrank == MASTER )
                OK = ( sscanf( data, " %lf %lf %lf %lf \n", &vsize[0], &vsize[1], &vsize[2], &vsize[3]) == 4 );
            MPI_Bcast(&OK, 1, MPI_INT, MASTER, comm);
            if( OK )
            {
                MPI_Bcast( vsize, 4, MPI_DOUBLE, MASTER, comm );
                OK = solver.SetVoxelSize( vsize[0], vsize[1], vsize[2], vsize[3] );
            }
            break;
        }
        case CMD_LOAD_MATERIALS : // LOAD MATERIALS FROM FILE
        {
            char materialsfn[MAX_FILENAME_LENGTH];
            int materialsfnlen(0);
            if( MPIrank == MASTER )
            {
                OK = ( sscanf( data, " %s \n", &materialsfn ) == 1 );
                materialsfnlen = strlen(materialsfn) + 1;
            }
            MPI_Bcast(&OK, 1, MPI_INT, MASTER, comm);
            if( OK )
            {
                MPI_Bcast( &materialsfnlen, 1, MPI_INT, MASTER, comm);
                MPI_Bcast( materialsfn, materialsfnlen, MPI_CHAR, MASTER, comm );
                OK = solver.LoadMaterials( materialsfn );
            }
            break;
        }
        case CMD_LOAD_MODEL : // SET VOXEL DIMENSIONS, LOAD MODEL
        {
            xyzType vdims[3];
            char modelfn[MAX_FILENAME_LENGTH];
            int modelfnlen(0);
            if( MPIrank == MASTER )
            {
                OK = ( sscanf( data, " %d %d %d %s \n", &vdims[0], &vdims[1], &vdims[2], &modelfn ) == 4 );
                modelfnlen = strlen(modelfn) + 1;
            }
            MPI_Bcast(&OK, 1, MPI_INT, MASTER, comm);
            if( OK )
            {
                MPI_Bcast( vdims, 3, MPI_DOUBLE, MASTER, comm );
                OK = solver.SetVoxelDimensions( vdims[0], vdims[1], vdims[2] );
            }
            if( OK )
            {
                MPI_Bcast( &modelfnlen, 1, MPI_INT, MASTER, comm);
                MPI_Bcast( modelfn, modelfnlen, MPI_CHAR, MASTER, comm );
                OK = solver.LoadModel( modelfn );
            }
            break;
        }
        case CMD_SELECTION_OF_NODES : // SELECTION OF NODES -- DOES NOTHING...
        {
            OK = true;
            break;
        }
            
        case CMD_LOAD_CONSTRAINTS : // LOAD CONSTRAINTS (FORCES, PRESERVED NODES)
        {
            char constraintsfn[MAX_FILENAME_LENGTH];
            int constraintsfnlen(0);
            if( MPIrank == MASTER )
            {
                OK = ( sscanf( data, " %s \n", &constraintsfn ) == 1 );
                constraintsfnlen = strlen(constraintsfn)+1;
            }
            MPI_Bcast(&OK, 1, MPI_INT, MASTER, comm);
            if( OK )
            {
                MPI_Bcast(&constraintsfnlen, 1, MPI_INT, MASTER, comm);
                MPI_Bcast( constraintsfn, constraintsfnlen, MPI_CHAR, MASTER, comm );
                OK = solver.LoadConstraints( constraintsfn );
            }
            break;
        }
            
        case CMD_SET_MAX_ITER : // MAX ITERATION
        {
            long int themaxiter;
            if( MPIrank == MASTER )
                OK = ( sscanf( data, " %d \n", &themaxiter ) == 1 );
            MPI_Bcast(&OK, 1, MPI_INT, MASTER, comm);
            if( OK )
            {
                MPI_Bcast( &themaxiter, 1, MPI_LONG_INT, MASTER, comm );
                OK = solver.SetMaxIter( themaxiter );
            }
            break;
        }
            
        case CMD_SET_TOLERANCE : // TOLERANCE
        {
            double thetolerance;
            if( MPIrank == MASTER )
                OK = ( sscanf( data, " %lf \n", &thetolerance ) == 1 );
            MPI_Bcast(&OK, 1, MPI_INT, MASTER, comm);
            if( OK )
            {
                MPI_Bcast( &thetolerance, 1, MPI_DOUBLE, MASTER, comm );
                OK =  solver.SetTolerance( thetolerance );
            }
            break;
        }
            
        // Doesn't do anything at present. check again list of  solver/pc types ?
        case CMD_SET_ALGORITHM_FEA : // SET SOLVERTYPE & PCTYPE
        {
            char thesolvertype[MAX_FILENAME_LENGTH];
            char thepctype[MAX_FILENAME_LENGTH];
            int solvertypelen;
            int pctypelen;
            int nargs; // need this to accommodate old solver
            if( MPIrank == MASTER )
            {
                nargs = sscanf( data, " %s  %s \n", &thesolvertype, &thepctype );
                
                if( nargs == 1 )
                {
                    strcpy(thepctype, "");
                    OK = true;
                }
                if( nargs == 2 )
                    OK = true;
                // if nargs 1 or 2, else default is OK = false
                
                solvertypelen = strlen(thesolvertype);
                pctypelen = strlen(thepctype);
                
            }// if master
            MPI_Bcast(&OK, 1, MPI_INT, MASTER, comm);
            if( OK )
            {
                MPI_Bcast( &solvertypelen, 1, MPI_INT, MASTER, comm );
                MPI_Bcast( &pctypelen, 1, MPI_INT, MASTER, comm );
                MPI_Bcast( thesolvertype, solvertypelen, MPI_CHAR, MASTER, comm );
                MPI_Bcast( thepctype, pctypelen, MPI_CHAR, MASTER, comm );
                OK = solver.SetAlgorithmFEA( thesolvertype, thepctype );
            }
            break;
        }
            
        case CMD_COMPUTE_SED : // COMPUTE SED -- DOES NOTHING...
        {
            OK = true;
            break;
        }
            
        case CMD_SOLVE : // SOLVE
        {
            OK = ( !solver.Solve() ); // PetscErrorCode
            break;
        }
            
        case CMD_PRINT_DISPLACEMENTS : // PRINT DISPLACEMENTS
        {
            // only MASTER process prints displacements (for now)
            if( MPIrank == MASTER )
            {
                if( solver.SOLVE_DONE )
                {
                    char outputfn[MAX_FILENAME_LENGTH];
                    OK = ( sscanf( data, " %s \n", &outputfn ) == 1 );
                    if( OK )
                        OK = solver.PrintDisplacements( outputfn );
               
                }// if solve done
                else
                {
                    printf( "%d: ERROR: cannot print displacements before solving! \n", MPIrank);
                    OK = false;
                }
            }// if MASTER rank
            
            MPI_Bcast( &OK, 1, MPI_INT, MASTER, comm );
            break;
        }
            
        case CMD_FINISH : // FINISHED READING SCRIPT
        {
            if( MPIrank == MASTER )
                printf( "%d: Finished reading script \n", MPIrank );
            // cleanup here
            solver.FINISH_DONE = true;
            OK =  true;
            
            break;
        }
        default:
            printf( "%d: ERROR: command %d not recognised\n", MPIrank, command );
    }
    
    return OK;
}// execute_command()







