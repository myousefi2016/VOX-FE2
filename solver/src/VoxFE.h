//
//  VoxFE.h
//  
//
//  Created  by N. Banglawala
//
//

#ifndef GUARDS_VOXFE_H
#define GUARDS_VOXFE_H

#include "pFESolver.h"

class VoxFE{

public:
    
    MPI_Comm    comm;
    int         MPIrank;
    int         MPIsize;
    
    char data[MAX_DATA_LENGTH];
    
    pFESolver   solver;
    
    typedef        std::map<std::string, int> MapStringInt;  // map of commands
    typedef        std::map<int, std::string> MapIntString;  // map of commands
    
    MapStringInt   CommandM; // mapping commands to long ints
    MapIntString   InputM; // mapping commands to long ints
    
    typedef        MapStringInt::iterator           CommandM_it;
    typedef        MapIntString::iterator             InputM_it;
    
    // measure time
    double startT, endT;
    
    VoxFE(const MPI_Comm mpicomm, const int mpirank);
    ~VoxFE();

    // create map < command string , command enum>
    void create_command_map();
    
    // MASTER ONLY
    void read_script_file(const char * filename);
    
    // parse line from script --> command, data
    void parse_input_line(std::string line);
    
    // COLLECTIVE
    bool execute_command(int command);
    
    // MASTER only
    // output final flag file to show solver has completed (successfully or otherwise)
    void output_end_flag(bool allOK);
    
    // MASTER only
    // output flag file once input has been read in
    void output_input_flag();
};



#endif /* defined( GUARDS_VOXFE_H ) */
