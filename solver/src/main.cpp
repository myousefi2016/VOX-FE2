#include "VoxFE.h"

using namespace std;

int main(int argc, char* argv[])
{
    
    //======================= MPI / PETSc =======================
    
    // Start up
    PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
    
    int rank;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);

    
    VoxFE *vox;
    vox = new VoxFE(comm, rank);
    
    if( rank == 0 )
        printf("%d: calling read_script_file \n", rank);
    vox -> read_script_file(argv[1]);
    
    // Once script has been executed, clean up and shut down
    
    delete vox;
    vox = 0;
    
   
    PetscFinalize();
    
    return 0;
}
