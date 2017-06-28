#ifndef GUARD_PFESOLVER_H
#define GUARD_PFESOLVER_H

//======== C++ HEADER FILES =================


#include <cstdio>
#include <cstdlib>

#include <iostream>
#include <fstream>
#include <exception>
#include <cerrno>

#include <vector>
#include <set>
#include <map>
#include <string>

#include "Common.h"
#include "Commands.h"
#include "Cons.h"
#include "Node.h"
#include "Elem.h"
#include "Material.h"

#include "IntegrationMtx.h"
#include "PropertyMtx.h"
#include "GradientMtx.h"
#include "LocalStiffnessMtx.h"

#include "petscversion.h"


//=========== COMPARE FUNCTION for Elem/Node Sets/Maps  ==========

template <class T, class P>
struct compare {
    bool operator() ( T *n, T *m) const{
        P nx = n->x; P mx = m->x;
        P ny = n->y; P my = m->y;
        P nz = n->z; P mz = m->z;
        if(nz == mz){
            if(ny == my){
                if(nx == mx){
                    return false;
                }else{
                    return (nx < mx); }
            }else{
                return (ny < my); }
        }else{
            return (nz < mz); }
    }
};

//========== pFESolver class

class pFESolver
{
    
public:
    
    // MPI processor rank and number of processors
    
    MPI_Comm    comm; // change
    int         MPIrank;
    int         MPIsize;
    
    //======================== MAPS / ITERATORS ======================
    
    typedef        std::set< Elem< xyzType > * , compare< Elem< xyzType >, xyzType > > ElementSet; // set of elements
    typedef        std::set< Node< xyzType > * , compare< Node< xyzType >, xyzType > > NodeSet;    // set of nodes
    typedef        std::map< midxType, Material >  MaterialMap;   // map of materials
    typedef        std::vector<std::string>        VecString;     // vector of strings

    
    // specific sets/maps
    ElementSet      ElemS;  // set of valid elements;
    NodeSet         NodeS;  // set of valid nodes;
    MaterialMap     MateM;  // map of all materials
    
    
    // iterators
    typedef        ElementSet::iterator         ElementSet_it;
    typedef        ElementSet::const_iterator   ElementSet_const_it;
    typedef        NodeSet::iterator            NodeSet_it;
    typedef        NodeSet::const_iterator      NodeSet_const_it;
    typedef        MaterialMap::iterator        MaterialMap_it;
    typedef        MaterialMap::const_iterator  MaterialMap_const_it;
    typedef        VecString::iterator          VecString_it;

    
    //============================= VARS ============================
    
    double          A_SIZE, B_SIZE, C_SIZE, SCALE_FACTOR;
    xyzType         DX, DY, DZ, DTOT, DXY;  // elements per dim etc.
    xyzType         NX, NY, NZ, NTOT, NXY;  // nodes per dim etc.
    long int        MAXITER;                // max number of iterations
    double          TOLERANCE;              // tolerance
    // solver options
    char   SOLVERTYPE[64];
    char   PCTYPE[64];
    
    // set filenames for reference later
    char    SCRIPTFILE[MAX_FILENAME_LENGTH];
    char    MATERIALSFILE[MAX_FILENAME_LENGTH];
    char    MODELFILE[MAX_FILENAME_LENGTH];
    char    CONSTRAINTSFILE[MAX_FILENAME_LENGTH];
    char    OUTPUTFILE[MAX_FILENAME_LENGTH];

    //================= Flags
  
    bool VOXELSIZE_DONE;
    bool MATERIALS_DONE;
    bool VOXELDIMS_DONE;
    bool MODEL_DONE;
    bool CONSTRAINTS_DONE;
  
    bool SOLVE_DONE;
    bool FINISH_DONE;
    
    char data[MAX_DATA_LENGTH];
    
    //================== CONSTRUCTOR / DESTRUCTOR
    
    // COLLECTIVE
    pFESolver(const MPI_Comm comm, const int rank);
    
    // Destructor, frees anything that was allocated
    // COLLECTIVE
    ~pFESolver();
  
    // initialise variables etc.
    void Init();
    
    // check if basic setup complete. Cannot solve without this being true
    bool IsSetupDone();
    
    
    //================== Build model / setters
    
    // Sets the model dimensions in voxels
    bool SetVoxelDimensions(const xyzType dimx, const xyzType dimy, const xyzType dimz);
    
    // Sets the voxel dimensions : SetSize
    bool SetVoxelSize(const double a, const double b, const double c, const double sf );
    //const double a, const double b, const double c, const double sf);
    
    // Read material file and build material map
    bool LoadMaterials(const char *filename);
    
    // Read model file and create Elem sets and Node sets
    bool LoadModel(const char *filename);
    
    // Read constraints/forces file and add constraints.
    bool LoadConstraints(const char *filename);
    
    // Sets parameters (e.g. solver accuracy, other stuff)
    bool SetMaxIter(const long int themaxiter);
    
    bool SetTolerance(const double thetolerance);
  
    bool SetAlgorithmFEA(const char *thesolvertype, const char *thepctype);
    
    //================ Solve
    
    // Solves the system, updating the node displacements
    // COLLECTIVE
    PetscErrorCode Solve();
    
    //================ Output
    
    bool PrintDisplacements(const char *filename);
    
    //================= Add / Remove / Set
    
    // Adds a material with given properties, computes and stores the LSM
    bool AddMaterial(const midxType index, const double ym, const double pr, const int remodelflag);
    
    // Adds an element to the model (along with its nodes)
    bool AddElement(const midxType material, const xyzType elx, const xyzType ely, const xyzType elz);
    
    // Adds an element to the model (along with its nodes)
    bool AddConstraint(const xyzType x, const xyzType y, const xyzType z, const double vx, const double vy, const double vz, const int cx, const int cy, const int cz, const int pf);
    
    //=================
    
    // Delete node --> constrain deletion so you do not produce island nodes
    bool RemoveNode(const xyzType nx, const xyzType ny, const xyzType nz);
    
    // Delete an element to the model (along with its nodes)
    bool RemoveElement(const midxType material, const xyzType elx, const xyzType ely, const xyzType elz);
    
    // Delete constraint
    bool RemoveConstraint(const xyzType x, const xyzType y, const xyzType z);
    
    //================ Getters
    
    // Gets a set of all the (local) elements
    ElementSet GetLocalElements();
    
    // Gets a set of all the (local) nodes
    NodeSet GetLocalNodes();
    
    // Gets a particular element by coords
    Elem<xyzType>* GetElement(const xyzType x, const xyzType y, const xyzType z);
    
    // Gets a particular node by coords
    Node<xyzType>* GetNode(const xyzType x, const xyzType y, const xyzType z);
    
    Cons<xyzType>* GetConstraint(const xyzType x, const xyzType y, const xyzType z);
    
    
protected:
    // Petsc error code
    PetscErrorCode  ierr;

    GradientMtx gradientMtx[NODES_PER_ELEMENT];
    
    // PETSc version
    char petscVersionStr[MAX_FILENAME_LENGTH];
    
    // global and local row/col sizes for GSM/RHS etc.
    PetscInt       globalgsmrows, globalgsmcols, localgsmrows, localgsmcols, diagonalterms, offdiagonalterms, localrhslength, localsollength;
    
    // force vector : rhs
    PetscScalar    *forcevalue;
    PetscInt       *forcecolidx;
    idxType        forcecount, totalrhs;
    
    // Petsc KSP solver
    KSP             ksp; // the solver context
    
    // vector solution etc.
    PetscScalar *solution; // solution data retrieved
    Vec vecout;
    
    //================ Build model
    
    bool RenumberNodes();
    
    //================ Matrices
    
    // COLLECTIVE
    PetscErrorCode ComputeGSM(Mat *gsm);
    
    // COLLECTIVE
    PetscErrorCode ComputeRHS(Vec *rhs);
    
    //================ Allocate Local Rows / Columns (matrix / vec)
    
    PetscErrorCode AllocateLocalMatrix(Mat *mat);
    
    PetscErrorCode AllocateLocalVec(Vec *vec);
    
    //================ Finish up
    
    PetscErrorCode cleanup();
    
    //================ Micro methods (inline)
    
    // MicroGetters
    
    void GetNodeCons(NodeSet_const_it cnitr, Cons<xyzType>* cons)
    {
        cons = ((*cnitr)->constraint);
    }
    
    Cons<xyzType>* GetNodeCons(NodeSet_const_it cnitr)
    {
        return (*cnitr)->constraint;
    }
    
    bool GetNodeConsX(NodeSet_const_it cnitr)
    {
        if( (*cnitr)->constraint )
            return (* ((*cnitr)->constraint) ).cx;
        else
            return true;
    }
    
    bool GetNodeConsY(NodeSet_const_it cnitr)
    {
        if( (*cnitr)->constraint )
            return (* ((*cnitr)->constraint) ).cx;
        else
            return true;
    }
    
    bool GetNodeConsZ(NodeSet_const_it cnitr)
    {
        if( (*cnitr)->constraint )
            return (* ((*cnitr)->constraint) ).cx;
        else
            return true;
    }
    
    idxType GetNodeIndex(NodeSet_const_it cnitr)
    {
        return (*cnitr)->idx;
    }
    
    
    Elem<xyzType>* GetNodeElement(NodeSet_const_it cnitr, int el)
    {
        return (*cnitr)->elems[el];
    }
    
    midxType GetElementMaterial(NodeSet_const_it cnitr, int el)
    {
        return (*( (*cnitr)->elems[el]) ).material;
    }
    
    idxType GetNodeNeighbourIndex(NodeSet_const_it cnitr, int el, int neighbour)
    {
        return ( (*((*cnitr)->elems[el])).nodes[neighbour] )->idx;
    }
    
    bool GetNodeNeighbourConsX(NodeSet_const_it cnitr, int el, int neighbour)
    {
        return (*(( (*((*cnitr)->elems[el])).nodes[neighbour] )->constraint) ).cx;
    }
    
    bool GetNodeNeighbourConsY(NodeSet_const_it cnitr, int el, int neighbour)
    {
        return (*(( (*((*cnitr)->elems[el])).nodes[neighbour] )->constraint) ).cy;
    }
    
    bool GetNodeNeighbourConsZ(NodeSet_const_it cnitr, int el, int neighbour)
    {
        return (*(( (*((*cnitr)->elems[el])).nodes[neighbour] )->constraint) ).cz;
    }
    
    Cons<xyzType>* GetNodeNeighbourCons(NodeSet_const_it cnitr, int el, int neighbour)
    {
        return  ((* ((*cnitr)->elems[el])).nodes[neighbour])->constraint;
    }
    
    double GetLSMValue(midxType midx, idxType index)
    {
        return (MateM[midx].lsm.getval(index));
    }
    
    
};


#endif
