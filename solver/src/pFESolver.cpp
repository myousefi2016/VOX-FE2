// VoxFE new PETSc solver class

#include "pFESolver.h"

//============ Constructor / Destructor

// COLLECTIVE
pFESolver::pFESolver(const MPI_Comm mpicomm, const int rank) : comm(mpicomm), MPIrank(rank)
{
    // set all flags to false
    Init();
}

// COLLECTIVE
pFESolver::~pFESolver()
{}

//============ Initialise

// COLLECTIVE
// initialise variables here
// part of solver, not script
void pFESolver::Init()
{
    // setup and rebuild flags
    VOXELSIZE_DONE = false;
    VOXELDIMS_DONE = false;
    MATERIALS_DONE = false;
    MODEL_DONE = false;
    CONSTRAINTS_DONE = false;
    
    SOLVE_DONE = false;
    FINISH_DONE = false;
    
    MPI_Comm_size(comm, &MPIsize);
    
}// Init()


//========== Check setup done

// check all minimal setup commands have been called. Basic sanity check before solve called
bool pFESolver::IsSetupDone()
{
    bool OK = VOXELSIZE_DONE && VOXELDIMS_DONE && MATERIALS_DONE && MODEL_DONE && CONSTRAINTS_DONE;
    return OK;
}// IsSetupDone()

//========= Setters

// Sets the voxel dimensions and scale factor
bool pFESolver::SetVoxelSize(const double a, const double b, const double c, const double sf )
{
    bool OK(false);
    if(a * b * c * sf)
    { // check if not zero
        
        A_SIZE = a;
        B_SIZE = b;
        C_SIZE = c;
        SCALE_FACTOR = sf;
        OK = true;
        VOXELSIZE_DONE = true;
    }
    else
    {
        printf("%d: ERROR: could not add voxelsizes and/or scale factor \n",MPIrank);
        OK = false;
    }
    return OK;
}// SetVoxelSize()


// Sets the model dimensions in voxels : no default
bool pFESolver::SetVoxelDimensions(const xyzType dimx, const xyzType dimy, const xyzType dimz)
{
    bool OK(true);
    if( dimx*dimy*dimz )
    {
        DX = dimx;
        DY = dimy;
        DZ = dimz;
        DTOT = DX*DY*DZ;
        DXY  = DX*DY;
        NX   = DX + 1;
        NY   = DY + 1;
        NZ   = DZ + 1;
        NTOT = NX*NY*NZ;
        NXY  = NX*NY;
        OK   = true;
        VOXELDIMS_DONE = true;
    }
    else
    {
        printf("%d: ERROR: voxel dimensions must be non-zero\n", MPIrank);
        OK = false;
    }
    return  OK;
}// SetVoxelDimensions()


bool pFESolver::SetMaxIter(const long int themaxiter)
{
    bool OK(true);
    MAXITER = themaxiter;
    return OK;
}// SetMaxIter()

// set tolerance
bool pFESolver::SetTolerance(const double thetolerance)
{
    bool OK(true);
    TOLERANCE = thetolerance;
    return OK;
}// SetTolerance()

bool pFESolver::SetAlgorithmFEA(const char *thesolvertype, const char *thepctype)
{
    bool OK(true);
    strcpy( SOLVERTYPE, thesolvertype );
    strcpy( PCTYPE, thepctype );
    return OK;
}

//================ Loaders

// Create Material map
bool pFESolver::LoadMaterials(const char *filename)
{
    strcpy( MATERIALSFILE, filename );
    if( MPIrank == MASTER )
        printf("%d: material file is %s \n", MPIrank, filename);
    
    FILE * materialFile;
    
    midxType material_index;
    double youngs_mod;
    double poissons_rat;
    int remodel_flag;
    bool OK(true);
    
    if(VOXELSIZE_DONE)
    {
        if(MPIrank == MASTER)
            printf("%d: Creating Material Map \n",MPIrank);
        
        // compute gradient matrices : needed for computing LSM of each material
        for(int n=0; n < NODES_PER_ELEMENT; ++n){
    
            gradientMtx[n] = GradientMtx(A_SIZE, B_SIZE, C_SIZE, n);
        }
        
        if( MPIrank == MASTER )
            printf("%d: Opening material file \n", MPIrank);
        materialFile = fopen(filename, "r");
        errno = 0;
        if(materialFile == NULL)
        {
            printf("%d: ERROR: could not open material file \n", MPIrank);
            OK = false;
            return OK;
        }else
        {
            if( MPIrank == MASTER )
                printf("%d: material file open \n", MPIrank);
            
            while( fscanf(materialFile, "%d %lf %lf %d \n", &material_index, &youngs_mod, &poissons_rat, &remodel_flag) != EOF && OK )
            {// while
                
                OK = AddMaterial(material_index, youngs_mod, poissons_rat, remodel_flag);
                
            }// while adding materials
            if(!OK){
                printf("ERROR: couldn't add material \n");
                fclose(materialFile);
                OK = false;
                return OK;
            }// if material not added
            
            fclose(materialFile);
            
            MATERIALS_DONE = true;
        }
        
    }// if voxel sizes specified
    else
    {
        printf("%d: ERROR: voxel sizes must be specified before loading materials \n", MPIrank);
        OK = false;
    }// voxel size or voxel dimensions not specified
    
    
    return OK;
}// LoadMaterials()


// Load elements, create element and node sets
bool pFESolver::LoadModel(const char *filename)
{
    strcpy( MODELFILE, filename );
    
    midxType                    material;
    xyzType                     element_num, x, y, z;
    FILE *                      modelFile;
    bool OK(true);
    
    OK = VOXELSIZE_DONE && VOXELDIMS_DONE && MATERIALS_DONE;
    if(!OK)
    {
        printf("%d: ERROR: voxel size, dimensions and materials must be specified before loading model \n", MPIrank);
        return OK;
    }
    
    if(MPIrank==MASTER)
        printf("%d: Creating Element and Node Sets\n",MPIrank);
    
    modelFile = fopen(filename,"r");
    
    if(modelFile == NULL)
    {
        printf ("%d: ERROR: could not open model file \n", MPIrank);
        OK = false;
        return OK;
    }
    else
    {
        // discard first two lines and continue
        fscanf (modelFile, " ");
        fscanf (modelFile, "%d", &element_num);
        
        while( fscanf(modelFile, "%d %d %d %d %d",&element_num, &material, &x, &y, &z)  != EOF && OK)
        {
            // script coords are 1-based, solver coords are 0-based
            --x; --y; --z;
            
	    OK = AddElement(material, x, y, z);
        }
        if(!OK)
        {
            return OK;
        }
    }// else if file exists
    fclose (modelFile);
    
    OK = RenumberNodes();
    if(!OK)
    {
        return OK;
    }else
    {
        MODEL_DONE = true;
    }
    
    return OK;
}// LoadModel()


// Read in constraint and force data
bool pFESolver::LoadConstraints(const char *filename){
    
    strcpy( CONSTRAINTSFILE, filename );
    
    bool OK(true);
    xyzType x, y, z;
    int cx, cy, cz;
    double vx, vy, vz;
    char str[80];
    
    VecString lines;
    VecString_it vit = lines.begin();
    NodeSet_it nodeptr = NodeS.begin();
    
    idxType totlines;
    idxType conscount(0);
    
    
    OK = ( VOXELSIZE_DONE && VOXELDIMS_DONE && MATERIALS_DONE && MODEL_DONE );
    if(!OK)
    {
        printf("%d: ERROR: model parameters, materials and model need to be specified before constraints \n", MPIrank);
        return OK;
    }// if voxel sizes etc. specified
    
    
    std::ifstream consFile(filename);
    if(MPIrank == MASTER)
        printf("%d: Reading constraints\n", MPIrank);
    try
    {
        size_t cstart;
        for( std::string line; getline(consFile, line); )
        {
            // check if empty line
            cstart = line.find_first_not_of(" \t\n\r");
            
            // only add line if non-white space characters found
            if( cstart != std::string::npos )  lines.push_back(line);
        }
    }// try
    catch(std::exception &e)
    {
        printf("%d: ERROR: could not open constraint file \n", MPIrank);
        OK = false;
        return OK;
    }// catch
    consFile.close();
    totlines = lines.size();

    
    Node<xyzType> * tmpnode = new Node<xyzType>();
    
    for(vit = lines.begin(); vit != lines.end(); ++vit)
    {
        sscanf ((*vit).c_str(), "%s %d %d %d %lf %lf %lf %d %d %d",&str, &x, &y, &z, &vx, &vy, &vz, &cx, &cy, &cz);
        
        if(cx || cy || cz)
        {
            --x; --y; --z;
   
            cx = !cx;  cy = !cy;  cz = !cz;

            // add constraint
            ++conscount;
            if( strcmp(str, "CMD_SELECT_NODE_3D") )

                OK = AddConstraint(x, y, z, vx, vy, vz, cx, cy, cz, 0);
            else if( strcmp(str, "CMD_PRESERVE_NODE_3D") )
        
                OK = AddConstraint(x, y, z, vx, vy, vz, cx, cy, cz, 1);
            else
            {
                printf("%d: ERROR: constraint command not recognised : %s \n", MPIrank, str);
                OK = false;
                return OK;
            }
            
            if(!OK)
            {
                printf("%d: ERROR: could not add constrained node %d : [%d, %d, %d] \n", MPIrank, conscount, x, y, z);
                return OK;
            }
        }else
            break;
    }// for line in lines
    
    totalrhs = DOF_3D*(totlines-conscount);
    forcevalue = new PetscScalar[totalrhs];
    forcecolidx = new PetscInt[totalrhs];
    
    conscount = 0; // not important here
    for(vit; vit != lines.end(); ++vit){// is a force
        sscanf ((*vit).c_str(), "%s %d %d %d %lf %lf %lf %d %d %d",&str, &x, &y, &z, &vx, &vy, &vz, &cx, &cy, &cz);
  
        --x; --y; --z;
        cx = !cx;  cy = !cy;  cz = !cz;
        
        // add constraint;
        ++conscount;
        
        if( strcmp(str, "CMD_SELECT_NODE_3D") )
       
            OK = AddConstraint(x, y, z, vx, vy, vz, cx, cy, cz, 0);
        else if( strcmp(str, "CMD_PRESERVE_NODE_3D") )
          
            OK = AddConstraint(x, y, z, vx, vy, vz, cx, cy, cz, 1);
        else
        {
            printf("%d: ERROR: constraint command not recognised : %s \n",MPIrank, str);
            OK = false;
            return OK;
        }
        if(!OK)
        {
            printf("%d: ERROR: could not add constrained node %d : [%d, %d, %d] \n", MPIrank, conscount, x, y, z);
            return OK;
        }
        
    }// for lines
    
    CONSTRAINTS_DONE = true;
    
    return OK;
}// LoadConstraints()

//==================== Output

bool pFESolver::PrintDisplacements(const char *filename)
{
    // read displacement from a given node. So cycle through nodelist and get displacement and print
    // get displacements from node and output -- currently fixed format
    
    strcpy( OUTPUTFILE, filename );
    
    bool OK(false);
    idxType globalID;
    //idxType localID;
    
    xyzType nx, ny, nz;
    PetscScalar dx, dy, dz;
    FILE * outFile;
    NodeSet_it nitr;
    
    try
    {
        outFile = fopen(filename,"w");
        if(outFile != NULL){
            
            printf("%d: output file is open for displacements\n", MPIrank);
        }else{
            printf("ERROR: output file NULL\n");
            return OK;
        }
        
        char materialFilename[] = "./axial_material.txt";
        fprintf(outFile,"nNodes: %d\n", NodeS.size());
        fprintf(outFile,"Nodes: %d  %d  %d\n", NX, NY, NZ);
        fprintf(outFile,"Materials: %s\n", MATERIALSFILE);
        
        for (nitr = NodeS.begin(); nitr != NodeS.end(); ++nitr)
        {
            nx = (*nitr)->x;
            ny = (*nitr)->y;
            nz = (*nitr)->z;
            dx = (*nitr)->dx;
            dy = (*nitr)->dy;
            dz = (*nitr)->dz;
	    globalID = nx + ny*NX + nz*NXY;      
      
	    fprintf(outFile, "%d   %d   %d   %d   %.9le   %.9le   %.9le  \n", globalID, nx, ny, nz, dx, dy, dz);
        }// for each node in Node Set
        fclose(outFile);

        OK = true;
    }
    catch(std::exception &e)
    {
        printf("ERROR: could not print displacements \n");
        fclose(outFile);
        OK = false;
        return OK;
    }
    return ( OK );
    
} // PrintDisplacements()


//==================== Adders

// Adds a material with given properties, computes and stores the LSM
bool pFESolver::AddMaterial(const midxType index, const double ym, const double pr, const int rf){
    bool OK(false);
    try
    {
        LocalStiffnessMtx lsm(A_SIZE, B_SIZE, C_SIZE, ym, pr, gradientMtx);
         
        Material * new_material = new Material(index, ym, pr, lsm, rf);
        
        MateM[index] = (*new_material);
        
        delete new_material;
        new_material = 0;
        
        OK = true;
    }
    catch (std::exception &e)
    {
        printf("ERROR: could not add material : %d %lf %lf %d\n", index, ym, pr, rf);
        OK = false;
    }
    return OK;
}


// Add element and associated nodes if they don't already exist
bool pFESolver::AddElement(const midxType material, const xyzType x, const xyzType y, const xyzType z){
    std::pair< ElementSet_it , bool>  elem_pair;
    std::pair< NodeSet_it    , bool>  node_pair;
    xyzType  nx, ny, nz;
    
    NodeSet_it                   nitr;
    ElementSet_it                eitr;
    
    idxType                      ecount(0), ncount(0);
    bool OK(true);
    
    // if material is in the MaterialMap i.e. an accepted material type
    try
    {
        // does material exist ?
        if( !(MateM.find(material) != MateM.end()) )
        {
            printf("ERROR: material %d does not exist \n",material);
            OK = false;
            return OK;
        }

        Elem<xyzType> * elem = new Elem<xyzType>();
        elem->x = x; elem->y = y; elem->z = z;
        elem->material = material;
        elem->idx = x + y*DX + z*DXY;
        ++ecount;
        
        elem_pair = ElemS.insert(elem);
        eitr = elem_pair.first;
        
        Node<xyzType>* tmpnode = new Node<xyzType>();
        
        for(int n=0; n<NODES_PER_ELEMENT; ++n){// for each node
            
            //calculate nodal coordinates/idx from element coordinates
            nx = x + (n%2);
            ny = y + (n/2)*(n<4) + (n/6)*(n>4);
            nz = z + (n/4);
            tmpnode->x = nx;
            tmpnode->y = ny;
            tmpnode->z = nz;
            
            nitr = NodeS.find(tmpnode); // try to find node in Node Map
            
            if(nitr == NodeS.end()){ // if node not in Node Map

                Node<xyzType> * newnode = new Node<xyzType>();
                newnode->x = nx;
                newnode->y = ny;
                newnode->z = nz;
                newnode->idx = nx + ny*NX + nz*NXY;
    
                ++ncount;
                node_pair = NodeS.insert(newnode);
                nitr = node_pair.first;
                
            }// if node exists
            
            (*eitr)->nodes[n] = (*nitr);
            (*nitr)->elems[n] = (*eitr);
            
        }// for each node
        
        delete tmpnode;
        
    }// try
    catch(std::exception &e){
        printf("ERROR: could not add element \n");
        OK = false;
    }
    return OK;
}// AddElement()


// Adds constraint to the node
bool pFESolver::AddConstraint(const xyzType x, const xyzType y, const xyzType z, const double vx, const double vy, const double vz, const int cx, const int cy, const int cz, const int pf)
{
    
    NodeSet_it nitr = NodeS.begin();
    idxType nidx;
    Node<xyzType> * tmpnode = new Node<xyzType>();
    Cons<xyzType> * cons = new Cons<xyzType>();
    bool OK(true);
    
    cons->cx = cx;
    cons->cy = cy;
    cons->cz = cz;
    cons->vx = vx;
    cons->vy = vy;
    cons->vz = vz;
    cons->preserve = pf;
    
    tmpnode->x = x; tmpnode->y = y; tmpnode->z = z;

    nitr = NodeS.find(tmpnode);

    
    if(nitr != NodeS.end())
    {
        cons->node = (*nitr);
        (*nitr)->constraint = cons;
        
        if(cx && cy && cz) // add force
        {
            nidx = (*nitr)->idx;
            
            forcevalue[forcecount*DOF_3D + 0] = vx;
            forcevalue[forcecount*DOF_3D + 1] = vy;
            forcevalue[forcecount*DOF_3D + 2] = vz;
            forcecolidx[forcecount*DOF_3D + 0] = nidx*DOF_3D + 0;
            forcecolidx[forcecount*DOF_3D + 1] = nidx*DOF_3D + 1;
            forcecolidx[forcecount*DOF_3D + 2] = nidx*DOF_3D + 2;
            
            ++forcecount;
        }// if force
    }// if node exists
    else
    {
        delete cons;
        delete tmpnode;
        printf("ERROR: node R[%d, %d, %d], V[%f, %f, %f], C[%d, %d, %d] does not exists, cannot add constraint \n", x, y, z,vx,vy,vz,cx,cy,cz);
        OK = false;
        return OK;
    }
    return OK;
}// AddConstraint()


// Delete node
bool pFESolver::RemoveNode(const xyzType nx, const xyzType ny, const xyzType nz)
{
    
    bool OK(true);
    return OK;
}

// Delete an element to the model (along with its nodes)
bool pFESolver::RemoveElement(const midxType material, const xyzType x, const xyzType y, const xyzType z)
{
    
    bool OK(true);
    return OK;
}

// Delete constraint
bool pFESolver::RemoveConstraint(const xyzType x, const xyzType y, const xyzType z)
{
    
    bool OK(true);
    return OK;
}

//========== Getters

// Get a particular node by coords
Node<xyzType>* pFESolver::GetNode(const xyzType x, const xyzType y, const xyzType z)
{
    Node<xyzType>* tmpnode = new Node<xyzType>();
    tmpnode->x = x;
    tmpnode->y = y;
    tmpnode->z = z;
    NodeSet_it nodeptr =  NodeS.find(tmpnode);
    
    return (*nodeptr);
}// GetNode()

// Gets a particular element by coords
Elem<xyzType>* pFESolver::GetElement(const xyzType x, const xyzType y, const xyzType z)
{
    return (*ElemS.begin());
}// GetElement()

// Gets a set of all the (local) elements
pFESolver::ElementSet pFESolver::GetLocalElements()
{
    return ElemS;
}// GetLocalElements()

// Gets a set of all the (local) nodes
pFESolver::NodeSet pFESolver::GetLocalNodes()
{
    return NodeS;
}// GetLocalNodes()


Cons<xyzType>* pFESolver::GetConstraint(const xyzType x, const xyzType y, const xyzType z)
{
    Cons<xyzType>* constraint = (GetNode(x,y,z))->constraint;
    return constraint;
}


//================== Allocate system matrix / vector rows

// allocate gsm rows and columns according to number of processes
PetscErrorCode pFESolver::AllocateLocalMatrix(Mat *gsm)
{
    if(MPIrank == MASTER)
        printf("%d: NodeS.size() = %d \n",MPIrank, NodeS.size());
    
    if(MPIsize == 1)
    {// IF SERIAL
        localgsmrows = globalgsmrows;
        localgsmcols = globalgsmcols;
        localrhslength = globalgsmcols;
        localsollength = globalgsmcols;
        
        ierr = MatCreateSeqAIJ(comm, globalgsmrows, globalgsmcols, NUM_TERMS, PETSC_NULL, gsm); CHKERRQ(ierr);
        ierr = MatSeqAIJSetPreallocation(*gsm, NUM_TERMS, PETSC_NULL); CHKERRQ(ierr);
    }
    else
    {// IF PARALLEL
        double gr(globalgsmrows);
        double gc(globalgsmcols);
        localgsmrows = ceil(gr/MPIsize);
        localgsmcols = ceil(gc/MPIsize);
     
        diagonalterms = NUM_TERMS;
        offdiagonalterms = NUM_TERMS;
        
        localrhslength = localgsmcols;
        localsollength = localgsmcols;
 
        ierr = MatCreate(comm, gsm); CHKERRQ(ierr);
        ierr = MatSetSizes(*gsm, localgsmrows, localgsmcols, globalgsmrows, globalgsmcols); CHKERRQ(ierr);
        ierr = MatSetType(*gsm, MATMPIAIJ); CHKERRQ(ierr);
        ierr = MatMPIAIJSetPreallocation(*gsm, diagonalterms, PETSC_NULL, offdiagonalterms, PETSC_NULL); CHKERRQ(ierr);
    }
    
    // check rows
    PetscInt grows, gcols;
    ierr = MatGetSize(*gsm, &grows, &gcols); CHKERRQ(ierr);
    
    if(MPIrank == MASTER)
        printf("%d: MatGetSize of GSM : rows=%d, cols=%d\n",MPIrank,grows, gcols);
    
    return 0;
}

// Allocate RHS rows according to number of processes
PetscErrorCode pFESolver::AllocateLocalVec(Vec *vec)
{
    if(MPIsize == 1)
    {// IF SERIAL
        ierr = VecSetSizes(*vec, globalgsmrows, globalgsmrows); CHKERRQ(ierr);
        ierr = VecSetType(*vec, VECSEQ); CHKERRQ(ierr);
    }
    else
    {// IF PARALLEL
        ierr = VecSetSizes(*vec, localgsmcols, globalgsmrows); CHKERRQ(ierr);
        ierr = VecSetType(*vec, VECMPI); CHKERRQ(ierr);
    }
    return 0;
}

//=============== System matrices

// Build Global Stiffness Matrix
PetscErrorCode pFESolver::ComputeGSM(Mat *GSM)
{
    if(MPIrank==MASTER)
        printf("%d: In GSM\n", MPIrank);
    
    globalgsmrows =  NodeS.size()*DOF_3D;
    globalgsmcols =  NodeS.size()*DOF_3D; // NUM_TERMS = max number of non-zero elements per row
    MatInfo matinfo; // use to get matrix info
    
    // allocate GSM rows
    pFESolver::AllocateLocalMatrix(GSM);
    
    ////////////
    // Currently, all PETSc parallel matrix formats are partitioned by
    // contiguous chunks of rows across the processors.  Determine which
    // rows of the matrix are locally owned.
    ///////////
    
    // inclusive of first, exclusive of last
    PetscInt localfirst, locallast;
    ierr = MatGetOwnershipRange(*GSM, &localfirst, &locallast); CHKERRQ(ierr);
    
    printf("%d: GSM local : first row=%d, last row=%d \n", MPIrank,localfirst, locallast-1);
    
    idxType gsmcolcount(0);
    std::map<idxType, idxType> tmp_gsmcolidx;
    idxType currentcol(0);
    PetscInt numcols(0);
    PetscInt gsmCol[NUM_TERMS];
    PetscInt gsmRowX[1], gsmRowY[1], gsmRowZ[1];
    
    PetscScalar gsmvalRowX[NUM_TERMS]; // row-centric storage etc.
    PetscScalar gsmvalRowY[NUM_TERMS];
    PetscScalar gsmvalRowZ[NUM_TERMS];
    
    bool consnodeX(1); // is curnode constrained in X dim? Default is 1 (true)
    bool consnodeY(1); // is curnode constrained in Y dim? Default is 1 (true)
    bool consnodeZ(1); // is curnode constrained in Z dim? Default is 1 (true)
    
    bool consnodeNeighbourX(1); // is curennode constrained in X dim? Default is 1  (true)
    bool consnodeNeighbourY(1); // is curennode constrained in X dim? Default is 1  (true)
    bool consnodeNeighbourZ(1); // is curennode constrained in X dim? Default is 1  (true)
    
    bool constraintX[DOF_3D] = {1,1,1}; // should X dim term be added to GSM? Default is 1 (yes)
    bool constraintY[DOF_3D] = {1,1,1}; // should Y dim term be added to GSM? Default is 1 (yes)
    bool constraintZ[DOF_3D] = {1,1,1}; // should Z dim term be added to GSM? Default is 1 (yes)
    
    ////////////////////////
    // BUILD GSM BY LOOPING THROUGH NODES IN NodeS
    ///////////////////////
    
    unsigned int lsmlen = NODES_PER_ELEMENT * DOF_3D; // 24
    
    idxType nodecount(0); // renumbering of nodes
    idxType nodes[NODES_PER_ELEMENT], renumNode(0), renumNodeNeighbour(0), nodeL(0), nodeNeighbourL(0);
    idxType nodeG(0), elemNeighbourG(0), nodeNeighbourG(0);
    idxType e, enn;
    midxType materialidx;
    xyzType c, colindex;
    
    idxType count(0);
    
    NodeSet_const_it cnitr;
    
    for(cnitr = NodeS.begin(); cnitr != NodeS.end(); ++cnitr)
    {

        renumNode = pFESolver::GetNodeIndex(cnitr); // renumbered or local index
        
        if(renumNode >= localfirst && renumNode < locallast){ // if this localrow is on rank
            
            // reset stuff...
            gsmcolcount = 0;
            tmp_gsmcolidx.clear();
            
            for(c=0; c<NUM_TERMS; ++c){
                gsmCol[c] = 0;
                gsmvalRowX[c]=0;
                gsmvalRowY[c]=0;
                gsmvalRowZ[c]=0;
            }
            
            gsmRowX[0] = renumNode*DOF_3D  + 0;
            gsmRowY[0] = renumNode*DOF_3D  + 1;
            gsmRowZ[0] = renumNode*DOF_3D  + 2;
            
            Cons<xyzType>* nodecons = pFESolver::GetNodeCons(cnitr);
            if(nodecons){
                consnodeX = nodecons->cx;
                consnodeY = nodecons->cy;
                consnodeZ = nodecons->cz;
            } else{
                consnodeX = 1;
                consnodeY = 1;
                consnodeZ = 1;
            }
            
            // enn = neighbour, e = elem
            for(int elem = 0; elem < NODES_PER_ELEMENT; ++elem){// for each element the node belongs to (max 8)
                
                if(pFESolver::GetNodeElement(cnitr, elem)){
                    materialidx = pFESolver::GetElementMaterial(cnitr, elem);
                    nodeL = elem;
                    
                    for(int neighbour=0; neighbour<NODES_PER_ELEMENT; ++neighbour){// for each neighbouring node on element e
                        renumNodeNeighbour = pFESolver::GetNodeNeighbourIndex(cnitr, elem, neighbour); //get renumbered index
                        
                        if(tmp_gsmcolidx.find(renumNodeNeighbour) == tmp_gsmcolidx.end()){ // if neighbouring doesn't already have a column number
                            tmp_gsmcolidx[renumNodeNeighbour] = gsmcolcount;
                            ++gsmcolcount;
                        }
                        currentcol = tmp_gsmcolidx[renumNodeNeighbour];
                        
                        Cons<xyzType>* neighbourcons = pFESolver::GetNodeNeighbourCons(cnitr, elem, neighbour);
                        if(neighbourcons){
                            consnodeNeighbourX = neighbourcons->cx; // constrained //
                            consnodeNeighbourY = neighbourcons->cy; // constrained //
                            consnodeNeighbourZ = neighbourcons->cz; // constrained //
                        }else{
                            consnodeNeighbourX = 1;
                            consnodeNeighbourY = 1;
                            consnodeNeighbourZ = 1;
                        }
                        
                        nodeNeighbourL = neighbour; // local index of neighbouring node on element e
                        
                        constraintX[0] = consnodeX && consnodeNeighbourX;
                        constraintX[1] = consnodeX && consnodeNeighbourY;
                        constraintX[2] = consnodeX && consnodeNeighbourZ;
                        
                        constraintY[0] = consnodeY && consnodeNeighbourX;
                        constraintY[1] = consnodeY && consnodeNeighbourY;
                        constraintY[2] = consnodeY && consnodeNeighbourZ;
                        
                        constraintZ[0] = consnodeZ && consnodeNeighbourX;
                        constraintZ[1] = consnodeZ && consnodeNeighbourY;
                        constraintZ[2] = consnodeZ && consnodeNeighbourZ;
                        
                        if(nodeL == nodeNeighbourL){
                            if(!consnodeX){
                                constraintX[0] = 1; constraintX[1] = 0; constraintX[2] = 0;
                            }// if consnodeX
                            if(!consnodeY){
                                constraintY[0] = 0; constraintY[1] = 1; constraintY[2] = 0;
                            }// if consnodeY
                            if(!consnodeZ){
                                constraintZ[0] = 0; constraintZ[1] = 0; constraintZ[2] = 1;
                            }// if consnodeZ
                        }// if diagonal element of GSM
                        
                        for(c=0; c<3; ++c){// for each x,y,z component/column
                            colindex = currentcol*DOF_3D + c;
                            gsmCol[colindex] = renumNodeNeighbour*DOF_3D + c;
                            
                            int xrowidx = (nodeL*DOF_3D + 0)*lsmlen + nodeNeighbourL*DOF_3D + c;
                            int yrowidx = (nodeL*DOF_3D + 1)*lsmlen + nodeNeighbourL*DOF_3D + c;
                            int zrowidx = (nodeL*DOF_3D + 2)*lsmlen + nodeNeighbourL*DOF_3D + c;
                            
                            
                            // add gsm value. Access correct LSM in MATERIALMAP using material index materialidx
                            gsmvalRowX[colindex] += pFESolver::GetLSMValue( materialidx, xrowidx) * constraintX[c]; // x row
                            gsmvalRowY[colindex] += pFESolver::GetLSMValue( materialidx, yrowidx ) * constraintY[c]; // y row
                            gsmvalRowZ[colindex] += pFESolver::GetLSMValue( materialidx, zrowidx ) * constraintZ[c]; // z row
                            
                        }// for each component
                        
                    }// for each neighbouring node enn on element e
                }// if valid element e
            }// for each element e of current node
            
            numcols = DOF_3D * gsmcolcount;
            
            ierr = MatSetValues(*GSM, 1, gsmRowX, numcols, gsmCol, gsmvalRowX, INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValues(*GSM, 1, gsmRowY, numcols, gsmCol, gsmvalRowY, INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValues(*GSM, 1, gsmRowZ, numcols, gsmCol, gsmvalRowZ, INSERT_VALUES); CHKERRQ(ierr);
            
            ++count;
            
        }// if on local rank
        
    }// for each node
    
    ierr = MatAssemblyBegin(*GSM, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*GSM, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    
    // Get info about GSM allocation etc.
    MatInfo info;
    double  mal, nz_a, nz_u, nz_un;
    
    ierr = MatGetInfo(*GSM, MAT_LOCAL, &info); CHKERRQ(ierr);
    mal  = info.mallocs;
    nz_a = info.nz_allocated;
    nz_u = info.nz_used;
    nz_un = info.nz_unneeded;
    
    printf("%d: GSM local info : mal = %lf, non-zero_allocated = %lf, non-zero_used = %lf, non-zero_unneeded = %lf \n", MPIrank, mal, nz_a, nz_u, nz_un);
    if(MPIrank == MASTER)
        printf("%d: Leaving GSM\n",MPIrank);
    
    return 0;
    
}// computeGSM()


// Build RHS force vector
PetscErrorCode pFESolver::ComputeRHS(Vec *rhs){
    
    PetscInt xs,ys,zs,nx,ny,nz;
    int ii,jj,kk;
    PetscInt size;
    
    if(MPIrank==MASTER)
        printf("%d: In computeRHS\n", MPIrank);
    
    ierr = VecCreate(comm, rhs); CHKERRQ(ierr);
    
    ierr = AllocateLocalVec(rhs);
    
    ierr = VecGetSize(*rhs, &size); CHKERRQ(ierr);
    
    if(MPIrank==MASTER)
        printf("%d: VecGetSize : size=%d\n",MPIrank,size);
    
    ierr = VecSet(*rhs,0); CHKERRQ(ierr);
    ierr = VecSetValues(*rhs, totalrhs, forcecolidx, forcevalue, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(*rhs); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(*rhs); CHKERRQ(ierr);
    
    
    if(MPIrank==MASTER)
        printf("%d: Leaving RHS \n", MPIrank);
    
    return 0;
    
}// ComputeRHS()


//========== Solve

// Solves the system, updating the node displacements
// COLLECTIVE

// only call this function if setup done and/or rebuild done

// Solve linear elastic problem
PetscErrorCode pFESolver::Solve(){
    
    idxType ii(0);
    xyzType xx(0), yy(0), zz(0);
    idxType  nodeidx;
    FILE * outFile;
    
    PC prec;              // preconditioner
    Mat GSM;              // GSM
    Vec sol;              // solution
    Vec rhs;              // rhs (force) vector
    PetscInt iters;       // number of iterations
    PetscReal norm;       // norm

    if( MPIrank == MASTER )
        printf("%d: In solve \n", MPIrank);
    
    ierr = KSPCreate(comm, &ksp); CHKERRQ(ierr);
    
    if( MPIrank == MASTER )
        printf("%d: Created KSP \n", MPIrank);
    
    ierr = pFESolver::ComputeGSM(&GSM); CHKERRQ(ierr);
    ierr = pFESolver::ComputeRHS(&rhs); CHKERRQ(ierr);
    
    ierr = VecCreate(comm, &sol); CHKERRQ(ierr);
    
    pFESolver::AllocateLocalVec(&sol);
    
#if PETSC_VERSION_LT(3,5,1)
    KSPSetOperators(ksp, GSM, GSM, DIFFERENT_NONZERO_PATTERN);
#else
    KSPSetOperators(ksp, GSM, GSM);
#endif
    
    // KSPCG
    ierr = KSPSetType(ksp,KSPCG); CHKERRQ(ierr);
    ierr = KSPGetPC(ksp,&prec); CHKERRQ(ierr);
    ierr = PCSetType(prec,PCJACOBI); CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp, TOLERANCE, PETSC_DEFAULT, PETSC_DEFAULT, MAXITER); CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(ksp, PETSC_FALSE); CHKERRQ(ierr);
    
    ierr = KSPSetUp(ksp); CHKERRQ(ierr);
    
    if(MPIrank == MASTER)
        printf("%d: Starting solver...\n",MPIrank);
    
    double startT = MPI_Wtime();
    
    ierr = KSPSolve(ksp, rhs, sol); CHKERRQ(ierr);

    double endT = MPI_Wtime();
    
    if(MPIrank == MASTER)
    {
        printf("%d: Done solving!\n", MPIrank);
        printf("%d: SOLVETIME = %.5le \n", MPIrank, endT - startT);
    }
    
    ierr = KSPGetSolution(ksp, &sol); CHKERRQ(ierr);
    
    ierr = KSPGetIterationNumber(ksp, &iters); CHKERRQ(ierr);
    ierr = KSPGetResidualNorm(ksp, &norm); CHKERRQ(ierr);
    
    
    if(MPIrank == MASTER)
        printf("%d: Converged to %f in %d iterations.\n", MPIrank, norm, iters);
    
    if(MPIrank == MASTER)
        printf("%d: About to print results\n",MPIrank);
    
    PetscScalar tmpx, tmpy, tmpz;
    PetscInt nlocal;
    ierr = VecGetLocalSize(sol, &nlocal); CHKERRQ(ierr);
    
    PetscInt vstart, vend;
    ierr = VecGetOwnershipRange(sol, &vstart, &vend); CHKERRQ(ierr);
    
    VecScatter vsctx;
    ierr = VecScatterCreateToZero(sol, &vsctx, &vecout); CHKERRQ(ierr);
    ierr = VecScatterBegin(vsctx, sol, vecout, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(vsctx, sol, vecout, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    
    printf("%d: vstart=%d, vend=%d \n", MPIrank, vstart, vend);
    
    if( MPIrank == MASTER ){
        ierr = VecGetArray(vecout, &solution); CHKERRQ(ierr);
        
        PetscInt vsize;
        ierr = VecGetSize(vecout,&vsize); CHKERRQ(ierr);
        printf("%d: vecout size = %d \n",MPIrank,vsize);
        
        // update nodes with final displacements
        NodeSet_it nitr = NodeS.begin();
        Node<xyzType>* nodeptr;
        
        for (nitr = NodeS.begin(); nitr != NodeS.end(); ++nitr){
            nodeptr = (*nitr);
            
            xx = nodeptr->x;
            yy = nodeptr->y;
            zz = nodeptr->z;
            ii = nodeptr->idx;
            // update node displacements
            nodeptr->dx = solution[ii*DOF_3D + 0];
            nodeptr->dy = solution[ii*DOF_3D + 1];
            nodeptr->dz = solution[ii*DOF_3D + 2];
            
        }// for each node in Node Set
        
        ierr = VecRestoreArray(vecout, &solution); CHKERRQ(ierr);
        
    }// if MASTER process
    
    
    ierr = VecScatterDestroy(&vsctx); CHKERRQ(ierr);
    
    delete[] forcevalue;
    delete[] forcecolidx;
    
    if(MPIrank==MASTER)
        printf("%d: Leaving Solve\n",MPIrank);
    
    
    SOLVE_DONE = true;
    
    return 0;
    
}// Solve()


//=========== Finish up

PetscErrorCode pFESolver::cleanup()
{
    // add other final operations here
    ierr = VecDestroy(&vecout); CHKERRQ(ierr);
    
    // FIX ME : check if set
    KSPDestroy(&ksp);
    if(MPIrank==MASTER)
        printf("%d: All done, bye, bye\n",MPIrank);
    
    return 0;
}// cleanup()


// Renumber nodes in NodeSet
bool pFESolver::RenumberNodes()
{
    bool OK(true);
    
    try
    {
        idxType ncount(0);
        idxType tmpidx(0);
        NodeSet_it nitr;
        
        for(nitr = NodeS.begin(); nitr != NodeS.end(); ++nitr)
        {
            tmpidx = (*nitr)->idx;
            (*nitr)->idx = ncount;
            ++ncount;
        }// for each node
    }
    catch(std::exception &e)
    {
        printf("ERROR: could not renumber node\n");
        OK = false;
    }
    return OK;
}

