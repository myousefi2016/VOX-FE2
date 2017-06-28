//
//  LocalStiffnessMtx.cpp
//  
//
//  Created by N. Banglawala
//
//

#include "LocalStiffnessMtx.h"


LocalStiffnessMtx::LocalStiffnessMtx() : a(1.0), b(1.0), c(1.0), youngsm(1.0), poissonsr(0.3)
{}

LocalStiffnessMtx::LocalStiffnessMtx(const double a_size, const double b_size, const double c_size, const double ym, const double pr) : a(a_size), b(b_size), c(c_size), youngsm(ym), poissonsr(pr)
{
    create();
}

LocalStiffnessMtx::LocalStiffnessMtx(const double a_size, const double b_size, const double c_size, const double ym, const double pr, GradientMtx *GradientMatrices) :
a(a_size), b(b_size), c(c_size), youngsm(ym), poissonsr(pr)
{
    create(GradientMatrices);
}

LocalStiffnessMtx::~LocalStiffnessMtx()
{}


PetscErrorCode LocalStiffnessMtx::create()
{
    GradientMtx gradientMtx[NODES_PER_ELEMENT];
    
    // create gradient matrices
    for(int n=0; n<NODES_PER_ELEMENT; ++n)
    {
        gradientMtx[n] = GradientMtx(a, b, c, n);
    }
    
    // create LSM
    ierr = create(gradientMtx); CHKERRQ(ierr);
    
    return 0;
}


PetscErrorCode LocalStiffnessMtx::create(GradientMtx * gradientMtx)
    {
        matrix = new PetscScalar();
        
        Mat LSM, propXgradMtx, tempMtx;
        int lsmlen = NODES_PER_ELEMENT * DOF_3D; // 24
        
        // Build Property Matrix
        PropertyMtx propMtx(youngsm, poissonsr);
        
        ierr = MatCreateSeqDense(PETSC_COMM_SELF, lsmlen, lsmlen, PETSC_NULL, &LSM); CHKERRQ(ierr);
        
        // assuming all 8 gradient matrices are available
        for (int n=0; n < NODES_PER_ELEMENT; ++n)
        {// per node n
            
            ierr = MatMatMult(*(propMtx.getmat()), *(gradientMtx[n].getmat()), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &propXgradMtx); CHKERRQ(ierr);
            ierr = MatMatMult(*(gradientMtx[n].getmat_t()), propXgradMtx, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &tempMtx); CHKERRQ(ierr);
            
            ierr = MatAXPY(LSM, 1.0, tempMtx, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
            
            ierr = MatDestroy(&propXgradMtx); CHKERRQ(ierr);
            ierr = MatDestroy(&tempMtx); CHKERRQ(ierr);
            
        }// per integration point i, gradient matrix
        
        ierr = MatScale(LSM, 0.125 * a * b * c); CHKERRQ(ierr);
        
        // Get LSM
        ierr = MatDenseGetArray(LSM, (&matrix)); CHKERRQ(ierr);
        
        return 0;
}


