//
//  GradientMtx.cpp
//
//
//  Created by N. Banglawala
//
//

#include "GradientMtx.h"

GradientMtx::GradientMtx()
{
    create();
}

GradientMtx::GradientMtx(const double a_size, const double b_size, const double c_size, const int n_num) : a(a_size), b(b_size), c(c_size), n(n_num)
{
    create();
}

GradientMtx::~GradientMtx()
{
    //destroy();
}

PetscErrorCode GradientMtx::destroy()
{
    ierr = MatDestroy(&matrix);  CHKERRQ(ierr);
    ierr = MatDestroy(&transpose); CHKERRQ(ierr);
    return 0;
}

PetscErrorCode GradientMtx::create()
{
    double temp;
    
    // now the gradient (or B) matrices
    ierr = MatCreateSeqDense(PETSC_COMM_SELF, 6, 24, PETSC_NULL, &matrix);  CHKERRQ(ierr);
    
    for (int i=0; i<8; i++)
    {// 8 integration points
        temp = 1.0 / a * SIGN(integMtx(i,0)) *
        ( 0.5 + integMtx(n,1) * SIGN(integMtx(i,1)) ) *
        ( 0.5 + integMtx(n,2) * SIGN(integMtx(i,2)) );
        
        ierr = MatSetValue(matrix, 0, 0+3*i, temp, INSERT_VALUES); CHKERRQ(ierr);
        ierr = MatSetValue(matrix, 3, 1+3*i, temp, INSERT_VALUES); CHKERRQ(ierr);
        ierr = MatSetValue(matrix, 5, 2+3*i, temp, INSERT_VALUES); CHKERRQ(ierr);
        
        temp = 1.0 / b * SIGN(integMtx(i,1)) *
        ( 0.5 + integMtx(n,0) * SIGN(integMtx(i,0)) ) *
        ( 0.5 + integMtx(n,2) * SIGN(integMtx(i,2)) );
        
        
        ierr = MatSetValue(matrix, 1, 1+3*i, temp, INSERT_VALUES); CHKERRQ(ierr);
        ierr = MatSetValue(matrix, 3, 0+3*i, temp, INSERT_VALUES); CHKERRQ(ierr);
        ierr = MatSetValue(matrix, 4, 2+3*i, temp, INSERT_VALUES); CHKERRQ(ierr);
        
        temp = 1.0 / c * SIGN(integMtx(i,2)) *
        ( 0.5 + integMtx(n,0) * SIGN(integMtx(i,0)) ) *
        ( 0.5 + integMtx(n,1) * SIGN(integMtx(i,1)) );
        
        
        ierr = MatSetValue(matrix, 2, 2+3*i, temp, INSERT_VALUES); CHKERRQ(ierr);
        ierr = MatSetValue(matrix, 4, 1+3*i, temp, INSERT_VALUES); CHKERRQ(ierr);
        ierr = MatSetValue(matrix, 5, 0+3*i, temp, INSERT_VALUES); CHKERRQ(ierr);
    }
    
    ierr = MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    
    // now create transpose
    MatTranspose(matrix, MAT_INITIAL_MATRIX, &transpose);
       
    return 0;
}


