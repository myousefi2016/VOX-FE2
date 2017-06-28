//
//  PropertyMtx.cpp
//  
//
//  Created by N. Banglawala
//
//

#include "PropertyMtx.h"

// use default values
PropertyMtx::PropertyMtx() : youngsm(1), poissonsr(0.3)
{
    create();
}

PropertyMtx::PropertyMtx(const double ym, const double pr) : youngsm(ym), poissonsr(pr)
{
    create();
}

PropertyMtx::~PropertyMtx()
{
    destroy();
}

PetscErrorCode PropertyMtx::destroy()
{
    ierr = MatDestroy(&matrix); CHKERRQ(ierr);
    return 0;
}

// Build Property Matrix
PetscErrorCode PropertyMtx::create()
{
    ierr = MatCreateSeqDense(PETSC_COMM_SELF, 6, 6, PETSC_NULL, &matrix); CHKERRQ(ierr);
    
    ierr = MatSetValue(matrix, 0, 0, 1.0, INSERT_VALUES); CHKERRQ(ierr);
    ierr = MatSetValue(matrix, 0, 1, poissonsr/(1.0-poissonsr), INSERT_VALUES); CHKERRQ(ierr);
    ierr = MatSetValue(matrix, 0, 2, poissonsr/(1.0-poissonsr), INSERT_VALUES); CHKERRQ(ierr);
    ierr = MatSetValue(matrix, 1, 0, poissonsr/(1.0-poissonsr), INSERT_VALUES); CHKERRQ(ierr);
    ierr = MatSetValue(matrix, 1, 1, 1.0, INSERT_VALUES); CHKERRQ(ierr);
    ierr = MatSetValue(matrix, 1, 2, poissonsr/(1.0-poissonsr), INSERT_VALUES); CHKERRQ(ierr);
    ierr = MatSetValue(matrix, 2, 0, poissonsr/(1.0-poissonsr), INSERT_VALUES); CHKERRQ(ierr);
    ierr = MatSetValue(matrix, 2, 1, poissonsr/(1.0-poissonsr), INSERT_VALUES); CHKERRQ(ierr);
    ierr = MatSetValue(matrix, 2, 2, 1.0, INSERT_VALUES); CHKERRQ(ierr);
    ierr = MatSetValue(matrix, 3, 3, (1.0-2.0*poissonsr)/(2.0*(1.0-poissonsr)), INSERT_VALUES); CHKERRQ(ierr);
    ierr = MatSetValue(matrix, 4, 4, (1.0-2.0*poissonsr)/(2.0*(1.0-poissonsr)), INSERT_VALUES); CHKERRQ(ierr);
    ierr = MatSetValue(matrix, 5, 5, (1.0-2.0*poissonsr)/(2.0*(1.0-poissonsr)), INSERT_VALUES); CHKERRQ(ierr);
    
    ierr = MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    
    ierr = MatScale(matrix, youngsm*(1.0-poissonsr)/((1.0+poissonsr)*(1.0-2.0*poissonsr)));  CHKERRQ(ierr);
    
    return 0;
}



