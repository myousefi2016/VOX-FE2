//
//  GradientMtx.h
//  
//
//  Created by N. Banglawala
//
//

#ifndef GUARD_GRADIENTMTX_H
#define GUARD_GRADIENTMTX_H

#include "Common.h"
#include "IntegrationMtx.h"


class GradientMtx
{
    
public:
    GradientMtx();
    GradientMtx(double a_size, double b_size, double c_size, int n_num);
    //GradientMtx(double a_size, double b_size, double c_size, int n_num, IntegrationMtx integMtx);
    ~GradientMtx();
    
    // inline
    // get underlying matrix and transpose
    Mat* getmat(){ return (&matrix); }
    Mat* getmat_t(){ return (&transpose); }

protected:
    
    int n; // nth node or integration point ?
    double a, b, c; // voxel sizes (x, y, z)
    
    IntegrationMtx integMtx;
    
    Mat matrix; // petsc matrix
    Mat transpose; // transpose of pmatrix
    PetscErrorCode ierr; // petsc errorcode
    
    
    PetscErrorCode create();
    PetscErrorCode destroy();

    
};

#endif /* defined(GUARD_GRADIENTMTX_H) */
