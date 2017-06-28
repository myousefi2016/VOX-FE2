//
//  LocalStiffnessMtx.h
//  
//
//  Created by N. Banglawala
//
//

#ifndef GUARD_LOCALSTIFFNESSMTX_H
#define GUARD_LOCALSTIFFNESSMTX_H

#include "Common.h"
#include "GradientMtx.h"
#include "PropertyMtx.h"


class LocalStiffnessMtx
{

public:
    double youngsm; // youngs modulus
    double poissonsr; // poissons ratio
    double a, b, c; // voxel sizes (x, y, z)
    
    PetscErrorCode ierr; // petsc errorcode
    
    PetscScalar *matrix; // lsm matrix --> actually, it's a petscscalar! (access to)
    
    LocalStiffnessMtx();
    LocalStiffnessMtx(const double a, const double b, const double c, const double ym, const double pr);
    LocalStiffnessMtx(const double a, const double b, const double c, const double ym, const double pr, GradientMtx * GradientMatrices);
    ~LocalStiffnessMtx();
    
    // inline
    // get matrix
    PetscScalar* getmat(){ return (matrix); }
    
    double getval(const unsigned int index)
    {
        return matrix[index];
    }
    
    
    
protected:
    
    PetscErrorCode create(GradientMtx *GradientMatrices);
    PetscErrorCode create();
    void destroy();

};

#endif /* defined( GUARD_LOCALSTIFFNESSMTX_H ) */
