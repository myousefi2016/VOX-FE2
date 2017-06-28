//
//  PropertyMtx.h
//  
//
//  Created by N. Banglawala
//
//

#ifndef GUARD_PROPERTYMTX_H
#define GUARD_PROPERTYMTX_H

#include "Common.h"

class PropertyMtx
{
    
public:
    
    PropertyMtx();
    PropertyMtx(const double ym, const double pr);
    ~PropertyMtx();

    
    // inline
    // get underlying matrix
    Mat* getmat() { return (&matrix); }
  

    
protected:
    double youngsm; // young's modulus
    double poissonsr; // poisson's ratio
    
    Mat matrix; // petsc matrix
    PetscErrorCode ierr; // petsc error code

    PetscErrorCode create();
    PetscErrorCode destroy();
    
    
};



#endif /* defined(GUARD_PROPERTYMTX_H) */
