//
//  IntegrationMtx.h
//  
//
//  Created by N. Banglawala
//
//

#ifndef GUARD_INTEGRATIONMTX_H
#define GUARD_INTEGRATIONMTX_H


#include "Common.h"


class IntegrationMtx
{

public:
    
    IntegrationMtx();
    ~IntegrationMtx();
    
    void create();
    
    // inline
    double operator()(const unsigned int idxcol, const unsigned int idxrow ) const
    { return  matrix[idxcol][idxrow]; }

private:

    double  matrix[8][DOF_3D];
};

#endif /* defined( GUARD_INTEGRATIONMTX_H ) */
