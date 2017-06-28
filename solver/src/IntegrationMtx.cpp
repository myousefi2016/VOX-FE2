//
//  IntegrationMtx.cpp
//  
//
//  Created by N. Banglawala
//
//

#include "IntegrationMtx.h"

// not a petsc matrix !

IntegrationMtx::IntegrationMtx() // a(1), b(1) c(1)
{
    create();
}

IntegrationMtx::~IntegrationMtx()
{}

void IntegrationMtx::create()
{
    
    // Create the Integration Matrix
    matrix[0][0] = -1; matrix[0][1] = -1; matrix[0][2] = -1;
    matrix[1][0] = +1; matrix[1][1] = -1; matrix[1][2] = -1;
    matrix[2][0] = -1; matrix[2][1] = +1; matrix[2][2] = -1;
    matrix[3][0] = +1; matrix[3][1] = +1; matrix[3][2] = -1;
    matrix[4][0] = -1; matrix[4][1] = -1; matrix[4][2] = +1;
    matrix[5][0] = +1; matrix[5][1] = -1; matrix[5][2] = +1;
    matrix[6][0] = -1; matrix[6][1] = +1; matrix[6][2] = +1;
    matrix[7][0] = +1; matrix[7][1] = +1; matrix[7][2] = +1;
    
    for (int j=0; j<DOF_3D; j++)
    {
        for (int i=0; i<8; i++)
        {// 8 integration points
            matrix[i][j] = matrix[i][j]/(2*sqrt(3)); // R in old solver...
        }
    }
}
