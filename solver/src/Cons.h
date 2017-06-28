//
// Constraint class
//
//  Cons.h
//  
//
//  Created by N. Banglawala
//
//

#ifndef GUARD_CONS_H
#define GUARD_CONS_H

#include "Common.h"

template <typename T> class Node;

template <typename T> class Cons
{
    
public:
    bool      cx, cy, cz;  // 0 if constrained, 1 if not
    double    vx, vy, vz;  // value of constraint
    bool      preserve; // preserve the node for remodelling
    // double  dx, dy, dz;  // RHS adjustment due to constraint (for fixed displacements)
    Node<T> * node;        // points to constrained node
    Cons();
    Cons(bool CX, bool CY, bool CZ);
    ~Cons();
};


#endif /* defined(GUARD_CONS_H) */

