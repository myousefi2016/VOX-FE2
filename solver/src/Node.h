//
// Node class
//
//
//  Node.h
//  
//
//  Created by N. Banglawala
//
//

#ifndef GUARD_NODE_H
#define GUARD_NODE_H

#include "Common.h"
#include "petscsys.h"

template <typename T> class Elem;
template <typename T> class Cons;

template<typename T> class Node
{
public:
    T           x,y,z;      // coordinates
    PetscScalar dx, dy, dz; // displacements
    idxType     idx;        // index
    Cons<T>*    constraint; // constraint pointer
    Elem<T>*    elems[NODES_PER_ELEMENT];
    Node();
    ~Node();


};


#endif /* defined(GUARD_NODE_H) */

