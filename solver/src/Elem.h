//
// Element class
//
//  Elem.h
//  
//
//  Created by N. Banglawala
//
//

#ifndef GUARD_ELEM_H
#define GUARD_ELEM_H

#include "Common.h"

template <typename T> class Node;

template <typename T>
class Elem{
public :
    T         x, y, z;
    idxType   idx;
    midxType  material; // material idx (only <= 2^8 = 256 expected)
    Node<T> * nodes[NODES_PER_ELEMENT];
    Elem();
    Elem(midxType M);
    ~Elem();
};

#endif /* defined(GUARD_ELEM_H) */

