//
// Node class
//
//
//  Node.cpp
//  
//
//  Created by N. Banglawala
//
//

#include "Node.h"

template <typename T> Node<T>::Node() : constraint(0), dx(0), dy(0), dz(0)
{
    for(int kk=0; kk<NODES_PER_ELEMENT; ++kk)   elems[kk] = 0;
}

template <typename T> Node<T>::~Node() { }


//=============================================================

template class Node<xyzType>;

