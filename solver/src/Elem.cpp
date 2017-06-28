//
// Element class
//
// Elem.cpp
//  
//
//  Created by N. Banglawala
//
//

#include "Elem.h"

//====================================================

template <typename T> Elem<T>::Elem()
{
    for(int kk=0; kk<NODES_PER_ELEMENT; ++kk) nodes[kk] = 0;
}

template <typename T> Elem<T>::Elem(midxType M) : material(M)
{
    for(int kk=0; kk<NODES_PER_ELEMENT; ++kk) nodes[kk] = 0;
}

template <typename T> Elem<T>::~Elem() { }

//====================================================

template class Elem<xyzType>;

