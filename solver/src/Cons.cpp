//
// Constraint class
//
//  Cons.cpp
//  
//
//  Created by N. Banglawala 
//
//

#include "Cons.h"


template <typename T> Cons<T>::Cons() : cx(1), cy(1), cz(1), node(0), preserve(1) { }

template <typename T> Cons<T>::Cons(bool CX, bool CY, bool CZ) : cx(CX), cy(CY), cz(CZ), node(0), preserve(1) { }

template <typename T> Cons<T>::~Cons() { }

//=================================================

template class Cons<xyzType>;

