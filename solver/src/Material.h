//
// Material class
//
//
//  Material.h
//  
//
//  Created by N. Banglawala 
//
//

#ifndef GUARD_MATERIAL_H
#define GUARD_MATERIAL_H

#include "Common.h"
#include "LocalStiffnessMtx.h"
#include "petscsys.h"

class Material
{
    
public:
    midxType        idx;
    double          youngsm;
    double          poissonsr;
    bool            remodel;
    LocalStiffnessMtx lsm;
    
    Material();
    Material(midxType id, double ym, double pr, LocalStiffnessMtx mat, bool rf);
    ~Material();
};


#endif /* defined(GUARD_MATERIAL_H) */

