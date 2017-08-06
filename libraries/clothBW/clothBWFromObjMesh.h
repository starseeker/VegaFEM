/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 3.0                               *
 *                                                                       *
 * "clothBW" library , Copyright (C) 2016 USC                            *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Andy Pierce, Yijing Li, Yu Yu Xu, Jernej Barbic         *
 * http://www.jernejbarbic.com/code                                      *
 *                                                                       *
 * Funding: National Science Foundation                                  *
 *          Zumberge Research and Innovation Fund at USC                 *
 *                                                                       *
 * This library is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this library in the file LICENSE.txt                    *
 *                                                                       *
 * This library is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

/*
  Constructs a cloth model from an obj mesh.
*/

#ifndef _CLOTHBWFROMOBJMESH_H_
#define _CLOTHBWFROMOBJMESH_H_

#include "clothBW.h"
#include "clothBWMT.h"
#include "objMesh.h"

class ClothBWFromObjMesh
{
public:
  
  // generate a cloth model from the given mesh (builds tensile, shear, and bending springs)
  // surface density and stiffnesses are the same for every triangle
  static ClothBW * GenerateClothBW(const ObjMesh * mesh, double surfaceDensity, const ClothBW::MaterialGroup & material,int addGravity=0);
  static ClothBWMT * GenerateClothBWMT(const ObjMesh * mesh, double surfaceDensity, const ClothBW::MaterialGroup & material, int numThreads = 1, int addGravity=0);

  // NOTE: materialGroup 0 is hard-coded as "default" in ObjMesh.cpp. So if a .mtl file 
  // specifies 2 materials (with 2 'usemtl' calls), there will actually be 3 material groups.
  // As a result, the density/stiffness arrays that specify a value for each material group
  // all must contain a value at the beginning for the default material (even if no vertices
  // belong to this default group).
  
  // generate a cloth model from the given mesh (builds tensile, shear, and bending springs)
  // user passes array of doubles to specify surface densities and stiffness values for each material group
  static ClothBW * GenerateClothBW(const ObjMesh * mesh, int numMaterialGroups, const double * groupSurfaceDensities, const ClothBW::MaterialGroup * materials, int addGravity=0);
  static ClothBWMT * GenerateClothBWMT(const ObjMesh * mesh, int numMaterialGroups, const double * groupSurfaceDensities, const ClothBW::MaterialGroup * materials,  int numThreads = 1, int addGravity=0);
  
protected:
};

#endif

