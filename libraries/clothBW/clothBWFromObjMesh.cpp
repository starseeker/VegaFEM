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

#include "macros.h"
#include "clothBWFromObjMesh.h"
using namespace std;

ClothBW * ClothBWFromObjMesh::GenerateClothBW(const ObjMesh * mesh, double surfaceDensity, const ClothBW::MaterialGroup & material, int addGravity)
{
  int numMaterialGroups = mesh->getNumGroups();
  
  // create arrays of density & stiffness values
  vector<ClothBW::MaterialGroup> materialGroups(numMaterialGroups);
  vector<double> densities(numMaterialGroups);
  
  // fill arrays all with the same value
  for (int i=0; i<numMaterialGroups; i++) 
  {
    materialGroups[i] = material;
    densities[i] = surfaceDensity;
  }
  
  // construct clothBW
  return GenerateClothBW(mesh, numMaterialGroups, densities.data(), materialGroups.data(), addGravity);
}

ClothBWMT * ClothBWFromObjMesh::GenerateClothBWMT(const ObjMesh * mesh, double surfaceDensity, const ClothBW::MaterialGroup & material, int numThreads,  int addGravity)
{
  int numMaterialGroups = mesh->getNumGroups();

  // create arrays of density & stiffness values
  vector<ClothBW::MaterialGroup> materialGroups(numMaterialGroups);
  vector<double> densities(numMaterialGroups);

  // fill arrays all with the same value
  for (int i=0; i<numMaterialGroups; i++) 
  {
    materialGroups[i] = material;
    densities[i] = surfaceDensity;
  }

  // construct clothBW
  return GenerateClothBWMT(mesh, numMaterialGroups, densities.data(), materialGroups.data(), numThreads, addGravity);
}


ClothBW * ClothBWFromObjMesh::GenerateClothBW(const ObjMesh * mesh, int numMaterialGroups, const double * groupSurfaceDensities,
    const ClothBW::MaterialGroup * material, int addGravity)
{
  int numVertices = 0, numTriangles = 0, numGroups = 0;
  double * restPositions = NULL;
  int * triangles = NULL, * triangleGroups = NULL;
  
  //exportGeometry returns triangulated results
  mesh->exportGeometry(&numVertices, &restPositions, &numTriangles, &triangles, &numGroups, &triangleGroups);

  if (numGroups != numMaterialGroups)
  {
    printf("Mismatch in the number of groups. Mesh has %d groups.\n", numGroups);
    free(restPositions);
    free(triangles);
    free(triangleGroups);
    return NULL;
  }
  
  vector<double> groupSurfaceMassDensity(numGroups), masses;
  memcpy(groupSurfaceMassDensity.data(), groupSurfaceDensities, sizeof(double) * numGroups);
  mesh->computeMassPerVertex(groupSurfaceMassDensity, masses);

  ClothBW * clothBW = new ClothBW(numVertices,  restPositions, masses.data(), numTriangles, triangles, triangleGroups,
      numGroups, material, addGravity);
  
  free(restPositions);
  free(triangles);
  free(triangleGroups);
  return clothBW;
}

ClothBWMT * ClothBWFromObjMesh::GenerateClothBWMT(const ObjMesh * mesh, int numMaterialGroups, const double * groupSurfaceDensities,
    const ClothBW::MaterialGroup * material, int numThreads, int addGravity)
{ 
  int numVertices = 0, numTriangles = 0, numGroups = 0;
  double * restPositions = NULL;
  int * triangles = NULL, * triangleGroups = NULL;
  
  mesh->exportGeometry(&numVertices, &restPositions, &numTriangles, &triangles, &numGroups, &triangleGroups);

  if (numGroups != numMaterialGroups)
  {
    printf("Mismatch in the number of groups. Mesh has %d groups.\n", numGroups);
    free(restPositions);
    free(triangles);
    free(triangleGroups);
    return NULL;
  }

  vector<double> groupSurfaceMassDensity(numGroups), masses;
  memcpy(groupSurfaceMassDensity.data(), groupSurfaceDensities, sizeof(double) * numGroups);
  mesh->computeMassPerVertex(groupSurfaceMassDensity, masses);

  ClothBWMT * clothBW = new ClothBWMT(numVertices, restPositions, masses.data(), numTriangles, triangles, triangleGroups,
      numGroups, material, addGravity, numThreads);

  free(restPositions);
  free(triangles);
  free(triangleGroups);
  return clothBW;
}
