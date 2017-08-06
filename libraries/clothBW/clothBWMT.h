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
 Multi-threaded version of the ClothBW class. 
 It uses the POSIX threads ("pthreads").
*/

#ifndef _CLOTHBWMT_H_
#define _CLOTHBWMT_H_

#include "clothBW.h"

class ClothBWMT : public ClothBW
{
public:
  
  // constructor that does not require triangleUVs input (it computes UVs automatically; note: the UVs are continuous only within each triangle; the UV map is not global (which is fine, provided one does not want to simulate anisotropic effects) )
  ClothBWMT(int numVertices, const double * restPositions, const double * masses,
      int numTriangles, const int * triangles, const int * triangleGroups,
      int numMaterialGroups, const MaterialGroup * materialGroups, int addGravity=0, int numThreads=1);
  
  // constructor with triangleUVs input
  ClothBWMT(int numVertices, const double * restPositions, const double * masses,
      int numTriangles, const int * triangles,  const double * triangleUVs, const int * triangleGroups,
      int numMaterialGroups, const MaterialGroup * materialGroups, int addGravity=0, int numThreads=1);

  // constructor from regular clothBW and thread count
  ClothBWMT(ClothBW & clothBW, int numThreads);

  // destructor
  virtual ~ClothBWMT();
  
  // multi-threaded computations

  // compute the energy, under the deformation u
  virtual double ComputeEnergy(const double * u);

  // compute the internal elastic force, under deformation u
  // note: the force has the same sign as an external force acting on the body (opposite sign as in the StVK class)
  virtual void ComputeForce(const double * u, double * f, bool addForce=false); // if addForce is "true", f will be not be reset to zero prior to adding the forces
  
  virtual void ComputeStiffnessMatrix(const double * u, SparseMatrix * K, bool addMatrix=false);

  virtual void ComputeForceAndMatrix(const double * u, double * f, SparseMatrix * K, bool addQuantity = false);

  // compute the damping force
  // unimplemented
  // you can use the damping available in the integrator class
  // virtual void ComputeDampingForce(const double * u, double * uvel, double * f, bool addForce=false);

protected:
  static void * ClothBWMT_WorkerThread(void * arg);

  int numThreads;
  
  std::vector<int> startTriangle;
  std::vector<int> endTriangle;
  
  std::vector<int> startQuad;
  std::vector<int> endQuad;
  
  std::vector<double> internalForceBuffer;
  std::vector<SparseMatrix *> sparseMatrixBuffer;
  std::vector<double> energyBuffer;
  
  void Initialize();

  void ComputeHelper(const double * u, double * energy, double * f, SparseMatrix * stiffnessMatrix, bool addQuantity);
  
};

#endif

