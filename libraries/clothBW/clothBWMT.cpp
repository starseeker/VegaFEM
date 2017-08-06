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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "clothBWMT.h"
#include <pthread.h>
using namespace std;

ClothBWMT::ClothBWMT(int numVertices_, const double * restPositions_, const double * masses_,
    int numTriangles_, const int * triangles_, const int * triangleGroups_,
    int numMaterialGroups_, const MaterialGroup * materialGroups_, int addGravity_, int numThreads_) :
    ClothBW(numVertices_, restPositions_, masses_, numTriangles_, triangles_, triangleGroups_,
        numMaterialGroups_, materialGroups_, addGravity_), numThreads(numThreads_)
{
  Initialize();
}

// constructor with triangleUVs input
ClothBWMT::ClothBWMT(int numVertices_, const double * restPositions_, const double * masses_,
    int numTriangles_, const int * triangles_, const double * triangleUVs_, const int * triangleGroups_,
    int numMaterialGroups_, const MaterialGroup * materialGroups_, int addGravity_, int numThreads_) :
    ClothBW(numVertices_, restPositions_, masses_, numTriangles_, triangles_, triangleUVs_, triangleGroups_,
        numMaterialGroups_, materialGroups_, addGravity_), numThreads(numThreads_)
{
  Initialize();
}

ClothBWMT::ClothBWMT(ClothBW & clothBW, int numThreads_) : ClothBW(clothBW), numThreads(numThreads_)
{
  Initialize();
}

// destructor
ClothBWMT::~ClothBWMT()
{
  for(size_t i = 0; i < sparseMatrixBuffer.size(); i++)
    delete sparseMatrixBuffer[i];
}

// data structure to hold data for each thread
struct ClothBWMT_threadArg
{
  ClothBWMT * clothBW_MT;
  const double * u;
  int rank;
  double * energy;
  double * internalForce;
  SparseMatrix * stiffnessMatrix;
};

// global function
void * ClothBWMT::ClothBWMT_WorkerThread(void * arg)
{
  // cast to struct
  struct ClothBWMT_threadArg * threadArgp = (struct ClothBWMT_threadArg*) arg;
  
  // copy info into local vars
  ClothBWMT * cloth = threadArgp->clothBW_MT;
  const double * u = threadArgp->u;
  int rank = threadArgp->rank;
  int startTriangle = cloth->startTriangle[rank];
  int endTriangle = cloth->endTriangle[rank];
  int startQuad = cloth->startQuad[rank];
  int endQuad = cloth->endQuad[rank];
  
  const bool * cond = cloth->cond;
  if (cond[0] || cond[2])
    cloth->AddStretchAndShear(u, startTriangle, endTriangle, NULL,
        (cond[0] ? threadArgp->internalForce : NULL), (cond[2] ? threadArgp->stiffnessMatrix : NULL));

  if (cond[1] || cond[3])
    cloth->AddBend(u, startQuad, endQuad, NULL,
        (cond[1] ? threadArgp->internalForce : NULL), (cond[3] ? threadArgp->stiffnessMatrix : NULL));

  return NULL;
}

// multi-threaded computations

double ClothBWMT::ComputeEnergy(const double * u)
{
  double energy = 0.0;
  ComputeHelper(u, &energy, NULL, NULL, false);
  return energy;
}
// compute the internal elastic force, under deformation u
// note: the force has the same sign as an external force acting on the body (opposite sign as in the StVK class)
// if addForce is "true", f will be not be reset to zero prior to adding the forces
void ClothBWMT::ComputeForce(const double * u, double * f, bool addForce)
{
  ComputeHelper(u, NULL, f, NULL, addForce);
}

//// compute the damping force
//// unimplemented
//void ClothBWMT::ComputeDampingForce(const double * u, double * uvel, double * f, bool addForce)
//{
//  
//}

void ClothBWMT::ComputeStiffnessMatrix(const double * u, SparseMatrix * K, bool addMatrix)
{
  ComputeHelper(u, NULL, NULL, K, addMatrix);
}

void ClothBWMT::ComputeForceAndMatrix(const double * u, double * f, SparseMatrix * K, bool addQuantity)
{
  ComputeHelper(u, NULL, f, K, addQuantity);
}

void ClothBWMT::Initialize()
{
  // size of our force buffer (3 x numVertices for each thread)
  internalForceBuffer.resize(numThreads * 3 * numVertices);
  energyBuffer.resize(numThreads);
  
  // generate skeleton matrices
  sparseMatrixBuffer.resize(numThreads, NULL);
  SparseMatrix * sparseMatrix = NULL;
  GenerateStiffnessMatrixTopology(&sparseMatrix);
  for(int i=0; i<numThreads; i++)
    sparseMatrixBuffer[i] = new SparseMatrix(*sparseMatrix);
  delete sparseMatrix;
  
  // split the workload
  startTriangle.resize(numThreads);
  endTriangle.resize(numThreads);
  
  startQuad.resize(numThreads);
  endQuad.resize(numThreads);
  
  int remainderTriangles = numTriangles % numThreads;
  int remainderQuads = numQuads % numThreads;
  
  // the first 'remainder' nodes will process one triangle + one quad more
  
  int jobSizeTriangles = numTriangles / numThreads;
  int jobSizeQuads = numQuads / numThreads;
  
  for(int rank=0; rank < numThreads; rank++)
  {
    // triangles
    if (rank < remainderTriangles)
    {
      startTriangle[rank] = rank * (jobSizeTriangles+1);
      endTriangle[rank] = (rank+1) * (jobSizeTriangles+1);
    }
    else
    {
      startTriangle[rank] = remainderTriangles * (jobSizeTriangles+1) + (rank-remainderTriangles) * jobSizeTriangles;
      endTriangle[rank] = remainderTriangles * (jobSizeTriangles+1) + ((rank-remainderTriangles)+1) * jobSizeTriangles;
    }
    
    // quads
    if (rank < remainderQuads)
    {
      startQuad[rank] = rank * (jobSizeQuads+1);
      endQuad[rank] = (rank+1) * (jobSizeQuads+1);
    }
    else
    {
      startQuad[rank] = remainderQuads * (jobSizeQuads+1) + (rank-remainderQuads) * jobSizeQuads;
      endQuad[rank] = remainderQuads * (jobSizeQuads+1) + ((rank-remainderQuads)+1) * jobSizeQuads;
    } 
  }
  
  printf("Total triangles: %d \n",numTriangles);
  printf("Total quads: %d \n",numQuads);
  
  printf("Num threads: %d \n",numThreads);
  printf("Canonical job size (triangles): %d \n",jobSizeTriangles);
  printf("Canonical job size (quads): %d \n",jobSizeQuads);
  
  printf("Num threads with job size augmented by one triangle: %d \n",remainderTriangles);
  printf("Num threads with job size augmented by one quad: %d \n",remainderQuads);
}

void ClothBWMT::ComputeHelper(const double * u, double * energy, double * internalForce, SparseMatrix * stiffnessMatrix, bool addQuantity)
{
  // launch threads
  int numVertices3 = 3*numVertices;
  vector<ClothBWMT_threadArg> threadArgv(numThreads);
  vector<pthread_t> tid(numThreads);
  
  for(int i=0; i<numThreads; i++)
  {
    threadArgv[i].clothBW_MT = this;
    threadArgv[i].u = u;
    threadArgv[i].rank = i;
    threadArgv[i].energy = (energy ? &energyBuffer[i] : NULL);
    threadArgv[i].internalForce = (internalForce ? &internalForceBuffer[i*numVertices3] : NULL);
    threadArgv[i].stiffnessMatrix = (stiffnessMatrix ? sparseMatrixBuffer[i] : NULL);
  }
  
  // initialize all buffers
  if (energy)
    energyBuffer.assign(numThreads, 0.0);

  if (internalForce)
    internalForceBuffer.assign(numVertices3 * numThreads, 0.0);

  if (stiffnessMatrix)
    for(int i=0; i<numThreads; i++)
      sparseMatrixBuffer[i]->ResetToZero();
  
  for(int i=0; i<numThreads; i++)  
  {
    if (pthread_create(&tid[i], NULL, ClothBWMT_WorkerThread, &threadArgv[i]) != 0)
    {
      printf("Error: unable to launch thread %d.\n", i);
      throw(1);
    }
  }
  
  for(int i=0; i<numThreads; i++)
  {
    if (pthread_join(tid[i], NULL) != 0)
    {
      printf("Error: unable to join thread %d.\n", i);
      throw(1);
    }
  }
  
  // assemble results
  if (energy)
  {
    if (!addQuantity)
      *energy = 0.0;

    for(int i = 0; i < numThreads; i++)
      *energy += energyBuffer[i];
  }

  if (internalForce)
  {
    if (!addQuantity)
      memset(internalForce, 0, sizeof(double) * numVertices3);

    for(int i=0; i<numThreads; i++)
    {
      double * source = &internalForceBuffer[i * numVertices3];
      for(int j=0; j<numVertices3; j++)
        internalForce[j] += source[j];
    }

    if (addGravity)
      ComputeGravity(internalForce);
  }

  if (stiffnessMatrix)
  {
    if (!addQuantity)
      stiffnessMatrix->ResetToZero();
    for(int i=0; i<numThreads; i++)
      *stiffnessMatrix += *(sparseMatrixBuffer[i]);
  }
}

