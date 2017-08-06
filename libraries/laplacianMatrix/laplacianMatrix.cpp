/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 3.0                               *
 *                                                                       *
 * "laplacianMatrix" library , Copyright (C) 2016 USC                    *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Hongyi Xu, Jernej Barbic                                *
 * http://www.jernejbarbic.com/code                                      *
 *                                                                       *
 * Research: Hongyi Xu, Jernej Barbic                                    *
 *                                                                       *
 * Funding: National Science Foundation                                  *
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

#include <vector>
#include <set>
using namespace std;
#include "laplacianMatrix.h"
#include "LagrangeMultiplierSolver.h"
#include "generateGradientMatrix.h"

// standard Laplace matrix on the tets
// dimension is T x T, T is #tets
// In each row i, the diagonal entry is the #tets neighboring to tet i. And the entries corresponding to the neighboring tet indices are set to -1.
SparseMatrix * LaplacianMatrix::ComputeTetLaplacianMatrix(const TetMesh * tetMesh, int biLaplace)
{
  int numElements = tetMesh->getNumElements();
  int numElementVertices = tetMesh->getNumElementVertices();
  int numVertices = tetMesh->getNumVertices();

  SparseMatrixOutline * LaplacianMatrixOutline = new SparseMatrixOutline(numElements);

  // build elements that neighbor each vertex
  vector<set<int> > vertexNeighbors(numVertices);

  for(int el = 0; el < numElements; el++)
  {
    for(int vtxIdx=0; vtxIdx<numElementVertices; vtxIdx++)
    {
      int vertexIndex = tetMesh->getVertexIndex(el, vtxIdx);
      vertexNeighbors[vertexIndex].insert(el);
    }
  }

  // build elements that neighbor each element, and assemble L
  for(int el=0; el<numElements; el++)
  {
    set<int> elementNeighbors;
    for(int vtxIdx=0; vtxIdx<numElementVertices; vtxIdx++)
    {
      int vertexIndex = tetMesh->getVertexIndex(el, vtxIdx);
      elementNeighbors.insert(vertexNeighbors[vertexIndex].begin(), vertexNeighbors[vertexIndex].end());
    }

    double counter = 0.0;
    for (set<int>::iterator it=elementNeighbors.begin(); it!=elementNeighbors.end(); it++)
    {
      if ((*it) != el)
      {
        LaplacianMatrixOutline->AddEntry(el, (*it), -1.0);
        counter++;
      }
    }

    LaplacianMatrixOutline->AddEntry(el, el, counter);
  }

  SparseMatrix * L;
  
  if (!biLaplace)
  {
    // compute standard Laplace
    L = new SparseMatrix(LaplacianMatrixOutline);
  }
  else 
  {
    // compute biLaplace
    SparseMatrix * LaplacianTetMatrix = new SparseMatrix(LaplacianMatrixOutline);
    SparseMatrix * Identity = SparseMatrix::CreateIdentityMatrix(numElements);
    L = new SparseMatrix(Identity->ConjugateMatrix(*LaplacianTetMatrix, 0, numElements));
    delete(LaplacianTetMatrix);
    delete(Identity);
  }

  delete(LaplacianMatrixOutline);

  return L;
}

SparseMatrix * LaplacianMatrix::ComputeTetVolumeWeightedLaplacianMatrix(const TetMesh * tetMesh)
{
  int numElements = tetMesh->getNumElements();
  int numElementVertices = tetMesh->getNumElementVertices();
  int numVertices = tetMesh->getNumVertices();

  SparseMatrixOutline * LaplacianMatrixOutline = new SparseMatrixOutline(numElements);

  // build elements that neighbor each vertex
  vector<set<int> > vertexNeighbors(numVertices);

  for(int el = 0; el < numElements; el++)
  {
    for(int vtxIdx=0; vtxIdx<numElementVertices; vtxIdx++)
    {
      int vertexIndex = tetMesh->getVertexIndex(el, vtxIdx);
      vertexNeighbors[vertexIndex].insert(el);
    }
  }

  double * vols = (double*) malloc (sizeof(double) * numElements);
  for(int el=0; el<numElements; el++)
    vols[el] = tetMesh->getElementVolume(el);

  // build elements that neighbor each element, and assemble L
  for(int el=0; el<numElements; el++)
  {
    set<int> elementNeighbors;
    for(int vtxIdx=0; vtxIdx<numElementVertices; vtxIdx++)
    {
      int vertexIndex = tetMesh->getVertexIndex(el, vtxIdx);
      elementNeighbors.insert(vertexNeighbors[vertexIndex].begin(), vertexNeighbors[vertexIndex].end());
    }
 
    double denom = 0.0;
    for (set<int>::iterator it=elementNeighbors.begin(); it!=elementNeighbors.end(); it++)
    {
      if ((*it) != el)
      {
        LaplacianMatrixOutline->AddEntry(el, (*it), -(vols[el] + vols[*it]));
        denom += (vols[el] + vols[*it]);
      }
    }

    LaplacianMatrixOutline->AddEntry(el, el, denom);
  }

  SparseMatrix * L = new SparseMatrix(LaplacianMatrixOutline);
    
  free(vols);
  delete(LaplacianMatrixOutline);

  return L;
}


SparseMatrix * LaplacianMatrix::ComputeTetFEMLaplacianMatrix(const TetMesh * tetMesh)
{
  int numVertices = tetMesh->getNumVertices();
  int numElements = tetMesh->getNumElements();
 
  double * sqrtvol = (double*) malloc (sizeof(double) * numElements);
  for(int el=0; el<numElements; el++)
    sqrtvol[el] = sqrt(tetMesh->getElementVolume(el));
 
  SparseMatrix * LaplacianVertexMatrix;
  GenerateGradientMatrix::GenerateForScalarField(tetMesh, &LaplacianVertexMatrix, sqrtvol);

  SparseMatrixOutline * smoothMatrixOutline = new SparseMatrixOutline(numElements);

  for(int i=0; i<numVertices; i++)
  {
    vector<int> elementVertices;
    elementVertices.push_back(i);

    vector<int> elementNeighbors;
    tetMesh->getElementsTouchingVertices(elementVertices, elementNeighbors);

    double volT = 0;
    for(unsigned int el = 0; el < elementNeighbors.size(); el++)
      volT += sqrtvol[elementNeighbors[el]] * sqrtvol[elementNeighbors[el]];

    for(unsigned int el = 0; el < elementNeighbors.size(); el++)
      smoothMatrixOutline->AddEntry(i, elementNeighbors[el], sqrtvol[elementNeighbors[el]] * sqrtvol[elementNeighbors[el]] / volT);
  }

  SparseMatrix * smoothMatrix = new SparseMatrix(smoothMatrixOutline);
  SparseMatrix * L = new SparseMatrix(LaplacianVertexMatrix->ConjugateMatrix(*smoothMatrix, 0, numElements));

  delete(LaplacianVertexMatrix);
  delete(smoothMatrixOutline);
  delete(smoothMatrix);
  free(sqrtvol);

  return L;
}

void LaplacianMatrix::ExtrapolateScalarField(const SparseMatrix * L, int numSelectedTets, const int * selectedTets, const double * Sbar, double * x, int numThreads)
{
  int numElements = L->GetNumRows();

  if (numSelectedTets == 0)
  {
    printf("Warning: no selected tets in \"ExtrapolateScalarField\". Output x is not modified.\n");
    return;
  }

  double * xExtended = (double*) malloc (sizeof(double) * (numElements + numSelectedTets));
  double * rhs = (double*) malloc (sizeof(double)  * (numSelectedTets + numElements));
  memset(rhs, 0, sizeof(double) * numElements);

  SparseMatrixOutline * JOutline = new SparseMatrixOutline(numSelectedTets);
  for(int i=0; i<numSelectedTets; i++)
  {
    JOutline->AddEntry(i, selectedTets[i], 1.0);
    rhs[numElements + i] = Sbar[i];
  }
  SparseMatrix * J = new SparseMatrix(JOutline);
  LagrangeMultiplierSolver * solver = new LagrangeMultiplierSolver(L, J, NULL, 0, NULL, numThreads);
  solver->SolveLinearSystem(xExtended, rhs);

  memcpy(x, xExtended, sizeof(double) * numElements);

  delete(solver);
  delete(J);
  delete(JOutline);
  free(rhs);
  free(xExtended);
}

