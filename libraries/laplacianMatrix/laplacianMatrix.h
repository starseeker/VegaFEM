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

/*
   Computes the Laplacian matrix of a 3D tetrahedral mesh. The
   computed matrices are element-based (as opposed to vertex-based) : 
   they transform scalar quantities defined on the elements. 
*/

#ifndef _LAPLACIANMATRIX_H_
#define _LAPLACIANMATRIX_H_

#include "sparseMatrix.h"
#include "tetMesh.h"

class LaplacianMatrix
{
public:

  // computes the ``umbrella'' discrete Laplacian matrix L on the tets
  // dimension is T x T, where T is the number of tets
  // In each row i, the diagonal entry is the #tets neighboring to tet i. 
  // And the entries corresponding to the neighboring tet indices are set to -1.
  // if biLaplace=1, the routine computes the bilaplace matrix, i.e., L^T L
  static SparseMatrix * ComputeTetLaplacianMatrix(const TetMesh * tetMesh, int biLaplace=0);

  // volume-weighted Laplacian on tets
  static SparseMatrix * ComputeTetVolumeWeightedLaplacianMatrix(const TetMesh * tetMesh);

  // "FEM-style" Laplacian on tets, obtained by conjugating the vertex-based Laplacian with 
  // a volume-weighted matrix converting tet quantities to vertex quantities
  static SparseMatrix * ComputeTetFEMLaplacianMatrix(const TetMesh * tetMesh);

  // Extrapolates a scalar field given at a few elements, to the entire mesh, in a smooth way
  // using the Laplace operator. Mathematically:
  //
  // minimize 1/2 <L * x, x>
  // subject to S x = Sbar
  //
  // input: L, S, Sbar
  // L is an arbitrary Laplacian matrix (computed, say, using the routines above)
  // S is the selection matrix. It is given by the array "selectedTets", of length "numSelectedTets".
  // Sbar are the prescribed values at the selected tets. 
  // numThreads is the number of employed threads. Value of 0 means single-threaded computation (equivalent to value of 1).
  // output: x (x is a vector of length n, where L is n x n); input value of x is ignored
  static void ExtrapolateScalarField(const SparseMatrix * L, int numSelectedTets, const int * selectedTets, const double * Sbar, double * x, int numThreads = 0);
};

#endif

