/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 3.0                               *
 *                                                                       *
 * "forceModel" library , Copyright (C) 2007 CMU, 2009 MIT, 2016 USC     *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Jernej Barbic                                            *
 * http://www.jernejbarbic.com/code                                      *
 *                                                                       *
 * Research: Jernej Barbic, Fun Shing Sin, Daniel Schroeder,             *
 *           Doug L. James, Jovan Popovic                                *
 *                                                                       *
 * Funding: National Science Foundation, Link Foundation,                *
 *          Singapore-MIT GAMBIT Game Lab,                               *
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
  Abstract class for f in Mu'' + Du' + f = f_ext .
  Serves as a connecting class between integrators and classes to calculate internal forces and tangent stiffness matrices.
*/

#ifndef _FORCEMODEL_H_
#define _FORCEMODEL_H_

#include <stdlib.h>
#include "sparseMatrix.h"

class ForceModel
{
public:
  virtual ~ForceModel();

  inline int Getr() { return r; }

  virtual double GetElasticEnergy(double * u) { return 0.0; }
  virtual void GetInternalForce(double * u, double * internalForces) = 0;
  virtual void GetTangentStiffnessMatrixTopology(SparseMatrix ** tangentStiffnessMatrix) = 0;
  virtual void GetTangentStiffnessMatrix(double * u, SparseMatrix * tangentStiffnessMatrix) = 0; 

  // sometimes computation time can be saved if we know that we will need both internal forces and tangent stiffness matrices:
  virtual void GetForceAndMatrix(double * u, double * internalForces, SparseMatrix * tangentStiffnessMatrix); 

  // reset routines
  virtual void ResetToZero() {}
  virtual void Reset(double * q) {}

protected:
  int r;
};

#endif

