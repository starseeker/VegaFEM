#ifndef _MODALMATRIXTRANSPOSED_H_
#define _MODALMATRIXTRANSPOSED_H_

/*
* Copyright (c) 2007, Carnegie Mellon University
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * Neither the name of Carnegie Mellon University, nor the
*       names of its contributors may be used to endorse or promote products
*       derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

  Code author: Jernej Barbic
  Research: Jernej Barbic, Doug L. James
  Funding: NSF, Link Foundation
*/

#include "matrixMacros.h"

template<class real>
class ModalMatrixTransposed 
{
public:
  
  // n is num vertices, matrix U is 3n x r
  ModalMatrixTransposed(int n, int r, real * U, int align=0); // makes an internal copy of UT
  ~ModalMatrixTransposed();

  inline int Getr() { return r; }
  inline int Getn() { return n; }

  inline real * GetMatrix() { return UT; }

  // computes: vreduced = U^T * v
  void ProjectVector(real * v, real * vreduced);
  void ProjectSparseVector(int numSparseEntries, real * sparseVector, int * sparseVectorIndices, real * vreduced);
  void ProjectSingleVertex(int vertex, real vx, real vy, real vz, real * vreduced);

  // vreduced += U^T * v
  void AddProjectVector(real * v, real * vreduced);
  void AddProjectSparseVector(int numSparseEntries, real * sparseVector, int * sparseVectorIndices, real * vreduced);
  void AddProjectSingleVertex(int vertex, real vx, real vy, real vz, real * vreduced);

  // computes u = U * q
  void AssembleVector(real * q, real * u);
  inline void AssembleSingleVertex(int vertex, real * q, real * ux, real * uy, real * uz);

  #if defined(WIN32) || defined(__GNUC__)
    void AssembleSingleVertex_SSE(int vertex, float * q, float * outputDefo); // q and outputDefo must be 16-byte aligned; function only works for real=float
  #endif

  inline void getMatrixEntry(int vertex, int mode, real * output);

  // u += U * q
  void AddAssembleVector(real * q, real * u);
  inline void AddAssembleSingleVertex(int vertex, real * q, real * ux, real * uy, real * uz);

protected:

  real * UT; // pointer to the deformation basis
  real * UTBuffer;

  int r; // number of columns
  int n; // number of vertices
};

template <class real>
inline void ModalMatrixTransposed<real>::getMatrixEntry(int vertex, int mode, real * output)
{
  output[0] = UT[ELT(r, mode, 3*vertex+0)];
  output[1] = UT[ELT(r, mode, 3*vertex+1)];
  output[2] = UT[ELT(r, mode, 3*vertex+2)];
}

// constructs the deformation of vertex 'vertex', given q
// result goes into (ux,uy,uz)
template <class real>
inline void ModalMatrixTransposed<real>::AssembleSingleVertex
  (int vertex, real * q, real * ux, real * uy, real * uz)
{
  real * posx = &UT[ELT(r,0,3*vertex)];
  real * posy = posx + r;
  real * posz = posy + r;

  real regx = 0;
  real regy = 0;
  real regz = 0;

  int j;
  for (j=0; j<r; j+=4) // over all columns of U
  {
    real q0 = q[j];
    real q1 = q[j+1];
    real q2 = q[j+2];
    real q3 = q[j+3];
    regx += posx[j] * q0 + posx[j+1] * q1 + posx[j+2] * q2 + posx[j+3] * q3;
    regy += posy[j] * q0 + posy[j+1] * q1 + posy[j+2] * q2 + posy[j+3] * q3;
    regz += posz[j] * q0 + posz[j+1] * q1 + posz[j+2] * q2 + posz[j+3] * q3;
  }

  for(; j<r; j++)
  {
    real q0 = q[j];
    regx += posx[j] * q0;
    regy += posy[j] * q0;
    regz += posz[j] * q0;
  }

  *ux = regx;
  *uy = regy;
  *uz = regz;
}

template <class real>
inline void ModalMatrixTransposed<real>::AddAssembleSingleVertex
  (int vertex, real * q, real * ux, real * uy, real * uz)
{
  real * posx = &UT[ELT(r,0,3*vertex)];
  real * posy = posx + r;
  real * posz = posy + r;

  real regx = 0;
  real regy = 0;
  real regz = 0;

  int j;
  for (j=0; j<r; j+=4) // over all columns of U
  {
    real q0 = q[j];
    real q1 = q[j+1];
    real q2 = q[j+2];
    real q3 = q[j+3];
    regx += posx[j] * q0 + posx[j+1] * q1 + posx[j+2] * q2 + posx[j+3] * q3;
	  regy += posy[j] * q0 + posy[j+1] * q1 + posy[j+2] * q2 + posy[j+3] * q3;
	  regz += posz[j] * q0 + posz[j+1] * q1 + posz[j+2] * q2 + posz[j+3] * q3;
  }

  for(; j<r; j++)
  {
    real q0 = q[j];
    regx += posx[j] * q0;
    regy += posy[j] * q0;
    regz += posz[j] * q0;
  }

  *ux += regx;
  *uy += regy;
  *uz += regz;
}

#endif

