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

#include <stdlib.h>
#include <stdio.h>
#include "lapack-headers.h"
#ifdef __APPLE__
  #include "TargetConditionals.h"
#endif
#include "modalMatrixTransposed.h"

template<bool C>
class _cblas_xgemv {};

template<>
class _cblas_xgemv<true> {
public:
    inline static void f (const  CBLAS_ORDER order,
                 const  CBLAS_TRANSPOSE trans, const int m, const int n,
                 const float alpha, const float * a, const int lda,
                 const float * x, const int incx, float beta,
                 float * y, const int incy)
     { cblas_sgemv(order, trans, m, n, alpha, a, lda, x, incx, beta, y, incy);}
};

template<>
class _cblas_xgemv<false> {
public:
    inline static void f (const  CBLAS_ORDER order,
                 const  CBLAS_TRANSPOSE trans, const int m, const int n,
                 const double alpha, const double * a, const int lda,
                 const double * x, const int incx, double beta,
                 double * y, const int incy)
     { cblas_dgemv(order, trans, m, n, alpha, a, lda, x, incx, beta, y, incy);}
};

template <class real>
ModalMatrixTransposed<real>::ModalMatrixTransposed(int n, int r, real * U, int align)
{
  this->n = n;

  if (align)
  {
    this->r = 4 * (((r - 1) / 4) + 1); 
  }
  else
    this->r = r;

  unsigned int blockSize = sizeof(real) * 3 * n * this->r;
  if (align)
    blockSize += 16;

  UTBuffer = (real*) malloc (blockSize);

  if (align)
  {
    // align UT to 16 bytes
    /*char * memPos = (char *) UTBuffer;
    unsigned int memPosValue = *(unsigned int*)(&memPos);
    memPosValue = (memPosValue & 0xFFFFFFF0) + 16;
    UT = *(real**)(&memPosValue);*/
    size_t alignedAddress = ((size_t)((unsigned char*)UTBuffer + sizeof(void*) - 1)) &
	    (~(16 - 1));
    alignedAddress += 16;
    UT = (real*) alignedAddress;
  }
  else
    UT = UTBuffer;

  for(int i=0; i<3*n; i++)
  {
    for(int j=0; j<r; j++)
    {
      UT[ELT(this->r,j,i)] = U[ELT(3*n,i,j)];
    }
    for(int j=r; j<this->r; j++)
      UT[ELT(this->r,j,i)] = 0;
  }
}

template <class real>
ModalMatrixTransposed<real>::~ModalMatrixTransposed()
{
  free(UTBuffer);
}

template <class real>
void ModalMatrixTransposed<real>::ProjectSingleVertex
  (int vertex, real vx, real vy, real vz, real * vreduced) 
{
  int j;
  for (j=0; j<r; j++) // over all columns of U
  {
    vreduced[j] = UT[ELT(r,j,3*vertex+0)] * vx +
           UT[ELT(r,j,3*vertex+1)] * vy +
           UT[ELT(r,j,3*vertex+2)] * vz;
  }
}

template <class real>
void ModalMatrixTransposed<real>::AddProjectSingleVertex
  (int vertex, real vx, real vy, real vz, real * vreduced) 
{
  int j;
  for (j=0; j<r; j++) // over all columns of U
  {
  	vreduced[j] += UT[ELT(r,j,3*vertex+0)] * vx +
             UT[ELT(r,j,3*vertex+1)] * vy +
             UT[ELT(r,j,3*vertex+2)] * vz;
  }
}

template <class real>
void ModalMatrixTransposed<real>::ProjectVector(real * v, real * vreduced) 
{
  // has to make inner product of vector f will all the columns of U
  // i.e. multiply U^T * f = q
  CBLAS_ORDER     order= CblasColMajor;
  CBLAS_TRANSPOSE trans= CblasNoTrans;
  int M = r;
  int N = 3*n;
  real alpha = 1;
  real * a = UT;
  int lda = r;
  real * x = v;
  int incx = 1;
  real beta = 0;
  real * y = vreduced; 
  int incy = 1;

  _cblas_xgemv<sizeof(real)==sizeof(float)>::f(order, trans, M, N, alpha, a, lda, x, incx, beta, y, incy);
  //cblas_dgemv(order, trans, M, N, alpha, a, lda, x, incx, beta, y, incy);
}

template <class real>
void ModalMatrixTransposed<real>::AddProjectVector(real * v, real * vreduced) 
{
  // has to make inner product of vector f will all the columns of U
  // i.e. multiply U^T * f = q
  CBLAS_ORDER     order= CblasColMajor;
  CBLAS_TRANSPOSE trans= CblasNoTrans;
  int M = r;
  int N = 3*n;
  real alpha = 1;
  real * a = UT;
  int lda = r;
  real * x = v;
  int incx = 1;
  real beta = 1;
  real * y = vreduced; 
  int incy = 1;

  _cblas_xgemv<sizeof(real)==sizeof(float)>::f(order, trans, M, N, alpha, a, lda, x, incx, beta, y, incy);
  //cblas_dgemv(order, trans, M, N, alpha, a, lda, x, incx, beta, y, incy);
}

template <class real>
void ModalMatrixTransposed<real>::ProjectSparseVector(int numSparseEntries, real * sparseVector, int * sparseVectorIndices, real * vreduced)
{
  // has to make inner product of vector f will all the columns of U
  int i,j;
  for (j=0; j<r; j++) // over all columns of U
  {
    // dot product of column j of U with vector f
    vreduced[j] = 0;
    for (i=0; i<numSparseEntries; i++)
      vreduced[j] += UT[ELT(r,j,sparseVectorIndices[i])] * sparseVector[i];
  }
}

template <class real>
void ModalMatrixTransposed<real>::AddProjectSparseVector(int numSparseEntries, real * sparseVector, int * sparseVectorIndices, real * vreduced)
{
  // has to make inner product of vector f will all the columns of U
  int i,j;
  for (j=0; j<r; j++) // over all columns of U
  {
    // dot product of column j of U with vector f
    for (i=0; i<numSparseEntries; i++)
      vreduced[j] += UT[ELT(r,j,sparseVectorIndices[i])] * sparseVector[i];
  }
}

template <class real>
void ModalMatrixTransposed<real>::AssembleVector(real * q, real * u) //u = U * q;
{
  CBLAS_ORDER     order= CblasColMajor;
  CBLAS_TRANSPOSE trans= CblasTrans;
  int M = r;
  int N = 3*n;
  real alpha = 1;
  real * a = UT;
  int lda = r;
  real * x = q;
  int incx = 1;
  real beta = 0;
  real * y = u; 
  int incy = 1;

  _cblas_xgemv<sizeof(real)==sizeof(float)>::f(order, trans, M, N, alpha, a, lda, x, incx, beta, y, incy);
  //cblas_dgemv(order, trans, M, N, alpha, a, lda, x, incx, beta, y, incy);
}

template <class real>
void ModalMatrixTransposed<real>::AddAssembleVector(real * q, real * u) //u = U * q;
{
  CBLAS_ORDER     order= CblasColMajor;
  CBLAS_TRANSPOSE trans= CblasTrans;
  int M = r;
  int N = 3*n;
  real alpha = 1;
  real * a = UT;
  int lda = r;
  real * x = q;
  int incx = 1;
  real beta = 1.0;
  real * y = u; 
  int incy = 1;

  _cblas_xgemv<sizeof(real)==sizeof(float)>::f(order, trans, M, N, alpha, a, lda, x, incx, beta, y, incy);
  //cblas_dgemv(order, trans, M, N, alpha, a, lda, x, incx, beta, y, incy);
}

#ifdef WIN32
template <>
void ModalMatrixTransposed<float>::AssembleSingleVertex_SSE
  (int vertex, float * q, float * outputDefo) // q and outputDefo must be 16-byte aligned
{
  long int rr = this->r;

  float * posx = &UT[ELT(rr,0,3*vertex)];
  float * posy = posx + rr;
  float * posz = posy + rr;
/*
  q[0] = 1.55;
  q[1] = -2.34;
  q[2] = 3.45;
  q[3] = 4.56;
  q[4] = 0.45;
  q[5] = 1.34;
  q[9] = 12.11;
  q[12] = 2.4;
  q[13] = 4.5;

  float bla;
*/

  #ifdef _M_IX86
  __asm{
    XOR ebx, ebx
    XORPS xmm0, xmm0
    XORPS xmm1, xmm1
    XORPS xmm2, xmm2

   loopStart: 
 
    MOV edx, ebx
    ADD edx, 4
    CMP edx, rr
    JG  outOfLoop

    // load q into xmm7
    MOV    edx, [q]
    MOVAPS xmm7, [edx]

    MOV    edx, posx
    MOVAPS xmm3, [edx] 
    MULPS  xmm3, xmm7
    ADDPS  xmm0, xmm3

    MOV    edx, posy
    MOVAPS xmm4, [edx] 
    MULPS  xmm4, xmm7
    ADDPS  xmm1, xmm4

    MOV    edx, posz
    MOVAPS xmm5, [edx] 
    MULPS  xmm5, xmm7
    ADDPS  xmm2, xmm5

    ADD    posx, 16
    ADD    posy, 16
    ADD    posz, 16
    ADD    q, 16
    ADD    ebx, 4
    JMP    loopStart

   outOfLoop:

    // make zero artifical fourth vector
    XORPS  xmm3, xmm3

    // must add entries in xmm0,1,2,3 store sums to u
    MOVAPS xmm4, xmm0
    SHUFPS xmm4, xmm1, 0x44
    MOVAPS xmm5, xmm0
    SHUFPS xmm5, xmm1, 0xEE
    ADDPS  xmm4, xmm5

    MOVAPS xmm6, xmm2
    SHUFPS xmm6, xmm3, 0x44
    MOVAPS xmm7, xmm2
    SHUFPS xmm7, xmm3, 0xEE
    ADDPS  xmm6, xmm7

    // now, xmm4 and xmm6 are set
    MOVAPS xmm5, xmm4
    MOVAPS xmm7, xmm4
    SHUFPS xmm5, xmm6, 0x88
    SHUFPS xmm7, xmm6, 0xDD
    ADDPS  xmm5, xmm7

    MOV    edx, [outputDefo]
    MOVAPS [edx], xmm5
  };
  #endif
}
#endif


#if ((defined __GNUC__) && (!TARGET_CPU_PPC))
template <>
void ModalMatrixTransposed<float>::AssembleSingleVertex_SSE
  (int vertex, float * q, float * outputDefo) // q and outputDefo must be 16-byte aligned
{
  long int rr = sizeof(float) * r;
  float * posx = &UT[r*3*vertex];
  long int i = 0;

  asm(
    //"xor %3, %3\n\t"
    "xorps %%xmm0, %%xmm0\n\t"
    "xorps %%xmm1, %%xmm1\n\t"
    "xorps %%xmm2, %%xmm2\n\t"
    "\n\t"

  "loopStart:\n\t" 
"\n\t" 
    "add $16, %3\n\t"
    "cmp %1, %3\n\t"
    "jg  outOfLoop\n\t"
"\n\t"
    "// load q into xmm7\n\t"
    "movaps (%2), %%xmm7\n\t"
"\n\t"
    "movaps (%0), %%xmm3\n\t" 
    "mulps  %%xmm7, %%xmm3\n\t"
    "addps  %%xmm3, %%xmm0\n\t"
"\n\t"
    "movaps (%0,%1,1), %%xmm4\n\t" 
    "mulps  %%xmm7, %%xmm4\n\t"
    "addps  %%xmm4, %%xmm1\n\t"
"\n\t"
    "movaps (%0,%1,2), %%xmm5\n\t" 
    "mulps  %%xmm7, %%xmm5\n\t"
    "addps  %%xmm5, %%xmm2\n\t"
"\n\t"
    "add    $16, %0\n\t"
    "add    $16, %2\n\t"
    "jmp    loopStart\n\t"
"\n\t"
   "outOfLoop:\n\t"
"\n\t"
    "// make zero artifical fourth vector\n\t"
    "xorps  %%xmm3, %%xmm3\n\t"
"\n\t"
    "// must add entries in xmm0,1,2,3 store sums to u\n\t"
    "movaps %%xmm0, %%xmm4\n\t"
    "shufps $0x44, %%xmm1, %%xmm4\n\t"
    "movaps %%xmm0, %%xmm5\n\t"
    "shufps $0xEE, %%xmm1, %%xmm5\n\t"
    "addps  %%xmm5, %%xmm4\n\t"
"\n\t"
    "movaps %%xmm2, %%xmm6\n\t"
    "shufps $0x44, %%xmm3, %%xmm6\n\t"
    "movaps %%xmm2, %%xmm7\n\t"
    "shufps $0xEE, %%xmm3, %%xmm7\n\t"
    "addps  %%xmm6, %%xmm7\n\t"
"\n\t"
    "// now, xmm4 and xmm6 are set\n\t"
    "movaps %%xmm4, %%xmm5\n\t"
    "movaps %%xmm4, %%xmm7\n\t"
    "shufps $0x88, %%xmm6, %%xmm5\n\t"
    "shufps $0xDD, %%xmm6, %%xmm7\n\t"
    "addps  %%xmm7, %%xmm5\n\t"
"\n\t"
    "movaps %%xmm5, (%4)\n\t"
    : /* no output registers */
    : "r" (posx), "r" (rr), "r" (q), "r" (i), "r" (outputDefo)
    : "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
      "%xmm4", "%xmm5", "%xmm6", "%xmm7", "memory"
   );
}
#endif

template ModalMatrixTransposed<double>::ModalMatrixTransposed(int n, int r, double * U, int align);
template ModalMatrixTransposed<float>::ModalMatrixTransposed(int n, int r, float * U, int align);

template ModalMatrixTransposed<double>::~ModalMatrixTransposed();
template ModalMatrixTransposed<float>::~ModalMatrixTransposed();

template void ModalMatrixTransposed<double>::ProjectSingleVertex
  (int vertex, double vx, double vy, double vz, double * vreduced);
template void ModalMatrixTransposed<float>::ProjectSingleVertex
  (int vertex, float vx, float vy, float vz, float * vreduced);

template void ModalMatrixTransposed<double>::AddProjectSingleVertex
  (int vertex, double vx, double vy, double vz, double * vreduced);
template void ModalMatrixTransposed<float>::AddProjectSingleVertex
  (int vertex, float vx, float vy, float vz, float * vreduced);

template void ModalMatrixTransposed<double>::ProjectVector(double * v, double * vreduced);
template void ModalMatrixTransposed<float>::ProjectVector(float * v, float * vreduced);

template void ModalMatrixTransposed<double>::AddProjectVector(double * v, double * vreduced);
template void ModalMatrixTransposed<float>::AddProjectVector(float * v, float * vreduced);

template void ModalMatrixTransposed<double>::ProjectSparseVector
  (int numSparseEntries, double * sparseVector, 
	 int * sparseVectorIndices, double * vreduced);
template void ModalMatrixTransposed<float>::ProjectSparseVector
  (int numSparseEntries, float * sparseVector, 
	 int * sparseVectorIndices, float * vreduced);

template void ModalMatrixTransposed<double>::AddProjectSparseVector
  (int numSparseEntries, double * sparseVector, 
	 int * sparseVectorIndices, double * vreduced);
template void ModalMatrixTransposed<float>::AddProjectSparseVector
  (int numSparseEntries, float * sparseVector, 
	 int * sparseVectorIndices, float * vreduced);

template void ModalMatrixTransposed<double>::AssembleVector(double * q, double * u);
template void ModalMatrixTransposed<float>::AssembleVector(float * q, float * u);

template void ModalMatrixTransposed<double>::AddAssembleVector(double * q, double * u);
template void ModalMatrixTransposed<float>::AddAssembleVector(float * q, float * u);
