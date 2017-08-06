/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 3.1                               *
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
#include <math.h>
#include <iostream>
#include <map>
#include <set>
#include <cassert>
#include "clothBW.h"
#include "macros.h"
#include "matrixMultiplyMacros.h"
#include "minivector.h"
using namespace std;

// constructor without triangleUVs
ClothBW::ClothBW(int numVertices_, const double * restPositions_, const double * masses_,
    int numTriangles_, const int * triangles_, const int * triangleGroups_,
    int numMaterialGroups_, const MaterialGroup * materialGroups_, int addGravity_)
  : addGravity(addGravity_), g(9.81)
{
  // build triangle UVs
  vector<double> triangleUVs(3 * 2 * numTriangles_);
  for (int i = 0 ; i < numTriangles_; i++)
  {	
    int vertexA = triangles_[3*i+0];
    int vertexB = triangles_[3*i+1];
    int vertexC = triangles_[3*i+2];

    triangleUVs[6*i+0] = 0.0; // vertex A, u coordinate
    triangleUVs[6*i+1] = 0.0; // vertex A, v coordinate

    Vec3d x0; // vector from A to B
    x0[0] = restPositions_[3*vertexA+0] - restPositions_[3*vertexB+0];
    x0[1] = restPositions_[3*vertexA+1] - restPositions_[3*vertexB+1]; 
    x0[2] = restPositions_[3*vertexA+2] - restPositions_[3*vertexB+2];  

    double lengthAB = len(x0);
    triangleUVs[6*i+2] = lengthAB;	// vertex B, u coordinate
    triangleUVs[6*i+3] = 0.0; // vertex B, v coordinate

    Vec3d x1; // vector from A to C
    x1[0] = restPositions_[3*vertexA+0] - restPositions_[3*vertexC+0];
    x1[1] = restPositions_[3*vertexA+1] - restPositions_[3*vertexC+1]; 
    x1[2] = restPositions_[3*vertexA+2] - restPositions_[3*vertexC+2];    

    Vec3d xn0 = norm(x0); // normalized vector from A to B
    double lengthAC = len(x1);
    triangleUVs[6*i+4] = dot(xn0, x1); // vertex C, u coordinate
    triangleUVs[6*i+5] = sqrt(lengthAC * lengthAC - triangleUVs[6*i+4] * triangleUVs[6*i+4]);	// vertex C, v coordinate
  }

  GenerateBW(numVertices_, restPositions_, masses_, numTriangles_, triangles_, triangleUVs.data(), 
    triangleGroups_, numMaterialGroups_, materialGroups_, addGravity_);
}

// constructor with trianglesUVs
ClothBW::ClothBW(int numVertices_, const double * restPositions_, const double * masses_,
    int numTriangles_, const int * triangles_, const double * triangleUVs_, const int * triangleGroups_,
    int numMaterialGroups_, const MaterialGroup * materialGroups_, int addGravity_) : addGravity(addGravity_), g(9.81)
{
  GenerateBW(numVertices_, restPositions_, masses_, numTriangles_, triangles_, triangleUVs_, 
    triangleGroups_, numMaterialGroups_, materialGroups_, addGravity_);
}

ClothBW::WuvInfo ClothBW::ComputeWuvInfo(const double triangleUV[6])
{
  // distance of neighboring vertices in planar coordinates. 
  // (delta_u1, delta_v1): planar vector from A to B,
  // (delta_u2, delta_v2): planar vector from B to A.
  double du1, du2, dv1, dv2;  
  du1 = triangleUV[2] - triangleUV[0];
  dv1 = triangleUV[3] - triangleUV[1];
  du2 = triangleUV[4] - triangleUV[0];
  dv2 = triangleUV[5] - triangleUV[1];
  double Delta = 1.0/(du1*dv2-du2*dv1);
  WuvInfo info; // compute derivatives of wu and wv with respect to vtx position x
  // 3x1 vector: wu = ( (x1-x0) dv2 - (x2-x0) dv1 ) / (du1 dv2 - du2 dv1), xi is a 3x1 vector for vtx i on a triangle<x0,x1,x2>
  // 3x1 vector: wv = (-(x1-x0) du2 + (x2-x0) du1 ) / (du1 dv2 - du2 dv1)

  info.pwupx[0] = (dv1 - dv2) * Delta;
  info.pwupx[1] = dv2 * Delta;
  info.pwupx[2] = -dv1 * Delta;
  info.pwvpx[0] = (du2 - du1) * Delta;
  info.pwvpx[1] = -du2 * Delta;
  info.pwvpx[2] = du1 * Delta;
  return info;
}

static int getVtxIndexInTriangle(int * tri, int vtx) 
{
  for(int i = 0; i < 3; i++)
    if (tri[i] == vtx) 
      return i;
  assert(0); // this should never be reached
  return 0;
}

double ClothBW::ComputeBendingStiffness(const BendInfo & info, const double * triangleUVs)
{
  // find the bend edge vtx
  int v[2] = {info.v[1], info.v[2]};
  int triA = info.tri[0], triB = info.tri[1];
  int g1 = triangleGroups[triA], g2 = triangleGroups[triB];  // material group index for triangle A, B

  Vec2d vtxUVA[2], vtxUVB[2]; // vtxUVA[2]: the UV coordinate of the two edge vtx (v1/v2)
  vtxUVA[0] = triangleUVs + 6*triA + 2*getVtxIndexInTriangle(&triangles[3*triA], v[0]);
  vtxUVA[1] = triangleUVs + 6*triA + 2*getVtxIndexInTriangle(&triangles[3*triA], v[1]);
  vtxUVB[0] = triangleUVs + 6*triB + 2*getVtxIndexInTriangle(&triangles[3*triB], v[0]);
  vtxUVB[1] = triangleUVs + 6*triB + 2*getVtxIndexInTriangle(&triangles[3*triB], v[1]);

  Vec2d duvA = vtxUVA[0] - vtxUVA[1]; // uv coordinate difference from edge vtx v1 to edge vtx v2, on triangle A
  Vec2d duvB = vtxUVB[0] - vtxUVB[1]; // uv coordinate difference from edge vtx v1 to edge vtx v2, on triangle B
  // on one triangle: k = (ku du^2 + kv dv^2) / (du^2 + dv^2)
  double stiffnessA = (materialGroups[g1].bendStiffnessU * duvA[0] * duvA[0] +  materialGroups[g1].bendStiffnessV * duvA[1] * duvA[1])
      / (duvA[0] * duvA[0] +  duvA[1] * duvA[1]);
  double stiffnessB = (materialGroups[g2].bendStiffnessU * duvB[0] * duvB[0] +  materialGroups[g2].bendStiffnessV * duvB[1] * duvB[1])
      / (duvB[0] * duvB[0] +  duvB[1] * duvB[1]);
  return (stiffnessA + stiffnessB) / 2.; // average bend stiffness from the two triangles
}

void ClothBW::GenerateBW(int numVertices_, const double * restPositions_, const double * masses_,
    int numTriangles_, const int * triangles_, const double * triangleUVs_, const int * triangleGroups_,
    int numMaterialGroups_, const MaterialGroup * materialGroups_, int addGravity_)
{
  bu = 1.0;
  bv = 1.0;
  // default to computing everything:
  cond[0] = true; // stretch&shear force
  cond[1] = true; // bend force
  cond[2] = true; // stretch&shear bend stiffness
  cond[3] = true; // bend stiffness matrix

  // default to taking rest angles into account
  useRestAnglesForBendingForces = 1;

  numVertices = numVertices_;
  restPositions.resize(numVertices);
  for(int i = 0; i < numVertices; i++)
    restPositions[i] = Vec3d(&restPositions_[3*i]);
  numTriangles = numTriangles_;
  triangles.resize(3 * numTriangles);
  memcpy(triangles.data(), triangles_, sizeof(int) * 3 * numTriangles);

  triangleGroups.resize(numTriangles);
  memcpy(triangleGroups.data(), triangleGroups_, sizeof(int) * numTriangles);

  masses.resize(numVertices);
  memcpy(masses.data(), masses_, sizeof(double) * numVertices);

  numMaterialGroups = numMaterialGroups_;
  materialGroups.resize(numMaterialGroups);
  memcpy(materialGroups.data(), materialGroups_, sizeof(MaterialGroup) * numMaterialGroups);

  // index all edges that can bend

  /* New ordering below
           B         D
           1---------3
          / \       /
         /   \  5  /
        /  4  \   /
       /       \ /
      0---------2
      A         C
   */
  // SORTED ensures that edges are paired in ascending order
  #define SORTED(i,j) ( (i) <= (j) ? make_pair((i),(j)) : make_pair((j),(i)) )
  typedef pair<int,int> Edge; // store sorted edge vtx index
  typedef pair<int, int> TriOpVtx; // store the triangle data (triangle index, the non opposite vtx on the triangle) for an edge
  map<Edge, TriOpVtx> edgeInfo; // map the edge to the <triangle index, opposite vtx> data
  map<Edge, TriOpVtx >::iterator iter;

  alphas.resize(numTriangles);
  wuvInfos.resize(numTriangles);
  for( int tri = 0 ; tri < numTriangles; tri++ )
  {
    const int * trivtx = &triangles[tri*3];
    assert(trivtx[0] >= 0 &&& trivtx[1] >= 0 && trivtx[2] >= 0);
    assert(trivtx[0] < numVertices && trivtx[1] < numVertices && trivtx[2] < numVertices);

    // compute the triangle surface area in UV coordinates (called "alpha" in [Baraff and Witkin 1998])
    double surfaceArea = GetTriangleSurfaceArea(restPositions[trivtx[0]], restPositions[trivtx[1]], restPositions[trivtx[2]]);

    // In [Baraff and Witkin 1998], alpha is set to be the surface area in UV coordinates
    //alphas[tri] = surfaceArea;
    // In David Pritchard's report (2012), it is reported that setting alpha = area^(3/4) is better if one wishes for cloth to have resolution-invariance
    alphas[tri] = pow(surfaceArea, 3.0 / 4.0);

    // compute wuv info
    wuvInfos[tri] = ComputeWuvInfo(&triangleUVs_[6*tri]);

    for( int j = 0 ; j < 3; j++)
    {
      Edge edge = SORTED(trivtx[j], trivtx[(j+1)%3]);
      int oppositeVtx = trivtx[(j+2)%3];

      if ( (iter = edgeInfo.find( edge )) == edgeInfo.end() ) // edge is not in edgeInfo, we insert it
        edgeInfo.insert( make_pair( edge, TriOpVtx(tri, oppositeVtx)));
      else // find one edge. We can form a bend edge with the new one and the old one
      {
        BendInfo info;
        TriOpVtx & previousTri = iter->second;
        // If we already visited this bend edge before and created a BendInfo,
        // there's an edge in the mesh that has three neighboring triangles.
        if (previousTri.first < 0) 
          throw 1;

        info.v[0] = previousTri.second; // index of the first triangle's vtx that is opposite to the edge
        info.v[1] = edge.first;         // index of the first edge vtx
        info.v[2] = edge.second;        // index of the second edge vtx
        info.v[3] = oppositeVtx;        // index of the second triangle's vtx that is opposite to the edge
        info.tri[0] = previousTri.first; // index of the first triangle
        info.tri[1] = tri;               // index of the second triangle
        quadComponentIndices.push_back(info);
        previousTri.first = -1; // mark this edge as visited
      }
    }
  }

  numQuads = quadComponentIndices.size();
  restAngles.resize(numQuads);
  bendStiffnesses.resize(numQuads);

  iter = edgeInfo.begin();

  // build the locations of non-zero entries in the stiffness matrix of size n x n
  SparseMatrixOutline skeletonOutline(numVertices);

  for (int i=0; i<numVertices; i++) // protection for isolated vertices
    skeletonOutline.AddEntry(i, i);

  for(int i=0; i<numTriangles; i++)
  {
    int * trivtx = &triangles[3*i];
    for(int j = 0; j < 3; j++)
      for(int k = 0; k < 3; k++)
        skeletonOutline.AddEntry(trivtx[j], trivtx[k]);
  }

  for(int i = 0 ; i < numQuads; i++)
  {
    int vertexA = quadComponentIndices[i].v[0];
    int vertexB = quadComponentIndices[i].v[1];
    int vertexC = quadComponentIndices[i].v[2];
    int vertexD = quadComponentIndices[i].v[3];

    skeletonOutline.AddEntry(vertexA, vertexD);
    skeletonOutline.AddEntry(vertexD, vertexA);

    // --- compute bend stiffness ---
    bendStiffnesses[i] = ComputeBendingStiffness(quadComponentIndices[i], triangleUVs_);

    // --- compute Resting Angles ---
    Vec3d restPositionA = restPositions[vertexA];
    Vec3d restPositionB = restPositions[vertexB];
    Vec3d restPositionC = restPositions[vertexC];
    Vec3d restPositionD = restPositions[vertexD];

    Vec3d vecBA = restPositionB - restPositionA;
    Vec3d vecCA = restPositionC - restPositionA;
    Vec3d vecBD = restPositionB - restPositionD;
    Vec3d vecCD = restPositionC - restPositionD;
    Vec3d vecBC = restPositionB - restPositionC;

    // new math ordering below
    Vec3d NA = cross(vecCA, vecBA); // normal for the first triangle
    Vec3d NB = cross(vecBD, vecCD);	// normal for the second triangle
    Vec3d E = vecBC; // edge vector

    // normalized normals and edge	
    Vec3d NAn = norm(NA);
    Vec3d NBn = norm(NB);
    Vec3d En = norm(E);

    double cosTheta = dot(NAn, NBn);
    double sinTheta;

    Vec3d NANB = cross(NAn, NBn); // storing the cross product of NAn and NBn

    sinTheta = dot(NANB, En);

    restAngles[i] = atan2(sinTheta, cosTheta);
  }

  // build inverse indices for stiffness matrix access
  SparseMatrix skeleton(&skeletonOutline);

  inverseIndicesStretchAndShear.resize(9 * numTriangles);
  for(int i=0; i<numTriangles; i++)
  {
    int * trivtx = &triangles[3*i];
    for(int j = 0, l = 0; j < 3; j++)
      for(int k = 0; k < 3; k++)
        inverseIndicesStretchAndShear[9*i+(l++)] = skeleton.GetInverseIndex(trivtx[j], trivtx[k]);
  }

  inverseIndicesQuad.resize(16 * numQuads);
  for(int i = 0 ; i < numQuads; i++)
  {
    const int * quadvtx = quadComponentIndices[i].v;
    for(int j = 0, l = 0; j < 4; j++)
      for(int k = 0; k < 4; k++)
        inverseIndicesQuad[16*i + (l++)] = skeleton.GetInverseIndex(quadvtx[j], quadvtx[k]);
  }
}

ClothBW::~ClothBW()
{
}

void ClothBW::GenerateMassMatrix(SparseMatrix **M, int expanded) const
{
  SparseMatrixOutline outline(expanded * numVertices);
  for(int i=0; i<numVertices; i++)
    for(int j=0; j<expanded; j++)
      outline.AddEntry(expanded*i+j, expanded*i+j, masses[i]); 
  *M = new SparseMatrix(&outline);	
}

double ClothBW::GetTriangleSurfaceArea(const Vec3d & p0, const Vec3d & p1, const Vec3d & p2)
{
  Vec3d s0 = p1 - p0;
  Vec3d s1 = p2 - p0;
  return 0.5 * len(cross(s0, s1));
}

void ClothBW::SetComputationMode(const bool conditions[4])
{
  memcpy(cond, conditions, sizeof(bool) * 4);
}

void ClothBW::ComputeGravity(double * f)
{
  for(int i=0; i<numVertices; i++)
    f[3*i+1] += g * masses[i];
}

void ClothBW::AddStretchAndShear(const double * u, int startTriangle, int endTriangle, double * energy, double * f, SparseMatrix * K)
{
  for(int tri = startTriangle; tri < endTriangle; tri++)
  {
    int group = triangleGroups[tri]; // group index

    int p0 = triangles[3*tri+0]; // triangle vtx indices
    int p1 = triangles[3*tri+1];
    int p2 = triangles[3*tri+2];

    Vec3d pos0 = restPositions[p0] + Vec3d(&u[3*p0]);
    Vec3d pos1 = restPositions[p1] + Vec3d(&u[3*p1]);
    Vec3d pos2 = restPositions[p2] + Vec3d(&u[3*p2]);

    Vec3d vec01 = pos1 - pos0;
    Vec3d vec02 = pos2 - pos0;

    const double * dwudx = wuvInfos[tri].pwupx;
    const double * dwvdx = wuvInfos[tri].pwvpx;
    Vec3d wu = dwudx[1] * vec01 + dwudx[2] * vec02; // = (dv2 * vec01 - dv1 * vec02) * Delta;
    Vec3d wv = dwvdx[1] * vec01 + dwvdx[2] * vec02; // = ( -du2 * vec01 + du1 * vec02) * Delta;
    double length_wu = len(wu);
    double length_wv = len(wv);
    Vec3d wun = wu / length_wu; //wu normalized
    Vec3d wvn = wv / length_wv; // wv normalized

    double alpha = alphas[tri];

    // --- compute stretch and shear energy ---
    // stretch energy Es = 0.5 Cs^T Cs
    // Cs = [Csu Csv]^T = alpha [ |wu|-bu |wv| -bv ]^T
    // shear energy E_h = 0.5 Ch^T Ch
    // Ch = alpha wu^T wv
    double Cstru = alpha * (length_wu - bu); // stretch energy in u: Csu
    double Cstrv = alpha * (length_wv - bv); // stretch energy in v: Csv
    double Cshear = alpha * dot(wu, wv);     // shear energy       : Ch

    if (energy)
    {
      *energy += 0.5 * materialGroups[group].tensileStiffness * (Cstru*Cstru + Cstrv * Cstrv);
      *energy += 0.5 * materialGroups[group].shearStiffness * Cshear * Cshear;
    }

    // --- compute the derivatives of Csu, Csv, Ch ---
    // For stretch:
    // 3x1 vector: pCsu/px = alpha pwu/px \hat{wu} ,  pCsv/px = alpha pwv/px \hat{wv} , \hat{W} means the normalized W

    // Define pwu/pxi = wuxi * I, pwv/pxi = wvxi * I, where wuxi, wvxi are both scalars
    // So,
    // pCsu/px0 = alpha * wux0 * \hat{Wu}
    // pCsu/px1 = alpha * wux1 * \hat{Wu}
    // pCsu/px2 = alpha * wux2 * \hat{Wu}
    // pCsv/px0 = alpha * wvx0 * \hat{Wv}
    // pCsv/px1 = alpha * wvx1 * \hat{Wv}
    // pCsv/px2 = alpha * wvx2 * \hat{Wv}

    // For shear:
    // 3x1 vector: pCh/pxi = alpha ( pwv/pxi wu + pwu/pxi wv) = alpha ( wvxi wu + wuxi wv )
    Vec3d pCsupx[3]; // stretch energy derivative on u: pCsu/px0, pCsu/px1, pCsupx2, <x0, x1, x2> vtx pos./disp. on one triangle
    Vec3d pCsvpx[3]; // stretch energy derivative on u: pCsv/px0, pCsv/px1, pCsvpx2
    Vec3d pChpx[3];  // shear energy derivative: pChpx0, pChpx1, pChpx2
    for(int vtx = 0; vtx < 3; vtx++)
    {
      pCsupx[vtx] = (alpha * dwudx[vtx]) * wun;
      pCsvpx[vtx] = (alpha * dwvdx[vtx]) * wvn;
      pChpx[vtx] = alpha * (dwudx[vtx] * wv + dwvdx[vtx] * wu);
    }

    // --- compute force ---
    // stretch energy Es = 0.5 [Csu Csv] [Csu Csv]^T
    // 3x1 vector: pEs/px = Csu pCsu/px + Csv pCsv/px
    // shear energy Eh = 0.5 Ch^2
    // 3x1 vector pEh/pxi = Ch pCh/pxi
    // force = pEspx + pEhpx
    if (f)
    {
      for(int vtx = 0; vtx < 3; vtx++)
      {
        int vtxid = triangles[3*tri+vtx];
        Vec3d force(0.0);
        force += materialGroups[group].tensileStiffness * (Cstru * pCsupx[vtx] + Cstrv * pCsvpx[vtx]);
        force += materialGroups[group].shearStiffness * Cshear * pChpx[vtx];
        force.addToArray(f + 3*vtxid);
      }
    }

    // --- compute stiffness matrix ---
    if (K)
    {
      for(int i = 0 ; i < 3 ; i++) // force vtx i
        for(int j = 0; j < 3; j++) // disp. vtx j
        {
          // second derivative of stretch Csu, Csv
          // p(pCsu/pxi)/pxj = (alpha / |Wu|) wuxi wuxj (I - \hat{Wu} \hat{Wu}^T)
          Mat3d p2Csupxij = ((alpha / length_wu) * dwudx[i] * dwudx[j]) * (Mat3d(1.) - tensorProduct(wun, wun));
          Mat3d p2Csvpxij = ((alpha / length_wv) * dwvdx[i] * dwvdx[j]) * (Mat3d(1.) - tensorProduct(wvn, wvn));
          // p2Es/pxipxj = p(pEs/pxi)/pxj = pCsu/pxi pCsu/pxj^T + Csu p(pCsu/pxi)/pxj + pCsv/pxi pCsv/pxj^T + Csv p(pCsv/pxi)/pxj
          Mat3d Kstrij = tensorProduct(pCsupx[i], pCsupx[j]) + Cstru * p2Csupxij
              + tensorProduct(pCsvpx[i], pCsvpx[j]) + Cstrv * p2Csvpxij;
          Kstrij *= materialGroups[group].tensileStiffness;

          // second derivative of shear Ch:
          // p2Eh/pxipxj = pCh/pxi pCh/pxj^T + Ch p2Ch/pxipxj
          // p2Ch/pxipxj = alpha ( wvxi wuxj + wuxi wvxj )
          Mat3d Kshij = tensorProduct(pChpx[i], pChpx[j]) +
              (Cshear * alpha * (dwvdx[i] * dwudx[j] + dwudx[i] * dwvdx[j])) * Mat3d(1.);
          Kshij *= materialGroups[group].shearStiffness;

          for(int m = 0 ; m < 3; m++)     // loop inside force vtx i component
            for(int n = 0 ; n < 3; n++)   // loop inside disp vtx j component
            {
              // Adding it to K topology
              int entry_row = 3 * triangles[3*tri+i] + m;
              int entry_col = 3 * inverseIndicesStretchAndShear[9*tri+i*3+j] + n;
              K->AddEntry(entry_row, entry_col, Kstrij[m][n] + Kshij[m][n]);
            }
        } // end i,j
    } // end if (K)
  }
}

void ClothBW::AddBend(const double * u, int startQuad, int endQuad, double * energy, double * f, SparseMatrix * K)
{
  for(int qua = startQuad ; qua < endQuad; qua++)
  {
    int p0 = quadComponentIndices[qua].v[0];
    int p1 = quadComponentIndices[qua].v[1];
    int p2 = quadComponentIndices[qua].v[2];
    int p3 = quadComponentIndices[qua].v[3];

    double kBend = bendStiffnesses[qua];

    // vtx Indices defined below
    //
    //       1---------3
    //      / \       /
    //     /   \  B  /
    //    /  A  \   /
    //   /       \ /
    //  0---------2

    Vec3d pos0 = restPositions[p0] + Vec3d(&u[3*p0]);
    Vec3d pos1 = restPositions[p1] + Vec3d(&u[3*p1]);
    Vec3d pos2 = restPositions[p2] + Vec3d(&u[3*p2]);
    Vec3d pos3 = restPositions[p3] + Vec3d(&u[3*p3]);

    Vec3d v01 = pos1 - pos0;
    Vec3d v02 = pos2 - pos0;
    Vec3d v31 = pos1 - pos3;
    Vec3d v32 = pos2 - pos3;
    Vec3d v21 = pos1 - pos2;

    Vec3d v12 = -v21;
    Vec3d v20 = -v02;
    Vec3d v13 = -v31;

    // normals for triangles A and B
    // nA = cross(x2-x0, x1-x0)
    // nB = cross(x1-x3, x2-x3) // force vtx m
    // e = x1 - x2
    Vec3d nA = cross(v02, v01);
    Vec3d nB = cross(v31, v32);
    Vec3d edge = v21;

    double nAL = len(nA);
    double nBL = len(nB);
    double eL = len(edge);
    double invnAL = 1.0 / nAL;
    double invnBL = 1.0 / nBL;
    assert(eL != 0.0);
    double inveL = 1.0 / eL;

    // now normalize the... normals
    Vec3d nAN = nA * invnAL;
    Vec3d nBN = nB * invnBL;
    Vec3d eN = edge * inveL;

    // --- compute bend energy Cb ---
    // bend energy E_b = 0.5 Cb^T Cb
    // Cb = theta - theta0
    // pEb/px = Cb pCb/pxi

    // calculate sin theta and cos theta
    double cosTheta = dot(nAN, nBN);
    double sinTheta = dot(cross(nAN, nBN), eN);
    double cbend = atan2(sinTheta, cosTheta); // store theta into cbend

    // account for the rest angles?
    if (useRestAnglesForBendingForces)
      cbend -= restAngles[qua];
    // now cbend stores Cb = theta - theta0

    if (energy)
      *energy += 0.5 * kBend * cbend * cbend;

    Vec3d qA[4] = { v12, v20, v01, Vec3d(0.0) };
    Vec3d qB[4] = { Vec3d(0.0), v32, v13, v21 };
    const double qe[4] = {0, 1, -1, 0};

    // --- compute derivatives of triangle A normal \hat{nA}, triangle B normal \hat{nB} and edge direction \hat{e} ---
    Mat3d phnApx[4], phnBpx[4], phepx[4];
    for(int vtx = 0; vtx < 4; vtx++)
    {
      // Here we compute the derivatives based on approximation that the normal and edge vectors have constant magnitude
      // phnA/pxms = 1.0 / |nA| * S_s(qAm)
      phnApx[vtx] = invnAL * skewSymmetricMatrix(qA[vtx]); // now phnA/pxms = phnApx[m][s]
      phnBpx[vtx] = invnBL * skewSymmetricMatrix(qB[vtx]);
      phepx[vtx] = inveL * qe[vtx] * Mat3d(1.0);
    }

    // --- compute derivatives of sin(theta), cos(theta) and theta ---
    Vec3d pSinpx[4], pCospx[4], pThetapx[4];
    for(int m = 0; m < 4; m++) // over 4 vtx
    {
      pCospx[m] = phnApx[m] * nBN + phnBpx[m] * nAN;
      for(int s = 0; s < 3; s++)
      {
        pSinpx[m][s] = dot(cross(phnApx[m][s], nBN) + cross(nAN, phnBpx[m][s]), eN) +
                       dot(cross(nAN, nBN), phepx[m][s]);
        pThetapx[m][s] = cosTheta * pSinpx[m][s] - sinTheta * pCospx[m][s];
      }

      // --- compute force ---
      if (f)
      {
        // Cb = theta - theta0
        // 3x1 vector: force = pEb/pxm = Cb pCb/pxm
        Vec3d force = kBend * cbend * pThetapx[m];
        int vtxidx = quadComponentIndices[qua].v[m];
        force.addToArray(f + 3*vtxidx);
      }
    }

    // --- compute stiffness matrix ---
    if (K)
    {
      // qA is used to compute pqAm/pxnt
      const double qA[4][4] = { {0,-1,1,0}, {1,0,-1,0}, {-1,1,0,0}, {0,0,0,0} };
      const double qB[4][4] = { {0,0,0,0}, {0,0,1,-1}, {0,-1,0,1}, {0,1,-1,0} };
      for(int m = 0; m < 4; m++)  // force vtx m
        for(int n = 0; n < 4; n++) // disp. vtx n
        {
          // compute K_mn,  3x3 matrix in K at rows 3*m -> 3*m+2, columns: 3*n -> 3*n+2

          // compute second derivatives of \hat{nA} and \hat{nB} given vtx index m and n fixed
          Mat3d p2hnApt[3], p2hnBpt[3];
          for(int t = 0; t < 3; t++) // loop over the component of disp. vtx n
          {
            Vec3d pqAmpxnt(0.0), pqBmpxnt(0.0);
            pqAmpxnt[t] = qA[m][n];              //pqAmpxnt = partial qAm / partial xnt
            pqBmpxnt[t] = qB[m][n];
            p2hnApt[t] = invnAL * skewSymmetricMatrix(pqAmpxnt); // p2hnApt[t][s] stores p2 hnA / pxms pxnt
            p2hnBpt[t] = invnBL * skewSymmetricMatrix(pqBmpxnt);
          }

          // compute second derivatives of cos(theta), sin(theta) and theta
          Mat3d p2Cospst, p2Sinpst, p2Thetapst;
          for(int s = 0; s < 3; s++)
            for(int t = 0; t < 3; t++)
            {
              p2Cospst[s][t] = dot(p2hnApt[t][s], nBN) + dot(phnBpx[n][t], phnApx[m][s])
                  + dot(phnApx[n][t], phnBpx[m][s]) + dot(nAN, p2hnBpt[t][s]);

              p2Sinpst[s][t] = dot(cross(p2hnApt[t][s], nBN) + cross(phnApx[m][s], phnBpx[n][t]) +
                  cross(phnApx[n][t], phnBpx[m][s]) + cross(nAN, p2hnBpt[t][s]), eN) +
                  dot(cross(phnApx[m][s], nBN) + cross(nAN, phnBpx[m][s]), phepx[n][t]) +
                  dot(cross(phnApx[n][t], nBN) + cross(nAN, phnBpx[n][t]), phepx[m][s]);

              p2Thetapst[s][t] = cosTheta * p2Sinpst[s][t] - sinTheta * p2Cospst[s][t] +
                  (sinTheta * sinTheta - cosTheta * cosTheta) * (pSinpx[m][s] * pCospx[n][t] + pCospx[m][s] * pSinpx[n][t])
                  + 2 * sinTheta * cosTheta * (pCospx[m][s] * pCospx[n][t] - pSinpx[m][s] * pSinpx[n][t]);
            }

          // Cb = theta - theta0
          // 3x1 vector: pCb/px
          // 3x1 vector: pEb/pxm = Cb pCb/pxm
          // p2Eb/pxmxn = pCb/pxm pCb/pxn^T + Cb p2Cb/pxmxn
          Mat3d p2Eb = kBend * ( tensorProduct(pThetapx[m], pThetapx[n]) + cbend * p2Thetapst );

          for(int s = 0; s < 3; s++)
            for(int t = 0; t < 3; t++)
            {
              int entry_row = 3 * quadComponentIndices[qua].v[m] + s;
              int entry_col = 3 * inverseIndicesQuad[qua*16+ m*4+ n] + t;
              K->AddEntry(entry_row, entry_col, p2Eb[s][t]);
            }
        } // end n 0 -> 4, // end m 0 -> 4
    } // end if K
  } // end quad
}

//void ClothBW::ComputeDampingForce(const double *u, double *uvel, double *f, bool addForce)
//{
//  // unimplemented
//}

//void ClothBW::AddDampingForce(double *uvel, double *f, int startTriangle, int endTriangle)
//{
//  // unimplemented
//}

void ClothBW::GenerateStiffnessMatrixTopology(SparseMatrix **K) const
{
  SparseMatrixOutline KOutline(3*numVertices);
  for(int vtx=0; vtx < numVertices; vtx++) // for all vtx. in case there's an isolated vtx
    KOutline.AddBlock3x3Entry(vtx, vtx);

  // for per triangle
  for(int i=0; i<numTriangles; i++)
  {
    const int * trivtx = &triangles[3*i];
    for(int j=0; j<3; j++)
      for(int k=0; k<3; k++)
        KOutline.AddBlock3x3Entry(trivtx[j], trivtx[k]);
  }

  // for per quad, in addition to per triangle topology
  for( int i = 0 ; i < numQuads; i++ )
  {
    int vertexA = quadComponentIndices[i].v[0];
    int vertexD = quadComponentIndices[i].v[3];

    KOutline.AddBlock3x3Entry(vertexA, vertexD);
    KOutline.AddBlock3x3Entry(vertexD, vertexA);
  }

  *K = new SparseMatrix(&KOutline);
}

double ClothBW::ComputeEnergy(const double * u)
{
  double energy = 0.0;
  if (cond[0])
    AddStretchAndShear(u, 0, numTriangles, &energy, NULL, NULL);
  if (cond[1])
    AddBend(u, 0, numQuads, &energy, NULL, NULL);
  return energy;
}

void ClothBW::ComputeForce(const double * u, double * f, bool addForce)
{
  if (!addForce) 
    memset(f, 0, sizeof(double) * 3 * numVertices);

  if (cond[0])
    AddStretchAndShear(u, 0, numTriangles, NULL, f, NULL);

  if (cond[1])
    AddBend(u, 0, numQuads, NULL, f , NULL);

  if (addGravity) 
    ComputeGravity(f);
}

void ClothBW::ComputeStiffnessMatrix(const double * u, SparseMatrix * K, bool addMatrix)
{
  if (!addMatrix) 
    K->ResetToZero();

  if (cond[2])
    AddStretchAndShear(u, 0, numTriangles, NULL, NULL, K);

  if (cond[3])
    AddBend(u, 0, numQuads, NULL, NULL, K);
}

void ClothBW::ComputeForceAndMatrix(const double * u, double * f, SparseMatrix * K, bool addQuantity)
{
  if (addQuantity == false)
  {
    memset(f, 0, sizeof(double) * 3 * numVertices);
    K->ResetToZero();
  }

  if (cond[0] || cond[2])
    AddStretchAndShear(u, 0, numTriangles, NULL, (cond[0] ? f : NULL), (cond[2] ? K : NULL));

  if (cond[1] || cond[3])
    AddBend(u, 0, numQuads, NULL, (cond[1] ? f : NULL), (cond[3] ? K : NULL));

  if (addGravity) 
    ComputeGravity(f);
}

