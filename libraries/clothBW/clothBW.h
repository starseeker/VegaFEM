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
 This class implements the cloth model presented in:

 Baraff and Witkin: Large Steps in Cloth Simulation, SIGGRAPH 1998. 
 
 It can compute the elastic energy, the internal elastic forces, and the tangent stiffness matrix,
 for stretch/shear and bend terms.
 Damping is not implemented, but you can use the damping in the Vega integrator class.

 Our implementation follows the original Baraff and Witkin paper, with the insights
 presented in the following report:

 David Pritchard: Implementing Baraff & Witkin's Cloth Simulation, May 2003, with minor updates April 2012
*/

#ifndef _CLOTHBW_H_
#define _CLOTHBW_H_

#include "sparseMatrix.h"
#include "vec3d.h"

class ClothBW
{
public:

  // cloth material parameters
  struct MaterialGroup
  {
    double tensileStiffness;
    double shearStiffness;
    double bendStiffnessU;
    double bendStiffnessV;
    double damping;
  };

  // creates the cloth elastic model, from a given triangle mesh
  // "masses" is an array of length "numVertices"
  // "restPositions" is an array of length 3x"numVertices"
  // "triangles" is an integer array of length 3x"numTriangles" (giving integer indices of the three particle forming a triangle)
  // "triangleGroups" is an integer array of length "numTriangles" (giving the integer index of the material group to which each triangle belongs)
  // "triangleUVs" is an double array of length 3x2x"numTriangles", indicating the uv for every vertex; this array can be user-provided, or the constructor can compute it automatically
  // all indices in this class are 0-indexed
  // constructor that does not require triangleUVs input (it computes UVs automatically; note: the UVs are continuous only within each triangle; the UV map is not global (which is fine, provided one does not want to simulate anisotropic effects) )
  ClothBW(int numVertices, const double * restPositions, const double * masses,
      int numTriangles, const int * triangles, const int * triangleGroups,
      int numMaterialGroups, const MaterialGroup * materialGroups, int addGravity=0);
  
  // constructor with triangleUVs input
  ClothBW(int numVertices, const double * restPositions, const double * masses,
      int numTriangles, const int * triangles, const double * triangleUVs, const int * triangleGroups,
      int numMaterialGroups, const MaterialGroup * materialGroups, int addGravity=0);
    
  virtual ~ClothBW();

  int GetNumVertices() const { return numVertices; }
  int GetNumTriangles() const { return numTriangles; }
  const int * GetTriangles() const { return triangles.data(); }
  const Vec3d * GetRestPositions() const { return restPositions.data();  }

  void SetRestUVStretchValues(double bu, double bv) { this->bu = bu; this->bv = bv; } // these are the b_u and b_v stretch values from Equation (10) in the BW paper (default values are 1.0)

  // computes the gravitational force (result goes into f)
  void SetGravity(bool addGravity, double g=9.81) { this->addGravity = addGravity; this->g = g; } // if addGravity is enabled, ComputeForces will add the gravity force to the elastic forces

  // this allows users to toggle on/off computation of the stretch/shear forces, bend forces,
  // stretch/shear stiffness matrix, and bend stiffness matrix
  // the function ComputeEnergy uses mode[0] and mode[1] flags
  // mode[0] = computeStretchAndShearForce
  // mode[1] = computeBendForce
  // mode[2] = computeStretchAndShearStiffnessMatrices
  // mode[3] = computeBendStiffnessMatrices
  void SetComputationMode(const bool mode[4]);

  // allows user to toggle on/off the use of rest angles in bend force/stiffness matrix
  // calculations; if set to 1, bend force/matrix will be computed in relation to the quad's rest
  // angle. if set to 0, bend force/matrix will be computed in relation to a flat angle of 0.0
  // default is 1.
  void UseRestAnglesForBendingForces(bool useRestAnglesForBend) { useRestAnglesForBendingForces = useRestAnglesForBend;}

  // creates the mass matrix (which is diagonal); each diagonal entry is expanded into a diagonal submatrix of size 'expanded' (typically, for 3D simulations, expanded should be 3)
  void GenerateMassMatrix(SparseMatrix ** M, int expanded=3) const;

  // === compute elastic energy, force and stiffness matrix ===

  // compute the elastic energy, under deformation u
  virtual double ComputeEnergy(const double * u); 

  // compute the internal elastic force, under deformation u
  // note: the force has the sign of the left side of the dynamic equation, Mu'' + Du' + f_int(u) = f_ext(t), i.e., f_int(u), that is, **opposite** to an external force f_ext(t) acting on the body 
  virtual void ComputeForce(const double * u, double * f, bool addForce=false); // if addForce is "true", f will be not be reset to zero prior to adding the forces

  // compute the tangent stiffness matrix of the elastic force
  // call once to establish the location of sparse entries of the stiffness matrix
  void GenerateStiffnessMatrixTopology(SparseMatrix ** K) const;

  virtual void ComputeStiffnessMatrix(const double * u, SparseMatrix * K, bool addMatrix=false);

  virtual void ComputeForceAndMatrix(const double * u, double * f, SparseMatrix * K, bool addQuantity = false);
    
  // compute the damping force
  // unimplemented
  // (use damping provided in the integrator class)
  // virtual void ComputeDampingForce(const double * u, double * uvel, double * f, bool addForce=false);

protected:

  struct WuvInfo // derivatives of wu, wv with regard to vtx displacements
  {
    double pwupx[3];
    double pwvpx[3];
  };

  struct BendInfo // store vtx and triangle indices for computing bend force
  {
    int v[4];
    int tri[2];
  };

  static double GetTriangleSurfaceArea(const Vec3d & p0, const Vec3d & p1, const Vec3d & p2);

  static WuvInfo ComputeWuvInfo(const double triangleUV[6]);

  double ComputeBendingStiffness(const BendInfo & bendInfo, const double * triangleUVs);

  void GenerateBW(int numVertices, const double * restPositions, const double * masses,
      int numTriangles, const int * triangles, const double * triangleUVs, const int * triangleGroups,
      int numMaterialGroups, const MaterialGroup * materialGroups, int addGravity=0);
    
  void AddStretchAndShear(const double * u, int startTriangle, int endTriangle, double * energy, double * f, SparseMatrix * K);

  void AddBend(const double * u, int startQuad, int endQuad, double * energy, double * f, SparseMatrix * K);

  // unimplemented
  // (use damping provided in the integrator class)
  //void AddDampingForce(double * uvel, double * f, int startTriangle, int endTriangle);

  void ComputeGravity(double * f);

  int numVertices;
  std::vector<double> masses;
  std::vector<Vec3d> restPositions;
  int numTriangles;
  std::vector<int> triangles;
//  std::vector<double> triangleUVs;

  std::vector<WuvInfo> wuvInfos;
  std::vector<int> inverseIndicesStretchAndShear;
  std::vector<int> triangleGroups;
  std::vector<double> alphas; // a scaling factor used in stretch and shear energy
    
  // Internal variables for bending computation:
  // "numQuads" is number of 'quads' made up of two adjacent triangles. Important because each quad contains a bendable edge (the edge shared by the two triangles).
  // "restAngles" is array of size numQuads containing the rest angles of edges that can bend (indexed by numQuads)
  // "inverseIndicesQuad" is numQuads x 16 integer array
  // "quadComponentIndices" is numQuads x 6 integer array
  int numQuads;
  std::vector<double> restAngles;       // rest angle for each quad
  std::vector<double> bendStiffnesses;  // bend stiffness for each quad
  std::vector<int> inverseIndicesQuad;

  std::vector<BendInfo> quadComponentIndices;
    
  int numMaterialGroups;

  std::vector<MaterialGroup> materialGroups;
    
  double bu; // stretch constraint in u (default = 1.0)
  double bv; // stretch constraint in v (default = 1.0)
    
  int addGravity;
  double g;
    
  bool cond[4];
  bool useRestAnglesForBendingForces;
};

#endif

