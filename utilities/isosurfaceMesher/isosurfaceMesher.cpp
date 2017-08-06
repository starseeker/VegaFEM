/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 3.0                               *
 *                                                                       *
 * "isosurfaceMesher" utility , Copyright (C) 2016 USC                   *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Danyong Zhao, Yijing Li, Jernej Barbic                  *
 * http://www.jernejbarbic.com/code                                      *
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
  Driver for the IsosurfaceMesher. See the mesher library for more details.

  Computes a quality triangular mesh of a distance field isosurface (level set),
  using the following paper:

  Steve Oudot, Laurent Rineau, Mariette Yvinec:
  Meshing Volumes Bounded by Smooth Surfaces
  Proceedings of the 14th International Meshing Roundtable 2005, pp 203-219
 */

#include <iostream>

#include "tetMesh.h"
#include "delaunayMesher.h"
#include "distanceFieldCreator.h"
#include "isosurfaceMesher.h"
#include "performanceCounter.h"
#include "getopts.h"
#include "objMeshOrientable.h"
#include "rVector.h"
#include "matrixIO.h"
using namespace std;

int main(int argc, char **argv)
{
  int numFixedArg = 4;
  if (argc < numFixedArg)
  {
    cout << "Computes a quality triangular mesh of a distance field isosurface (level set)" << endl;
    cout << argv[0] << ": <distance field> <triangleSize> <out obj mesh> (-i<isoValue> -a<angle limit(in degrees, 0<a<=30)> "
        "-e<epsilon> -q<query type: 0:double, 1:exact, 2:filtered>, -m<input is instead a surface mesh>, "
        "-c <check Delaunay> -e <enforce to output a 2-manifold mesh> -n <input a narrowband distance field>)" << endl;
    return 0;
  }
  char * fieldFilename = argv[1];
  char * surfaceMeshFilename = argv[1];
  double triangleSize = atof(argv[2]);
  char * objMeshFilename = argv[3];
  double isoValue = 0;
  bool inputSurfaceMesh = false;
  double angleLimit = 30;
  double epsilon = 1e-6;
  int queryType = 0;
  bool checkDelaunay = false;
  char angleLimitString[4096] = "30.0";
  char isoValueString[4096] = "0.0";
  char epsilonString[4096] = "1e-6";
  bool narrowBand = false;
  bool enforceManifold = true;
  opt_t opttable[] =
  {
    { "a", OPTSTR, angleLimitString },
    { "i", OPTSTR, isoValueString },
    { "m", OPTBOOL, &inputSurfaceMesh },
    { "e", OPTSTR, epsilonString },
    { "q", OPTINT, &queryType },
    { "c", OPTBOOL, &checkDelaunay },
    { "n", OPTBOOL, &narrowBand },
    { "e", OPTBOOL, &enforceManifold },
    { NULL, 0, NULL }
  };

  argv += 3;
  argc -= 3;
  int optup = getopts(argc, argv, opttable);
  if (optup != argc)
  {
    printf("Error parsing options. Error at option %d: %s.\n", optup, argv[optup + numFixedArg - 1]);
    return 1;
  }

  angleLimit = strtod(angleLimitString, NULL);
  isoValue = strtod(isoValueString, NULL);
  epsilon = strtod(epsilonString, NULL);
  PerformanceCounter pc;

  DistanceFieldBase * field;
  if (narrowBand)
    field = new DistanceFieldNarrowBand;
  else
    field = new DistanceField;
  ObjMesh * surfaceMesh = NULL;

  if (angleLimit > 30)
    angleLimit = 30;
  else if (angleLimit < 0)
    angleLimit = 1;
  angleLimit = angleLimit * M_PI / 180.;
  double radiusLimit = triangleSize;
  if (radiusLimit < 0)
    radiusLimit = 2 * field->diagonal();

  if (inputSurfaceMesh == false)
  {
    int ret = field->load(fieldFilename);
    if (ret != 0)
    {
      cout << "Failed to load distance field from " << fieldFilename << endl;
      return 1;
    }
    cout << "Time to load distance field: " << pc.GetElapsedTime() << endl;

    //  double gridX,gridY,gridZ;
    //  field.getGridSpacing(&gridX, &gridY, &gridZ);
    //  double gridSize = max(gridX, max(gridY, gridZ))pc.GetElapsedTime();;
    //  cout << "gridSize: " << gridSize << endl;
    //  //double triangleSize = 1.;
    //
    //  double absTriangleLimit = triangleSize * gridSize;
    //  cout << "triangle size (world unit): " << absTriangleLimit << endl;

    Vec3d bmin = field->bmin();
    Vec3d bmax = field->bmax();

    // convert to absolute value
    cout << "radiusLimit: " << radiusLimit << endl;
    cout << "isoValue: " << isoValue << endl;
  }
  else
  {
    surfaceMesh = new ObjMesh(surfaceMeshFilename);
    Vec3d bmin, bmax;
    surfaceMesh->getBoundingBox(1.0, &bmin, &bmax);
    Vec3d grid = (bmax - bmin) / 1024;
  }

  bool outputManfold = false;
  pc.StartCounter();
  ObjMesh * mesh = NULL;

  Query::Type type = Query::DOUBLE;
  if (queryType == 1)
    type = Query::RATIONAL;
  else if (queryType == 2)
    type = Query::FILTERED;

  if (inputSurfaceMesh)
  {
    IsosurfaceMesher mesher(surfaceMesh);
    mesher.checkDelaunayAfterComputation(checkDelaunay);

    mesher.compute(isoValue, angleLimit, radiusLimit, 20, epsilon, type, 10000);
    mesh = mesher.getMesh(outputManfold);
  }
  else
  {
    IsosurfaceMesher mesher(field);
    mesher.checkDelaunayAfterComputation(checkDelaunay);
    mesher.compute(isoValue, angleLimit, radiusLimit, 20, epsilon, type);
    mesh = mesher.getMesh(outputManfold);
  }
  pc.StopCounter();
  if (mesh == NULL)
  {
    cout << "Fail to generate isosurface mesh" << endl;
    return 1;
  }
  cout << "Time to generate mesh: " << pc.GetElapsedTime() << endl;

  if (enforceManifold)
    IsosurfaceMesher::enforceManifoldnessAndOrientNormals(mesh);
  else
  {
    try
    {
      ObjMeshOrientable orientedMesh(mesh);
      cout << "Mesh is manifold" << endl;
      if (orientedMesh.hasBoundary())
      {
        cout << "Mesh has boundary" << endl;
      }
      else
      {
        cout << "Mesh has no boundary" << endl;
      }
    }
    catch (int)
    {
      cout << "Fail to orient the mesh" << endl;
    }
  }
  mesh->save(objMeshFilename);
  cout << "Save to " << objMeshFilename << endl;
  //delete mesh;
  delete surfaceMesh;
  delete mesh;
  delete field;

  return 0;
}

