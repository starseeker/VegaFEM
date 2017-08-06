/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 3.0                               *
 *                                                                       *
 * "tet mesher" driver , Copyright (C) 2016 USC                          *
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
  Driver for constrained 3D Delaunay tet mesher with refinement.
  See the mesher library for more info.

  The input is a manifold and self-intersection-free triangle mesh,
  and the output is a quality tet mesh of the volume enclosed by the input mesh
  that conforms to the input mesh.
  This mesh is generated using constrained Delaunay tetrahedralization,
  and by inserting new vertices (``Steiner vertices'') into the mesh as needed
  to satisy the quality conditions.

  Note: The tet mesher is experimental. It works in most cases, but
  may generate sliver tets with degenerate inputs, such as all input points lying
  on the same sphere in 3D.
*/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fstream>
#include <iostream>
using namespace std;

#include "getopts.h"
#include "objMesh.h"
#include "tetMesher.h"

int main( int argc, char** argv )
{
  if ( argc < 3 ) 
  {
    std::cout << "Creates a tet mesh, by meshing the space occupied by the given input manifold surface mesh." << std::endl;
    std::cout << "Usage: " << argv[0] << " [input surface mesh file] [output tet mesh filename] [-a alpha] [-q refinement quality] [-s max Steiner vertices] [-d minimal dihedral angles the tetmesh can have] [-t maximal seconds the refine process takes]" << std::endl;
    std::cout << "Refinement quality is a scalar. It must be greater or equal to 1.0. The lower the number, the more the mesh will be refined. Default is 1.1." << std::endl;
    std::cout << "Max Steiner points controls the maximum number of inserted vertices. Default: unbounded." << std::endl;
    std::cout << "If a point is closer than alpha times the average edge length (in input mesh), skip adding this Steiner point." << std::endl;
    return 1;
  }

  char * inputSurfaceMeshFilename = argv[1];
  char * ouputTetMeshFilename = argv[2];
  char refinementQualityString[4096] = "1.1";
  char alphaString[4096] = "1.0";
  int maxSteinerVertices = -1;
  char minDihedralString[4096] = "0.0";
  char maxTimeSecondsString[4096] = "-1.0";

  opt_t opttable[] =
  {
    { "a", OPTSTR, alphaString },
    { "q", OPTSTR, refinementQualityString },
    { "s", OPTINT, &maxSteinerVertices },
    { "d", OPTSTR, minDihedralString },
    { "t", OPTSTR, maxTimeSecondsString },
    { NULL, 0, NULL }
  };

  argv += 2;
  argc -= 2;
  int optup = getopts(argc,argv,opttable);
  if (optup != argc)
  {
    printf("Error parsing options. Error at option %s.\n",argv[optup]);
  }

  double refinementQuality = strtod(refinementQualityString, NULL);
  cout << "Refinement quality is: " << refinementQuality << std::endl;

  double alpha = strtod(alphaString, NULL);
  cout << "Alpha is: " << alpha << std::endl;

  double minDihedral = strtod(minDihedralString, NULL);
	cout << "Minimal dihedral is: " << minDihedral << endl;

	double maxTimeSeconds = strtod(maxTimeSecondsString, NULL);

  ObjMesh * inputSurfaceMesh = new ObjMesh(inputSurfaceMeshFilename);

  printf("Running the tet mesher...\n");fflush(NULL);
  TetMesh * tetMesh = TetMesher().compute(inputSurfaceMesh, refinementQuality, alpha, minDihedral, maxSteinerVertices, maxTimeSeconds);
  
  if (tetMesh == NULL)
  {
    printf("Error running the mesher.\n");
  }
  else
  {
    printf("Saving the output mesh to %s.\n", ouputTetMeshFilename);
    tetMesh->save(ouputTetMeshFilename);
  }
  delete tetMesh;
  delete inputSurfaceMesh;
  return(0);
}

