/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 3.0                               *
 *                                                                       *
 * "distance field" driver , Copyright (C) 2007 CMU, 2016 USC            *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Jernej Barbic, Hongyi Xu, Yijing Li                     *
 * http://www.jernejbarbic.com/code                                      *
 *                                                                       *
 * Research: Jernej Barbic, Hongyi Xu, Doug L. James                     *
 *                                                                       *
 * Funding: National Science Foundation, Link Foundation,                *
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

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <fstream>
#include <iostream>
#include <float.h>
#include <vector>
using namespace std;

#include "distanceFieldCreator.h"
#include "closestPointField.h"
#include "objMesh.h"
#include "getopts.h"

/*
  Computes signed or unsigned distance field for a given triangle mesh,
  optionally using the pipeline from:

  Hongyi Xu, Jernej Barbic:
  Signed Distance Fields for Polygon Soup Meshes, Graphics Interface 2014
  Montreal, Canada
*/

int main( int argc, char** argv )
{
  bool printDetailedHelp = (argc >= 2) && (strcmp(argv[1],"-h") == 0);

  if ( argc < 5 ) 
  {
    cout << "Computes signed or unsigned distance field for a given triangle mesh. Use -h for more help." << endl;
    cout << "Usage: " << argv[0] << " [-h<help>] [obj file] [resolutionX] [resolutionY] [resolutionZ] "
        "[-n<narrow band>] [-s<signed field>] [-o<output file>] [-T<num threads>] [-m<signed field mode>] "
        "[-b<xmin,ymin,zmin,xmax,ymax,zmax>] [-c<cube box>] [-e<box expansion ratio>] "
        "[-d<max octree depth>] [-t<max #triangles per octree cell>] "
        "[-w<band width>] "
        "[-g<sigma>] [-r<do not subtract sigma>] [-i<precomputed unsigned field>] "
        "[-v<output voronoi diagram file>] [-p<compute closest filed>] " << endl;
    if (!printDetailedHelp)
      return 1;
  }

  if (printDetailedHelp)
  {
    cout << "  -n: compute only narrow band distance field; default: false" << endl;
    cout << "  -s: compute signed distance field (requires manifold surface with no boundary); default: unsigned" << endl;
    cout << "  -o: specify output file; default: out.dist " << endl;
    cout << "  -T: num of threads used to compute distance field (not used in narrow band)" << endl;
    cout << "  -m: signed field computation mode: " << endl;
    cout << "      0: BASIC assume the input obj mesh is manifold and self-intersection-free " << endl;
    cout << "      1: POLYGONSOUP handle non-manifold and/or self-intersecting meshes, using the SignedDistanceFieldFromPolygonSoup pipeline,\n"
            "         as published in:\n"
            "         Hongyi Xu, Jernej Barbic:\n"
            "         Signed Distance Fields for Polygon Soup Meshes, Graphics Interface 2014, "
        "Montreal, Canada" << endl;
    cout << "      default: BASIC" << endl;
    cout << "      2: AUTO use BASIC if mesh is manifold, otherwise POLYGONSOUP" << endl;

    cout << "  ============== Bounding Box =============" << endl;
    cout << "  -b: specify scene bounding box; default: automatic box with cubic cells" << endl;
    cout << "  -c: force bounding box into a cube; default: automatic box with cubic cells " << endl;
    cout << "  -e: expansion ratio for box (when automatic); default: 1.5" << endl;

    cout << "  ================= Octree ================" << endl;
    cout << "  -d: max octree depth to use (default: 10)" << endl;
    cout << "  -t: max num triangles per octree cell (default: 15)" << endl;

    cout << "  =============== Narrow Band =============" << endl;
    cout << "  -w: the band width for narrow band distance field" << endl;

    cout << "  ========= Polygon Soup Pipeline =========" << endl;
    cout << "  -g: sigma value" << endl;
    cout << "  -r: do not substract sigma when creating signed field (default: false)" << endl;
    cout << "  -i: precomputed unsigned field for creating signed field (default: none)" << endl;

    cout << "  ================ Extra ==================" << endl;
    cout << "  -v: also compute voronoi diagram (defaut: not computed)" << endl;
    cout << "  -p: also compute closest points (defaut: not computed)" << endl;

    cout << "  ***  Output file is binary. Format: " << endl;
    cout << "    - resolutionX,resolutionY,resolutionZ (three signed 4-byte integers (all equal), 12 bytes total)" << endl;
    cout << "    - bminx,bminy,bminz (coordinates of the lower-left-front corner of the bounding box: (three double precision 8-byte real numbers , 24 bytes total)" << endl;
    cout << "    - bmaxx,bmaxy,bmaxz (coordinates of the upper-right-back corner of the bounding box: (three double precision 8-byte real numbers , 24 bytes total)" << endl;
    cout << "    - distance data (in single precision; data alignment: [0,0,0],...,[resolutionX,0,0],[0,1,0],...,[resolutionX,resolutionY,0],[0,0,1],...,[resolutionX,resolutionY,resolutionZ]; total num bytes: sizeof(float)*(resolutionX+1)*(resolutionY+1)*(resolutionZ+1))" << endl;

    return 0;
  }

  int resolutionX = strtol(argv[2], NULL, 10);
  int resolutionY = strtol(argv[3], NULL, 10);
  int resolutionZ = strtol(argv[4], NULL, 10);
  if ((resolutionX <= 0) || (resolutionY <= 0) || (resolutionZ <= 0))
  {
    printf("Invalid resolution.\n");
    return 1;
  }

  char * objMeshname = argv[1];

  bool doublePrecision = false;
  bool narrowBandField = false;
  bool signedField = false;
  char outputFile[4096] = "out.dist";
  int numThreads = 0;
  int mode = 0;

  bool useCubeBox = false;
  char bboxString[4096] = "__default";
  double expansionRatio = 1.5;
  char expansionRatioString[4096] = "__none";

  int maxTriCount = 15;
  int maxDepth = 10;

  double bandWidth = 5;
  char bandWidthString[4096] = "__none";

  char precomputedUnsignedDistanceFieldFilename[4096] = "__none";
  double sigma = 0;
  char sigmaString[1028] = "__none";
  bool nonSubtractSigma = false;

  char outputVoronoiDiagramFile[4096] = "__none";
  bool closestPointField = false;

  opt_t opttable[] =
  {
    { "n", OPTBOOL, &narrowBandField},
    { "m", OPTINT, &mode },
    { "w", OPTSTR, bandWidthString},
    { "b", OPTSTR, bboxString },
    { "r", OPTBOOL, &nonSubtractSigma},
    { "c", OPTBOOL, &useCubeBox },
    { "d", OPTINT, &maxDepth },
    { "e", OPTSTR, expansionRatioString},
    { "o", OPTSTR, outputFile },
    { "p", OPTBOOL, &closestPointField },
    { "s", OPTBOOL, &signedField },
    { "t", OPTINT, &maxTriCount },
    { "T", OPTINT, &numThreads },
    { "v", OPTSTR, &outputVoronoiDiagramFile },
    { "i", OPTSTR, precomputedUnsignedDistanceFieldFilename},
    { "g", OPTSTR, sigmaString },
    { NULL, 0, NULL }
  };

  argv += 4;
  argc -= 4;
  int optup = getopts(argc,argv,opttable);
  if (optup != argc)
  {
    printf("Error parsing options. Error at option %s.\n",argv[optup]);
    return 1;
  }

  if (mode < 0 || mode > 2)
  {
    printf("Invalid sign field computation mode: %d.\n", mode);
    return 1;
  }

  if (strcmp(expansionRatioString, "__none") != 0)
    expansionRatio = strtod(expansionRatioString, NULL);
  if (expansionRatio < 1.0)
  {
    printf("Invalid expansion ratio: %G.\n", expansionRatio);
    return 1;
  }
  printf("Expansion ratio is: %G.\n", expansionRatio);

  if (strcmp(sigmaString, "__none") != 0)
  {
    sigma = strtod(sigmaString, NULL);
    if (sigma <= 0.0)
    {
      printf("Invalid sigma: %G.\n", sigma);
      return 1;
    }
    printf("sigma is: %G.\n", sigma);
  }
  else if (mode == 1 && signedField) // no sigma given
  {
    printf("No sigma given when computing signed filed with POLYGON mode.\n");
    return 1;
  }

  if (narrowBandField)
  {
    if (strcmp(bandWidthString, "__none") != 0)
    {
      bandWidth = strtod(bandWidthString, NULL);
      if (bandWidth <= 0)
      {
        printf("Invalid narrow band width: %G.\n", bandWidth);
        return 1;
      }

      if (mode == 1 && signedField && bandWidth <= sigma)
      {
        printf("Error: band width %G, smaller than sigma %G.\n", bandWidth, sigma);
        return 1;
      }
    }
    else
    {
      printf("Error: No band width is provided.\n");
      return 1;
    }
  }

  ObjMesh * objMesh = NULL;
  try
  {
    objMesh = new ObjMesh(objMeshname);
  }
  catch(int exceptionCode)
  {
    printf("Error: could not load obj file %s (exception code %d).\n", objMeshname, exceptionCode);
    exit(1);
  }

  printf("Data output format is: IEEE %s precision.\n", doublePrecision ? "double (64-bit)" : "single (32-bit)");
  Vec3d * bmin = NULL;
  Vec3d * bmax = NULL;
  if (strcmp(bboxString,"__default") != 0)
  {
    printf("bbox string: %s\n", bboxString);
    double xmin,ymin,zmin,xmax,ymax,zmax;
    sscanf(bboxString,"%lf,%lf,%lf,%lf,%lf,%lf",&xmin,&ymin,&zmin,&xmax,&ymax,&zmax);
    bmin = new Vec3d(xmin,ymin,zmin);
    bmax = new Vec3d(xmax,ymax,zmax);
    cout << "Scene bounding set to: " << bmin << " " << bmax << endl;
    if ((xmin >= xmax) || (ymin >= ymax) || (zmin >= zmax))
    {
      cout << "Invalid bounding box.\n";
      return 1;
    }
  }

  DistanceFieldCreator * distanceFieldCreator = new DistanceFieldCreator(objMesh, expansionRatio, useCubeBox, bmin, bmax);

  char * inputUnsignedFieldFilename = NULL;
  if (strcmp(precomputedUnsignedDistanceFieldFilename, "__none") != 0)
    inputUnsignedFieldFilename = precomputedUnsignedDistanceFieldFilename;

  DistanceFieldCreator::SignedFieldCreationMode creatorMode;
  if (mode == 0)
    creatorMode = DistanceFieldCreator::BASIC;
  else if (mode == 1)
    creatorMode = DistanceFieldCreator::POLYGONSOUP;
  else
    creatorMode = DistanceFieldCreator::AUTO;

  if (narrowBandField)
  {
    DistanceFieldNarrowBand * field = distanceFieldCreator->ComputeDistanceFieldNarrowBand(resolutionX, resolutionY, resolutionZ,
        bandWidth, signedField, creatorMode, sigma, !nonSubtractSigma, maxTriCount, maxDepth, inputUnsignedFieldFilename);
    cout << "Saving the distance field to " << outputFile << " ." << endl;
    field->save(outputFile, doublePrecision);
  }
  else
  {
    bool computeVoronoiDiagram = (strcmp(outputVoronoiDiagramFile, "__none") != 0);
    DistanceField * field = distanceFieldCreator->ComputeDistanceField(resolutionX, resolutionY, resolutionZ,
        signedField, creatorMode, sigma, !nonSubtractSigma, numThreads, computeVoronoiDiagram, maxTriCount, maxDepth, closestPointField,
        inputUnsignedFieldFilename);
    cout << "Saving the distance field to " << outputFile << " ." << endl;
    field->save(outputFile, doublePrecision);

    if (computeVoronoiDiagram)
    {
      cout << "Saving Voronoi diagram to " << outputVoronoiDiagramFile << " ." << endl;
      field->saveVoronoiDiagram(outputVoronoiDiagramFile); 
    }
  }
 
  return(0);
}

