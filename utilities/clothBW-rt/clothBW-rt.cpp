/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 3.0                               *
 *                                                                       *
 * "Cloth simulator" application,                                        *
 * Copyright (C) 2007 CMU, 2009 MIT, 2016 USC                            *
 *                                                                       *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Andy Pierce, Yijing Li, Yu Yu Xu, Jernej Barbic         *
 * http://www.jernejbarbic.com/code                                      *
 *                                                                       *
 * Funding: National Science Foundation                                  *
 *          Zumberge Research and Innovation Fund at USC                 *
 *                                                                       *
 * This utility is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this library in the file LICENSE.txt                    *
 *                                                                       *
 * This utility is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

#include <math.h>
#include <vector>
#include "GL/glui.h"
#include "initGraphics.h"
#include "performanceCounter.h"
#include "objMesh.h"
#include "clothBWMT.h"
#include "clothBWFromObjMesh.h"
#include "clothBWForceModel.h"
#include "implicitNewmarkSparse.h"
#include "implicitBackwardEulerSparse.h"
#include "configFile.h"
#include "loadList.h"
#include "sceneObjectDeformable.h"
#include "saveScreenShot.h"

/*
  Driver for cloth simulation. Cloth is implemented using:

  Baraff and Witkin: Large Steps in Cloth Simulation, SIGGRAPH 1998. 
  David Pritchard: Implementing Baraff & Witkin's Cloth Simulation, May 2003, with minor updates April 2012
*/

#define Mm_PI 3.1415926
#define  MAX_FILE 4096

// simulation parameters
float timeStep;               //default: 0.01
float dampingMass;            //default: 0.00
float dampingStiffness;       //default: 0.001
float surfaceDensity;         //default: 0.25;
float tensileStiffness;       //default: 8500.0; 
float shearStiffness;         //default: 100.0;
float bendStiffnessU;         //default: 0.0002;
float bendStiffnessV;         //default: 0.0002;
float mouseForceCoeff;         //default: 0.05;
float wind_x;                 //default: 0.0;
float wind_y;                 //default: 0.0;
float wind_z;                 //default: 0.0;
float gravityForce;           //default: 9.81
int addGravity;               //default: 1;
int useRestAngles;            //default: 1

// computation parameters (default to computing everything)
int computeStretchShearForce, computeStretchShearStiffness, computeBendForce, computeBendStiffness;
int numInternalForceThreads; // default = 0 (not multi-threaded)
int numSolverThreads; // default = 0 (not multi-threaded)
int renderFixedVertices = 1;

void initScene();

//glui
GLUI * glui = NULL;
GLUI_StaticText * systemSolveStaticText;
GLUI_StaticText * forceAssemblyStaticText;

int loadNewConfigFile = 0;

// graphics
char windowTitleBase[4096] = "Cloth Simulator";
int windowID;
int windowWidth = 800;
int windowHeight = 600;

//interactive control
double zNear = 0.01;               //default: 0.01
double zFar = 10.0;                //default:10.0;
double cameraRadius;
double focusPositionX, focusPositionY, focusPositionZ;
double cameraLongitude, cameraLatitude;
int g_iMenuId;			// mouse activity
int g_vMousePos[2];
int g_iLeftMouseButton,g_iMiddleMouseButton,g_iRightMouseButton;
double forceAssemblyTime = 0.0;
double systemSolveTime = 0.0;

// start out paused, wire-frame view, scene unlocked (you can drag to add forces)
int runSimulation=0, renderWireframe=1, saveScreenToFile=0, dragForce = 0, pulledVertex = -1, lockScene = 0, axis = 0, pin = 0, sprite=0, renderNormals = 0, displayMenu = 0, useTextures = 1;

int shiftPressed=0;
int altPressed=0;
int ctrlPressed=0;

int dragStartX, dragStartY;
int loadOption = 0;
int cloth_width = 30;
int cloth_length =45;
std::vector<int> pin_points;
int graphicsFrame = 0;

//simulation objects
ClothBW * clothBW = NULL;
ClothBWMT * clothBWMT = NULL;

// rendering
SceneObjectDeformable * sceneObjDeform = NULL;
SceneObject * extraSceneGeometry = NULL;
Lighting * light = NULL;

// solver
IntegratorBase * integratorBase = NULL;
ImplicitNewmarkSparse * implicitNewmarkSparse = NULL;
ImplicitBackwardEulerSparse * implicitBackwardEulerSparse = NULL;
ForceModel * forceModel = NULL;
ClothBWForceModel * clothBWForceModel = NULL;
SparseMatrix * massMatrix = NULL;

// camera
SphericalCamera * camera = NULL;

// load cloth from .obj
ObjMesh * objMesh = NULL;
ClothBWFromObjMesh * clothBWFromObjMesh = NULL;

PerformanceCounter timer;
double FPS = 0.0; 

PerformanceCounter display_timer;
double display_time = 0.0;
double dotimestep_time = 0.0;

int * fixedVertices;
int numFixedVertices;

int n;						// record number of particles
double * u = NULL;			// record current deformation
double * f_ext = NULL;		// record current external forces

// files
char configFilename[MAX_FILE];
char fixedVerticesFilename[MAX_FILE];
char objMeshname[MAX_FILE] = "skirt.obj";
char lightingFilename[MAX_FILE];
char extraSceneGeometryFilename[MAX_FILE];

// This function specifies parameters (and default values, if applicable) for the
// configuration file. It then loads the config file and parses the options contained
// within. After parsing is complete, a list of parameters is printed to the terminal.
void initConfigurations()
{
  printf("Parsing configuration file %s...\n", configFilename);
  ConfigFile configFile;
  
  // specify the entries of the config file
  
  // get camera info (optional)
  configFile.addOptionOptional("focusPositionX", &focusPositionX, 0.0);
  configFile.addOptionOptional("focusPositionY", &focusPositionY, 10.0);
  configFile.addOptionOptional("focusPositionZ", &focusPositionZ, 0.0);
  configFile.addOptionOptional("cameraRadius", &cameraRadius, 6.0);
  configFile.addOptionOptional("cameraLongitude", &cameraLongitude, -10.0);
  configFile.addOptionOptional("cameraLatitude", &cameraLatitude, 45.0);
  configFile.addOptionOptional("zBufferNear", &zNear, 0.01);
  configFile.addOptionOptional("zBufferFar", &zFar, 10.0);
  configFile.addOptionOptional("renderWireframe", &renderWireframe, renderWireframe);
  
  // stiffness and damping parameters (optional)
  configFile.addOptionOptional("tensileStiffness", &tensileStiffness, 6500.0f);
  configFile.addOptionOptional("shearStiffness", &shearStiffness, 100.0f);
  configFile.addOptionOptional("bendStiffnessU", &bendStiffnessU, 0.0002f);
  configFile.addOptionOptional("bendStiffnessV", &bendStiffnessV, 0.0002f);
  configFile.addOptionOptional("dampingMass", &dampingMass, 0.00f);
  configFile.addOptionOptional("dampingStiffness", &dampingStiffness, 0.001f);
  configFile.addOptionOptional("surfaceDensity", &surfaceDensity, 0.25f);
  
  // force parameters and gravity (optional)
  configFile.addOptionOptional("mouseForceCoeff", &mouseForceCoeff, 0.05f); // how powerful your mouse is
  configFile.addOptionOptional("wind_x", &wind_x, 0.0f); // default is no wind
  configFile.addOptionOptional("wind_y", &wind_y, 0.0f);
  configFile.addOptionOptional("wind_z", &wind_z, 0.0f);
  configFile.addOptionOptional("gravityForce", &gravityForce, 9.81f);
  configFile.addOptionOptional("addGravity", &addGravity, 1); // 1 for gravity, 0 for no gravity
  
  // computation parameters
  configFile.addOptionOptional("computeStretchShearForce", &computeStretchShearForce, 1);
  configFile.addOptionOptional("computeStretchShearStiffness", &computeStretchShearStiffness, 1);
  configFile.addOptionOptional("computeBendForce", &computeBendForce, 1);
  configFile.addOptionOptional("computeBendStiffness", &computeBendStiffness, 1);
  
  configFile.addOptionOptional("useRestAngles", &useRestAngles, 1);
  
  //set to 0 to disable multi-threading (or 1 to only launch a single thread...)
  configFile.addOptionOptional("numInternalForceThreads", &numInternalForceThreads, 0);
  configFile.addOptionOptional("numSolverThreads", &numSolverThreads, 0);
  
  // timestep (optional)
  configFile.addOptionOptional("timeStep", &timeStep, 0.02f);
  
  // ===== files =====
  
  // need .obj file
  configFile.addOption("objMeshname", objMeshname);
  
  // need file for fixed vertices
  configFile.addOption("fixedVerticesFilename", fixedVerticesFilename);
  
  // optional file for lighting
  //configFile.addOptionOptional("lightingFilename", lightingFilename, "/Users/Andy/Downloads/lighting-v1.0/example.lighting");
  configFile.addOption("lightingFilename", lightingFilename);

  // statis scene geometry (optional)
  configFile.addOptionOptional("extraSceneGeometryFilename", extraSceneGeometryFilename, "__none");
  
  // parse the configuration file
  if (configFile.parseOptions(configFilename) != 0)
  {
    printf("Error parsing options.\n");
    exit(1);
  }
  
  // the config variables have now been loaded with their specified values
  
  // informatively print the variables (with assigned values) that were just parsed
  configFile.printOptions();
}
  
void Sync_GLUI()
{
  glui->sync_live();
}

void drawString(const char * str) 
{
  glPushAttrib(GL_LIGHTING_BIT | GL_CURRENT_BIT); // lighting and color mask
  glDisable(GL_LIGHTING);     // need to disable lighting for proper text color
  
  glColor3f(1.0, 1.0, 1.0); // set text color
  
  // loop all characters in the string
  while(*str)
  {
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, *str);
    ++str;
  }
  
  glEnable(GL_LIGHTING);
  glPopAttrib();
}

void display_menu(double fps)
{
  glDisable(GL_LIGHTING);
  glColor3f(1.0,0.0,0.0);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix(); //push
  glLoadIdentity();
  glMatrixMode(GL_PROJECTION);
  glPushMatrix(); //push
  glLoadIdentity();
  gluOrtho2D (0, windowWidth, 0, windowHeight);
  
  int dist_small = 15;
  int dist_big = 25;
  int base = 15;
 
  glRasterPos2f(0,windowHeight-base); 
  char fps_string[20];
  sprintf(fps_string,"FPS:%G",fps);
  drawString(fps_string);
  
  
  base = base + dist_big;
  glRasterPos2f(0,windowHeight-base); 
  char simulation_stat[50];
  if(runSimulation == 0)
    sprintf(simulation_stat,"Simulation[p]: Paused");
  else
    sprintf(simulation_stat,"Simulation[p]: Running");
  drawString(simulation_stat);
  
  base = base + dist_big;
  glRasterPos2f(0,windowHeight-base); 
  char scene_lock_status[100];
  if(lockScene == 1)
  {
    sprintf(scene_lock_status,"Scene Status[l]:Locked! Cannot add interactive force.");
    drawString(scene_lock_status);
  }
  else
  {
    sprintf(scene_lock_status,"Scene Status[l]:Unlocked!");
    drawString(scene_lock_status);
    if(runSimulation == 1)
    {
      base = base + dist_small;
      glRasterPos2f(20,windowHeight-base); 
      sprintf(scene_lock_status,"Press left mouse button to select a vertex");
      drawString(scene_lock_status);
      base = base + dist_small;
      glRasterPos2f(20,windowHeight-base); 
      sprintf(scene_lock_status,"Press right mouse button to finish adding the force");
      drawString(scene_lock_status);
    }
    else
    {
      base = base + dist_small;
      glRasterPos2f(20,windowHeight-base); 
      sprintf(scene_lock_status,"Cannot add interactive force unless running the simulation!");
      drawString(scene_lock_status);			
    }
  }
  
  base = base + dist_big;
  glRasterPos2f(0,windowHeight-base); 
  char axis_status[50];
  if(axis == 1)
  {
    sprintf(axis_status,"Axis Display[a]:ON");
    drawString(axis_status);
  }
  else
  {
    sprintf(axis_status,"Axis Display[a]:OFF");
    drawString(axis_status);
  }
  
  base = base + dist_big;
  glRasterPos2f(0,windowHeight-base); 
  char constraint_status[100];
  if(pin == 0)
  {
    sprintf(constraint_status,"Constraint[c]:DEFAULT");
    drawString(constraint_status);		
  }
  else
  {
    sprintf(constraint_status,"Constraint[c]:RANDOM");
    drawString(constraint_status);
    base = base + dist_small;
    glRasterPos2f(20,windowHeight-base); 
    sprintf(constraint_status,"Click On the Vertex you want to PIN, Reset and RUN!");
    drawString(constraint_status);	
  }
  
  // adding force/stiffness matrix toggling functionality
  base = base + dist_big;
  glRasterPos2f(0,windowHeight-base); 
  
  char stshStatus[100];
  char bendStatus[100];
  
  // first forces
  
  // stretch & shear
  if (computeStretchShearForce)
    sprintf(stshStatus, "Stretch/Shear Force[s]: ENABLED");
  else
    sprintf(stshStatus, "Stretch/Shear Force[s]: DISABLED");
  drawString(stshStatus);
  
  base = base + dist_small;
  glRasterPos2f(0,windowHeight-base);
  
  // bend
  if (computeBendForce)
    sprintf(bendStatus, "Bend Force[b]: ENABLED");
  else
    sprintf(bendStatus, "Bend Force[b]: DISABLED");
  drawString(bendStatus);
  
  base = base + dist_small;
  glRasterPos2f(0,windowHeight-base);
  
  // now stiffness matrices
  
  // stretch & shear
  if (computeStretchShearStiffness)
    sprintf(stshStatus, "Stretch/Shear Stiffness[d]: ENABLED");
  else
    sprintf(stshStatus, "Stretch/Shear Stiffness[d]: DISABLED");
  drawString(stshStatus);
  
  base = base + dist_small;
  glRasterPos2f(0,windowHeight-base);
  
  // bend
  if (computeBendStiffness)
    sprintf(bendStatus, "Bend Stiffness[n]: ENABLED");
  else
    sprintf(bendStatus, "Bend Stiffness[n]: DISABLED");
  drawString(bendStatus);
  
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glEnable(GL_LIGHTING);
}

// this function does the following:
// (1) Clears display
// (2) Points camera at scene
// (3) Draws axes (if applicable) and sphere surrounding scene
// (4) Displays GUI menu
// (5) Sets cloth deformations
// (6) Sets lighting conditions
// (7) Builds surface normals for cloth (take this out to increase performance)
// (8) Renders cloth
// (9) Render pulled vertex in different color (if applicable)
void displayFunction()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
  glMatrixMode(GL_MODELVIEW); 
  glLoadIdentity();	
  camera->Look();
  glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
  glStencilFunc(GL_ALWAYS, 0, ~(0u));

  // render any extra scene geometry
  glStencilFunc(GL_ALWAYS, 0, ~(0u));
  if (extraSceneGeometry != NULL)
    extraSceneGeometry->Render();

  if(axis)
  {
    glDisable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);

    glBegin(GL_LINES);
    for (int i = 0; i < 3; i++)
    {
      float color[3] = { 0, 0, 0 };
      color[i] = 1.0;
      glColor3fv(color);

      float vertex[3] = {0, 0, 0};
      vertex[i] = 1.0;
      glVertex3fv(vertex);
      glVertex3f(0, 0, 0);
    }
    glEnd();
    glEnable(GL_LIGHTING);
  }
 
  if (displayMenu)
    display_menu(FPS);
  
  sceneObjDeform->SetVertexDeformations(u);
  
  // render cloth
  if (clothBW != NULL)
  { 
    glLineWidth(1.0);
    glStencilFunc(GL_ALWAYS, 1, ~(0u));
    //glEnable(GL_COLOR_MATERIAL);
    
    //sceneObjDeform->SetVertexDeformations(u);
    sceneObjDeform->SetLighting(light);
    
    //sceneObjDeform->SetUpTextures(textureMode);
    //sceneObjDeform->BuildNormals();
    sceneObjDeform->BuildNormalsFancy();

    if (renderNormals)
    {
      glDisable(GL_LIGHTING);
      glColor3f(0,0,1);
      sceneObjDeform->RenderNormals();
      glEnable(GL_LIGHTING);
    }

    // render fixed vertices
    glDisable(GL_LIGHTING);
    if (renderFixedVertices)
    {
      for(int i=0; i<numFixedVertices; i++)
      {
        glColor3f(1,0,0);
        double fixedVertexPos[3];
        sceneObjDeform->GetSingleVertexRestPosition(fixedVertices[i],
            &fixedVertexPos[0], &fixedVertexPos[1], &fixedVertexPos[2]);

        glEnable(GL_POLYGON_OFFSET_POINT);
        glPolygonOffset(-1.0,-1.0);
        glPointSize(12.0);
        glBegin(GL_POINTS);
        glVertex3f(fixedVertexPos[0], fixedVertexPos[1], fixedVertexPos[2]);
        glEnd();
        glDisable(GL_POLYGON_OFFSET_FILL);
      }
    }
    
    glEnable(GL_LIGHTING);
    sceneObjDeform->Render();

    if (renderWireframe)
    {
      glDisable(GL_LIGHTING);
      glColor3f(0,0,0);
      sceneObjDeform->RenderEdges();
      glEnable(GL_LIGHTING);
    }
  }
  else
  {
    printf("Error: cloth object has not been created.\nExiting...\n");
    exit(1);
  }
  
  // if pulling a vertex, render that vertex in diff color
  if (pulledVertex != -1)
  {
    glDisable(GL_LIGHTING);
    glDisable(GL_TEXTURE_2D);
    
    //sceneObjDeform->HighlightVertex(pulledVertex);
    
    glColor3f(1.0, 0.0, 0.0);
    glPointSize(8.5);
    glBegin(GL_POINTS);

    const Vec3d * restPos = clothBW->GetRestPositions();
    double vertexCoords[3];
    vertexCoords[0] = restPos[pulledVertex][0] + u[3*pulledVertex + 0];
    vertexCoords[1] = restPos[pulledVertex][1] + u[3*pulledVertex + 1];
    vertexCoords[2] = restPos[pulledVertex][2] + u[3*pulledVertex + 2];
   
    //printf("Vertex: (%G, %G, %G)\n", vertexCoords[0], vertexCoords[1], vertexCoords[2]);
    
    glVertex3f(vertexCoords[0], vertexCoords[1], vertexCoords[2]);
    
    glEnd();
    glEnable(GL_TEXTURE_2D);
    glEnable(GL_LIGHTING);
  }
  
  glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
  glutSwapBuffers();
}

// this function does the following:
// (1) Takes a screenshot (if applicable)
// (2) Updates & displays the fps
// (3) Applies mouse forces
// (4) Applies wind forces
// (5) Sets computation mode for clothBW
// (6) Performs a timestep
void idleFunction(void)
{
  static int timeStepCount = 0;
  float timestepsPerSecond = 1.0/timeStep;
  
  glutSetWindow(windowID);
  
  // recording
  char s[70]="picxxxx.ppm\0";
  
  s[40] = 48 + (sprite / 1000);
  s[41] = 48 + (sprite % 1000) / 100;
  s[42] = 48 + (sprite % 100 ) / 10;
  s[43] = 48 + sprite % 10;
  
  if (saveScreenToFile==1 && timeStepCount<timestepsPerSecond*10)
  {
    Screenshot::SaveScreenshot(s, ImageIO::FORMAT_PNG, windowWidth, windowHeight);
    //saveScreenToFile=0; // save only once, change this if you want continuous image generation (i.e. animation)
    sprite++;
  }
  if (timeStepCount>timestepsPerSecond*10 && saveScreenToFile)
  {
    printf("DONE!\n");
    runSimulation = 0;
    saveScreenToFile = 0;
  }
  if (sprite >= 300) // allow only 300 snapshots
    exit(0);	
  
  // external force & time step
  if (runSimulation == 1) // if unpaused
  {
    // clear old external force
    integratorBase->SetExternalForcesToZero();
    
    // external forces
    
    // user mouse forces
    if((g_iLeftMouseButton && pulledVertex != -1 && lockScene == 0))
    {
      double forceX = (g_vMousePos[0] - dragStartX);
      double forceY = -(g_vMousePos[1] - dragStartY);
      
      //add external force here
      double externalForce[3];
      camera->CameraVector2WorldVector_OrientationOnly3D(forceX, forceY, 0, externalForce);	
      for(int i = 0 ; i < 3; i++)
        externalForce[i] *= mouseForceCoeff;
      //printf("fx: %G fy: %G | %G %G %G\n", forceX, forceY, externalForce[0], externalForce[1], externalForce[2]);
      
      memset(f_ext, 0.0, 3*n);
      for(int i = 0 ; i < 3; i++)
        f_ext[3*pulledVertex+i] += externalForce[i];
      
      integratorBase->AddExternalForces(f_ext);
    }

    // wind forces
    for(int i = 0 ; i < clothBW->GetNumVertices(); i++)
    {
      f_ext[3*i+0] = wind_x;
      f_ext[3*i+1] = wind_y;
      f_ext[3*i+2] = wind_z;
    }
    integratorBase->AddExternalForces(f_ext);
    
    // adding force/stiffness matrix toggling ability
    bool parameters[4]; 
    parameters[0] = computeStretchShearForce;
    parameters[1] = computeBendForce;
    parameters[2] = computeStretchShearStiffness;
    parameters[3] = computeBendStiffness;
    
    clothBW->SetComputationMode(parameters);
    clothBW->UseRestAnglesForBendingForces(useRestAngles);
    
    integratorBase->DoTimestep(); // the big timestep
    timeStepCount++;

    memcpy(u, integratorBase->Getq(), sizeof(double) * 3 * n);
    
    sceneObjDeform->BuildFaceNormals();
  }

  // fps
  timer.StopCounter();
  double elapsedTime = timer.GetElapsedTime();  
  if (elapsedTime >= 1.0 / 4)  
  {
    timer.StartCounter(); 
    FPS = graphicsFrame / elapsedTime;

    forceAssemblyTime = implicitNewmarkSparse->GetForceAssemblyTime();
    systemSolveTime = implicitNewmarkSparse->GetSystemSolveTime();

    char ptext[96];
    sprintf(ptext, "Force assembly: %G", forceAssemblyTime);
    forceAssemblyStaticText->set_text(ptext);
    sprintf(ptext, "System solve: %G", systemSolveTime);
    systemSolveStaticText->set_text(ptext);
    Sync_GLUI();

    graphicsFrame = 0;
  }
 
  graphicsFrame++;
  
  glutPostRedisplay();
}

void reshape(int x,int y)
{
  glViewport(0,0,x,y);
  
  windowWidth = x;
  windowHeight = y;
  
  glMatrixMode(GL_PROJECTION); // Select The Projection Matrix
  glLoadIdentity(); // Reset The Projection Matrix
  
  // gluPerspective(90.0,1.0,0.01,1000.0);
  gluPerspective(60.0f, 1.0 * windowWidth / windowHeight, zNear, zFar);
  
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

void exit_buttonCallBack(int code)
{
  // free memory
  
  // cloth
  if (clothBW != NULL)
    delete clothBW;
  if (clothBWMT != NULL)
    delete clothBWMT;
  if (clothBWFromObjMesh != NULL)
    delete clothBWFromObjMesh;
  
  // solver
  if (implicitBackwardEulerSparse != NULL)
    delete implicitBackwardEulerSparse;
  
  // model
  if (clothBWForceModel != NULL)
    delete clothBWForceModel;
  
  // render objs
  if (objMesh != NULL)
    delete objMesh;
  if (sceneObjDeform != NULL)
    delete sceneObjDeform;
  if (light != NULL)
    delete light;
  
  // other
  if (u)
    free(u);
  if (f_ext)
    free(f_ext);
  
  exit(0);
}

void keyboardFunction(unsigned char key, int x, int y)
{
  switch (key)
  {
    case 27:
      exit_buttonCallBack(0);
      break;
      
    case 'a':
      axis = 1 - axis;
      break;
      
    case '\\':
      camera->Reset();
      break;
      
    case 'w':
      renderWireframe = !renderWireframe;
      break;

    case 't':
      useTextures = !useTextures;
      if (useTextures)
        sceneObjDeform->EnableTextures();
      else
        sceneObjDeform->DisableTextures();
      break;
      
    case 'p':
      runSimulation = 1 - runSimulation;
      break;
      
    case 'l':
      lockScene = 1 - lockScene;
      pulledVertex = -1;
      break;
      
    case 'i':
    {
      double focusPos[3];
      camera->GetFocusPosition(focusPos);
      double cameraX,cameraY,cameraZ;
      camera->GetAbsWorldPosition(cameraX,cameraY,cameraZ);
      printf("Camera is positioned at: %G %G %G\n", cameraX,cameraY,cameraZ);
      printf("Camera radius is: %G\n", camera->GetRadius());
      printf("Camera Phi is: %G\n", 180.0/M_PI*camera->GetPhi());
      printf("Camera Theta is: %G\n", 180.0/M_PI*camera->GetTheta());
      printf("Camera focus is: %G %G %G\n", focusPos[0], focusPos[1], focusPos[2]);
    }

    case 'c':
      pin = 1 - pin;
      break;
      
    case ' ':
      saveScreenToFile = 1 - saveScreenToFile;
      break;
      
    case 's':
      computeStretchShearForce = !computeStretchShearForce;
      printf("Computing stretch/shear force is now: %s.\n", computeStretchShearForce ? "ON" : "OFF");
      break;

    case 'S':
      computeStretchShearStiffness = !computeStretchShearStiffness;
      printf("Computing stretch/shear stiffness matrix is now: %s.\n", computeStretchShearStiffness ? "ON" : "OFF");
      break;

    case 'b':
      computeBendForce = !computeBendForce;
      printf("Computing bend force is now: %s.\n", computeBendForce ? "ON" : "OFF");
      break;

    case 'B':
      computeBendStiffness = !computeBendStiffness;
      printf("Computing bend stiffness matrix is now: %s.\n", computeBendStiffness ? "ON" : "OFF");
      break;

    case 'F':
      renderFixedVertices = !renderFixedVertices;
      break;

    case 'm':
      displayMenu = !displayMenu;
      break;

    case 'N':
      renderNormals = !renderNormals;
      break;
  }	
}

// reacts to pressed "special" keys
void specialFunction(int key, int x, int y)
{
  switch (key)
  {
    case GLUT_KEY_LEFT:
      camera->MoveFocusRight(+0.1 * camera->GetRadius());
    break;

    case GLUT_KEY_RIGHT:
      camera->MoveFocusRight(-0.1 * camera->GetRadius());
    break;

    case GLUT_KEY_DOWN:
      camera->MoveFocusUp(+0.1 * camera->GetRadius());
    break;

    case GLUT_KEY_UP:
      camera->MoveFocusUp(-0.1 * camera->GetRadius());
    break;

    case GLUT_KEY_PAGE_UP:
      break;

    case GLUT_KEY_PAGE_DOWN:
      break;

    case GLUT_KEY_HOME:
      break;

    case GLUT_KEY_END:
      break;

    case GLUT_KEY_INSERT:
      break;

    default:
      break;
  }
}

void mouseMotion (int x, int y)
{
  g_vMousePos[0] = x;
  g_vMousePos[1] = y;
}

void mouseButtonActivityFunction(int button, int state, int x, int y)
{
  switch (button)
  {
    case GLUT_LEFT_BUTTON:
      g_iLeftMouseButton = (state==GLUT_DOWN);
      shiftPressed = (glutGetModifiers() == GLUT_ACTIVE_SHIFT);
      altPressed = (glutGetModifiers() == GLUT_ACTIVE_ALT);
      ctrlPressed = (glutGetModifiers() == GLUT_ACTIVE_CTRL);
      if (g_iLeftMouseButton)
      {
        GLdouble model[16];
        glGetDoublev (GL_MODELVIEW_MATRIX, model);
        
        GLdouble proj[16];
        glGetDoublev (GL_PROJECTION_MATRIX, proj);
        
        GLint view[4];
        glGetIntegerv (GL_VIEWPORT, view);
        
        int winX = x;
        int winY = view[3]-1-y;
        
        float zValue;
        glReadPixels(winX,winY,1,1, GL_DEPTH_COMPONENT, GL_FLOAT, &zValue); 
        
        GLubyte stencilValue;
        glReadPixels(winX, winY, 1, 1, GL_STENCIL_INDEX, GL_UNSIGNED_BYTE, &stencilValue);
        
        GLdouble worldX, worldY, worldZ;
        gluUnProject (winX, winY, zValue, model, proj, view, &worldX, &worldY, &worldZ);
        
        //printf("x:%d y:%d zValue:%f stencil:%d world: %f %f %f\n",
        //       winX,winY,zValue,stencilValue,worldX,worldY,worldZ);
        
        if(lockScene == 0)
        {
          if (stencilValue == 1)
          {
            dragStartX = x;
            dragStartY = y;
            //pulledVertex = renderCloth->findClosestVertex(clothBW, u, worldX, worldY, worldZ);
            Vec3d pos(worldX, worldY, worldZ);
            pulledVertex = sceneObjDeform->GetClosestVertex(pos, u);
            printf("Clicked on vertex: %d (0-indexed)\n", pulledVertex);
          }
          else
          {
            printf("Clicked on empty stencil: %d.\n", stencilValue);
            pulledVertex = -1;
          }
        }
        if(pin == 1)
        {
          if (stencilValue == 1)
          {
            //int vertex = renderCloth->findClosestVertex(clothBW, u, worldX, worldY, worldZ);
            Vec3d pos(worldX, worldY, worldZ);
            int vertex = sceneObjDeform->GetClosestVertex(pos, u);
            printf("Pinned on vertex: %d (0-indexed)\n", vertex);
            pin_points.push_back(vertex);
          }
        }
        if (!g_iLeftMouseButton)
        {
          pulledVertex = -1;
        }
      }
      if (!g_iLeftMouseButton)
        pulledVertex = -1;
      break;
    case GLUT_MIDDLE_BUTTON:
      g_iMiddleMouseButton = (state==GLUT_DOWN);
      break;
    case GLUT_RIGHT_BUTTON:
      g_iRightMouseButton = (state==GLUT_DOWN);
      break;
  }
  
  g_vMousePos[0] = x;
  g_vMousePos[1] = y;	
}

void mouseMotionFunction(int x, int y)
{
  int mouseDeltaX = x-g_vMousePos[0];
  int mouseDeltaY = y-g_vMousePos[1];
  
  g_vMousePos[0] = x;
  g_vMousePos[1] = y;
  
  if (g_iRightMouseButton) // handle camera rotations
  {
    const double factor = 0.1;
    camera->MoveRight(factor * mouseDeltaX);
    camera->MoveUp(factor * mouseDeltaY);
  }	
  
  if ((g_iMiddleMouseButton) || (g_iLeftMouseButton && altPressed)) // handle zoom in/out
  {
    const double factor = 0.1;
    camera->ZoomIn(cameraRadius * factor * mouseDeltaY);
  }
}

// this function is called when the user
// changes the timestep parameter from 
// within the simulator
void update_timestep(int code)
{
  integratorBase->SetTimestep(timeStep);
  glui->sync_live();
}

// this function is called when a button in the
// simulator is pressed
void update_simulation_status(int code)
{
  switch (code)
  {
    case 0: // initialize new parameters (need to re-initialize scene)
      initScene();
      //runSimulation = 0;
      break;
      
    case 1: // run
      runSimulation = 1;
      break;
      
    case 2: // pause
      runSimulation = 0;
      break;
      
    case 3: // reset cloth
      integratorBase->ResetToRest();
      memset(u, 0, sizeof(double) * 3 * n);
      sceneObjDeform->ResetDeformationToRest();
      //runSimulation = 0;
      break;
      
    case 4: // load new config file
      initConfigurations();
      initScene();
      initGraphics(windowWidth, windowHeight);
      //runSimulation = 0;
      break;
      
    default:
      printf("Error. Invalid simulator update option.\n");
      break;
  }
  
  glui->sync_live();
}

// this function does the following:
// (1) Creates the camera
// (2) Loads a .obj file
// (3) Creates a clothBW from the .obj file
// (4) Loads file containing constrained vertices
// (5) Creates an object to render the clothBW
// (6) Creates a lighting object to light the scene
// (7) Takes care of internal threading (if applicable)
// (8) Creates numerical integrator
// (9) Initializes textures on the cloth rendering object
void initScene()
{
  if(camera != NULL)	
    delete(camera);
  
  double virtualToPhysicalPositionFactor = 1.0;
  int numConstrainedDOFs = 0; 
  int * constrainedDOFs = NULL;
  
  initCamera(cameraRadius, cameraLongitude, cameraLatitude,
             focusPositionX, focusPositionY, focusPositionZ,
             1.0 / virtualToPhysicalPositionFactor,
             &zNear, &zFar, &camera);
  
  printf("Loading the OBJ mesh.\n");    
  
  if (objMesh != NULL)
    delete objMesh;
  objMesh = new ObjMesh(objMeshname);
  
  if (clothBWFromObjMesh != NULL)
    delete clothBWFromObjMesh;
  clothBWFromObjMesh = new ClothBWFromObjMesh();

  // the below region tests the assignment of different stiffness values 
  // to different material groups within the same cloth model
  /*
  // TESTING PURPOSES ONLY
  
  int numMaterialGroups = objMesh->numMaterials();
  
  double * surfaceDensityArray = (double*) malloc (sizeof(double) * numMaterialGroups);
  double * tensileStiffnessArray = (double*) malloc (sizeof(double) * numMaterialGroups);
  double * shearStiffnessArray = (double*) malloc (sizeof(double) * numMaterialGroups);
  double * bendStiffnessUArray = (double*) malloc (sizeof(double) * numMaterialGroups);
  double * bendStiffnessVArray = (double*) malloc (sizeof(double) * numMaterialGroups);
  double * dampingArray = (double*) malloc (sizeof(double) * numMaterialGroups);
  int * buArray = (int*) malloc (sizeof(int) * numMaterialGroups);
  int * bvArray = (int*) malloc (sizeof(int) * numMaterialGroups);
  
  for (int i=0; i<numMaterialGroups; i++)
  {
    surfaceDensityArray[i] = surfaceDensity;
    tensileStiffnessArray[i] = tensileStiffness;
    shearStiffnessArray[i] = shearStiffness;
    bendStiffnessUArray[i] = bendStiffnessU;
    bendStiffnessVArray[i] = bendStiffnessV;
    dampingArray[i] = 0.0;
    buArray[i] = 1;
    bvArray[i] = 1;
  }
  
  tensileStiffnessArray[2] = 99999.0f;
  shearStiffnessArray[2] =   99999.0f;
  
  clothBWFromObjMesh->GenerateClothBW(objMesh, &clothBW, surfaceDensityArray, tensileStiffnessArray,
                                  shearStiffnessArray, bendStiffnessUArray,
                                  bendStiffnessVArray, dampingArray, buArray, bvArray, addGravity);
  
  // free memory we allocated
  free(surfaceDensityArray);
  free(tensileStiffnessArray);
  free(shearStiffnessArray);
  free(bendStiffnessUArray);
  free(bendStiffnessVArray);
  free(dampingArray);
  free(buArray);
  free(bvArray);
  
  // END TESTING PURPOSES ONLY
  */

  if (clothBW != NULL)
    delete clothBW;
  ClothBW::MaterialGroup material;
  material.tensileStiffness = tensileStiffness;
  material.shearStiffness = shearStiffness;
  material.bendStiffnessU = bendStiffnessU;
  material.bendStiffnessV = bendStiffnessV;
  clothBW = clothBWFromObjMesh->GenerateClothBW(objMesh, surfaceDensity, material, addGravity);
  clothBW->SetGravity(addGravity, gravityForce);
  
  // constrain vertices from file
  LoadList loadlist;
  int offset = 1;
  if (loadlist.load(fixedVerticesFilename, &numFixedVertices, &fixedVertices, offset) != 0)
  {
    printf("Error reading fixed vertices.\n");
    exit(1);
  }
  //loadlist.print(numFixedVertices,fixedVertices);
  
  numConstrainedDOFs = 3*numFixedVertices;
  if (constrainedDOFs)
    free(constrainedDOFs);
  constrainedDOFs = (int*) malloc (sizeof(int) * numConstrainedDOFs);
  for (int i=0; i<numFixedVertices; i++)
  {
    constrainedDOFs[3*i+0] = fixedVertices[i] * 3 + 0;
    constrainedDOFs[3*i+1] = fixedVertices[i] * 3 + 1;
    constrainedDOFs[3*i+2] = fixedVertices[i] * 3 + 2;
  }
  
  // now create renderer
  if (sceneObjDeform != NULL)
    delete sceneObjDeform;
  sceneObjDeform = new SceneObjectDeformable(objMeshname);
  
  if (light != NULL)
    delete light;
  light = new Lighting(lightingFilename);
  
  // handle multi-threading
  if (numInternalForceThreads > 0)
  {
    printf("Launching threaded force/stiffness matrix evaluation: %d internal threads, %d solver threads.\n", numInternalForceThreads, numSolverThreads);
    
    // create multi-threaded version of clothBW
    ClothBWMT * clothBWMT = new ClothBWMT(*clothBW, numInternalForceThreads);
    
    // delete clothBW
    delete(clothBW);
    
    // and have old pointer now point to multi-threaded version
    clothBW = clothBWMT;
  }
  
  // initialize the integrator (ImplicitBackwardEulerSparse)
  if (clothBWForceModel != NULL)
    delete clothBWForceModel;
  clothBWForceModel = new ClothBWForceModel(clothBW);	
  forceModel = clothBWForceModel;
  SparseMatrix * massMatrix;
  clothBW->GenerateMassMatrix(&massMatrix);
  double totalMass = massMatrix->SumEntries() / 3.0;
  printf("Total cloth mass: %G\n", totalMass);
  
  int numVertices = clothBW->GetNumVertices();
  
  if (implicitBackwardEulerSparse != NULL)
    delete implicitBackwardEulerSparse;
  implicitBackwardEulerSparse = new ImplicitBackwardEulerSparse(3 * numVertices, timeStep, massMatrix, forceModel, numConstrainedDOFs, constrainedDOFs, dampingMass, dampingStiffness, 1, 1E-5, numSolverThreads);
  implicitNewmarkSparse = implicitBackwardEulerSparse;
  integratorBase = implicitNewmarkSparse;
  
  // initializing deformation
  n = numVertices;
  
  if (u)
    free(u);
  if (f_ext)
    free (f_ext);
  u = (double*) malloc (sizeof(double) * 3 * numVertices);
  f_ext = (double*) malloc (sizeof(double) * 3 * numVertices);
  for(int i = 0 ; i < 3 * numVertices; i++)
  {
    u[i] = 0.0;
    f_ext[i] = 0.0;
  }
  integratorBase->SetState(u);
  integratorBase->SetExternalForces(f_ext);
  
  pin_points.clear();

  sceneObjDeform->SetUpTextures(SceneObject::MODULATE);
  sceneObjDeform->BuildNeighboringStructure();

  // load any external geometry file (e.g. some static scene for decoration; usually there will be none)
  if (strcmp(extraSceneGeometryFilename,"__none") != 0)
  {
    extraSceneGeometry = new SceneObject(extraSceneGeometryFilename);
    extraSceneGeometry->BuildNormals(85.0);
  }
  else
    extraSceneGeometry = NULL;

  runSimulation = 1;
}

// Create the GUI 
void initGLUI()
{
  // generate the UI
  glui = GLUI_Master.create_glui( "Controls",0,windowWidth + 50, 0);

  // scene option
  GLUI_Panel * scene_option_panel = glui->add_panel("Scene&Parameters",GLUI_PANEL_EMBOSSED);
  scene_option_panel->set_alignment(GLUI_ALIGN_LEFT);

  glui->add_checkbox_to_panel(scene_option_panel,"WireFrame", &renderWireframe);
  
  //GLUI_RadioGroup * loadoption = glui->add_radiogroup_to_panel(scene_option_panel, &loadOption);
  
  GLUI_EditText * config_filename_text = glui->add_edittext_to_panel(scene_option_panel,"Config File Name",GLUI_EDITTEXT_TEXT,configFilename);
  config_filename_text->w = config_filename_text->w + 15;
  
  glui->add_button_to_panel(scene_option_panel,"Load New Config File",4, update_simulation_status);
  
  glui->add_separator_to_panel(scene_option_panel);
  glui->add_edittext_to_panel(scene_option_panel,"tensileStiffness",GLUI_EDITTEXT_FLOAT,&tensileStiffness);
  glui->add_edittext_to_panel(scene_option_panel,"shearStiffness",GLUI_EDITTEXT_FLOAT,&shearStiffness);
  glui->add_edittext_to_panel(scene_option_panel,"bendStiffnessU",GLUI_EDITTEXT_FLOAT,&bendStiffnessU);
  glui->add_edittext_to_panel(scene_option_panel,"bendStiffnessV",GLUI_EDITTEXT_FLOAT,&bendStiffnessV);
  glui->add_edittext_to_panel(scene_option_panel,"dampingMass",GLUI_EDITTEXT_FLOAT,&dampingMass);
  glui->add_edittext_to_panel(scene_option_panel,"dampingStiffness",GLUI_EDITTEXT_FLOAT,&dampingStiffness);
  glui->add_edittext_to_panel(scene_option_panel,"gravityForce",GLUI_EDITTEXT_FLOAT,&gravityForce);
  glui->add_edittext_to_panel(scene_option_panel,"mousePower",GLUI_EDITTEXT_FLOAT,&mouseForceCoeff);
  glui->add_checkbox_to_panel(scene_option_panel,"useRestAngles", &useRestAngles);
  glui->add_button_to_panel(scene_option_panel,"Upload Parameters",0, update_simulation_status);
  
  glui->add_separator();
  
  // external force: wind
  GLUI_Panel * extForce_panel = glui->add_panel("ExternalForce",GLUI_PANEL_EMBOSSED);
  extForce_panel->set_alignment(GLUI_ALIGN_LEFT);
  GLUI_Spinner * wind_x_spinner = glui->add_spinner_to_panel(extForce_panel, "WindX", GLUI_SPINNER_FLOAT, &wind_x);
  GLUI_Spinner * wind_y_spinner = glui->add_spinner_to_panel(extForce_panel, "WindY", GLUI_SPINNER_FLOAT, &wind_y);
  GLUI_Spinner * wind_z_spinner = glui->add_spinner_to_panel(extForce_panel, "WindZ", GLUI_SPINNER_FLOAT, &wind_z);
  wind_x_spinner->set_float_limits(-10.0, 10.0);
  wind_y_spinner->set_float_limits(-10.0, 10.0);
  wind_z_spinner->set_float_limits(-10.0, 10.0);
  wind_x_spinner->set_speed(0.1);
  wind_y_spinner->set_speed(0.1);
  wind_z_spinner->set_speed(0.1);
  
  // simulation control 
  GLUI_Panel * simulation_panel = glui->add_panel("Simulation",GLUI_PANEL_EMBOSSED);
  simulation_panel->set_alignment(GLUI_ALIGN_LEFT);
  simulation_panel->draw_name(0,0);
  glui->add_button_to_panel(simulation_panel,"Restart",3, update_simulation_status);
  glui->add_button_to_panel(simulation_panel,"Run",1, update_simulation_status);
  glui->add_button_to_panel(simulation_panel,"Pause",2, update_simulation_status);
  glui->add_edittext_to_panel(simulation_panel,"Timestep [sec]",GLUI_EDITTEXT_FLOAT,&timeStep,0,update_timestep);

  systemSolveStaticText = glui->add_statictext("System solve: ");
  forceAssemblyStaticText = glui->add_statictext("Force assembly: ");

  glui->add_separator();
  
  // end of the control
  glui->add_button("Quit", 0, exit_buttonCallBack);
  glui->sync_live();
  glui->set_main_gfx_window( windowID );
}

int main(int argc, char* argv[])
{
  int numFixedArgs = 2;
  
  if ( argc != numFixedArgs )
  {
    printf("=== Cloth Simulator ===\n");
    printf("Usage: %s [config file]\n", argv[0]);
    printf("Please specify a configuration file\n");
    return 1;
  }
  else
  {
    strncpy(configFilename, argv[1], strlen(argv[1]));
  }
  
  // make window and size it properly
  initGLUT(argc, argv, windowTitleBase, windowWidth, windowHeight, &windowID);
  
  // define background texture, set some openGL parameters
  initGraphics(windowWidth, windowHeight);

  // load info from config file
  initConfigurations();
  
  initScene();
  
  initGLUI();
  
  glutMainLoop(); 
  return 0;
}

