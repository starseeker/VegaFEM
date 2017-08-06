/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 3.0                               *
 *                                                                       *
 * "sceneObject" library , Copyright (C) 2007 CMU, 2009 MIT, 2016 USC    *
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
   Generates unique integer IDs. Makes it possible to release IDs.
*/

#ifndef _UNIQUEINTEGERID_
#define _UNIQUEINTEGERID_

#include <set>

class UniqueIntegerID
{
public:
  UniqueIntegerID(unsigned int startID=0);

  // get a unique ID
  unsigned int Get();

  // registers an already existing ID
  void Register(unsigned int ID);

  // returns 0 on success, 1 if ID does not exist
  int Release(unsigned int ID);

  void GetIDs(std::set<unsigned int> & IDs);

  void Clear(unsigned int startID=0); // clears all IDs

protected:
  unsigned int maxID;
  std::set<unsigned int> activeIDs;
  std::set<unsigned int> deletedIDs;
};

#endif

