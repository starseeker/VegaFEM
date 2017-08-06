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

#include "uniqueIntegerID.h"
using namespace std;

UniqueIntegerID::UniqueIntegerID(unsigned int startID) 
{
  maxID = startID;
}

unsigned int UniqueIntegerID::Get()
{
  if (deletedIDs.size() > 0)
  {
    set<unsigned int> :: iterator iter = deletedIDs.begin();
    unsigned int value = *iter;
    deletedIDs.erase(iter);
    activeIDs.insert(value); 
    return value;
  }
  else
  {
    maxID++;
    activeIDs.insert(maxID - 1);
    return maxID - 1;
  }
}

int UniqueIntegerID::Release(unsigned int ID)
{
  set<unsigned int> :: iterator iter = activeIDs.find(ID);
  if (iter == activeIDs.end())
    return 1;

  unsigned int value = *iter;
  deletedIDs.insert(value);
  activeIDs.erase(iter);

  return 0;
}

void UniqueIntegerID::GetIDs(set<unsigned int> & IDs)
{
  IDs = activeIDs;
}

void UniqueIntegerID::Register(unsigned int ID)
{
  // if ID is already active, do nothing
  if (activeIDs.find(ID) != activeIDs.end())
    return;

  // if ID is in the deleted set, remove it from deleted set, and add it to active set
  if (deletedIDs.find(ID) != deletedIDs.end())
  {
    deletedIDs.erase(ID);
    activeIDs.insert(ID);
    return;
  }

  // now, we know ID >= maxID (unless the user used startID and then passed a small value as ID, which can be handled by ignoring it)
  // put ID into activeIDs
  // move the range from maxID to ID-1 into deleteIDs
  // set maxID to ID + 1
  if (ID >= maxID)
  {
    activeIDs.insert(ID);
    for(unsigned int i=maxID; i < ID; i++)
      deletedIDs.insert(i);
    maxID = ID + 1;
  }
}

void UniqueIntegerID::Clear(unsigned int startID)
{
  activeIDs.clear();
  deletedIDs.clear();
  maxID = startID;
}

