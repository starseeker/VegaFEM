/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 3.0                               *
 *                                                                       *
 * "finite difference tester" utility , Copyright (C) 2016 USC           *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Yijing Li                                                *
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

#ifndef VECTORHELPER_H
#define VECTORHELPER_H

#include <vector>
#include <cassert>
#include <cstring>
#include <iostream>
#include <cmath>
#include <cassert>

template<class T>
inline void memset(std::vector<T> & v) 
{
  memset(v.data(), 0, sizeof(T) * v.size());
}

template<class T>
inline void memcpy(std::vector<T> & v, const T * src) 
{
  assert(src);
  memcpy(v.data(), src, sizeof(T) * v.size());
}

// free any memory allocated in the vector
template<class T>
inline void reset(std::vector<T> & v) 
{
  std::vector<T> tmp;
  v.swap(tmp);
}

/////////////////////////////////////////////////////
//                 Norm Computation                //
/////////////////////////////////////////////////////

// the Euclidean norm of the vector to the origin

inline double squaredEuclideanNorm(const std::vector<double> & v) 
{
  double sum = 0.;
  for(size_t i = 0; i < v.size(); i++)
    sum += v[i] * v[i];
  return sum;
}

inline double EuclideanNorm(const std::vector<double> & v) 
{
  return sqrt(squaredEuclideanNorm(v));
}

// the Euclidean distance between the two vector of the same size
inline double squaredEuclideanDistance(const std::vector<double> & v1, const std::vector<double> & v2) 
{
  double sum = 0;
  assert(v1.size() == v2.size());
  for(size_t i = 0; i < v1.size(); i++) 
  {
    double value = (v1[i] - v2[i]);
    sum += value * value;
  }
  return sum;
}

inline double EuclideanDistance(const std::vector<double> & v1, const std::vector<double> & v2) 
{
  return sqrt(squaredEuclideanDistance(v1,v2));
}

inline double EuclideanNorm(size_t r, const double * v) 
{
  double sum = 0;
  for(size_t i = 0; i < r; i++)
    sum += v[i] * v[i];
  return sqrt(sum);
}

/////////////////////////////////////////////////////
//         SAVE & LOAD IMPLEMENTATIONS             //
/////////////////////////////////////////////////////
// for use in saving/loading a class which includes vector as member var.
// Same save/load file format as used in vega mayaPlugin
template<class T>
bool saveToAscii(const std::vector<T> & v, std::ostream & out);

template<class T>
bool loadFromAscii(std::vector<T> & v, std::istream & in);

#endif

