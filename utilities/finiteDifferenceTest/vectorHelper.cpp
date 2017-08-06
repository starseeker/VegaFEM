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

#include "vectorHelper.h"
using namespace std;

template<class T>
bool saveToAscii(const std::vector<T> & v, std::ostream & out) 
{
  out << v.size() << " ";
  if(out.fail()) 
    return false;
  for(size_t i = 0; i < v.size(); i++) 
  {
    out << v[i] << " ";
    if(out.fail()) 
      return false;
  }
  return true;
}

template<class T>
bool loadFromAscii(std::vector<T> & v, std::istream & in) 
{
  size_t num = 0;
  in >> num;
  if(in.fail())  
    return false;
  v.reserve(v.size() + num);
  for(size_t i = 0; i < num; i++) 
  {
    T k = T();
    in >> k;
    if(in.fail()) 
      return false;
    v.push_back(k);
  }
  return true;
}

template bool saveToAscii<int>(const std::vector<int> & v, std::ostream & out);
template bool saveToAscii<double>(const std::vector<double> & v, std::ostream & out);

template bool loadFromAscii<int>(std::vector<int> & v, std::istream & in);
template bool loadFromAscii<double>(std::vector<double> & v, std::istream & in);

