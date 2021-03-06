/*===================================================================================

  File: debug.hh

  ===================================================================================

    Copyright (C) 2016 Vicente Palazón-González

    This program is free software; you can redistribute it and/or modify it under the
    terms of the GNU General Public License as published by the Free Software
    Foundation; either version 3 of the License, or (at your option) any later
    version.

    This program is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
    PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with this
    program (file LICENSE.txt); if not, write to the Free Software Foundation, Inc.,
    675 Mass Ave, Cambridge, MA 02139, USA.

    You can contact the author at:

    palazon@uji.es
    
  =================================================================================== */

#ifndef _DEBUG_HH
#define _DEBUG_HH

#include <vector>
#include <list>
#include <set>

using namespace std;

void showVector(vector<float> &v);

void showVector(vector<int> &v);

void showVectorSpaced(vector<float> &v);

void showMatrix(vector< vector<int> > &m);

void showMatrix(vector< vector<float> > &m);

void showList(list<float> &l);

void showList(list<int> &l);

void showSet(set<int> &l);

#endif

