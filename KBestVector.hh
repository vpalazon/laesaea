/*===================================================================================

  File: KBestVector.hh

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

#ifndef _KBESTVECTOR_HH
#define _KBESTVECTOR_HH

#include <vector>
#include <utility>

using namespace std;

typedef pair <unsigned, float> i_cost;

class KBestVector {
private:
    vector<i_cost> v;
public:
    KBestVector(unsigned n);
    float top(void);
    void insert(unsigned i, float cost);
    void insertExists(unsigned i, float cost);
    void get_sort(vector<int> &knn);
    void get_sort(vector< pair<unsigned, float> > &knn);
};

#endif

