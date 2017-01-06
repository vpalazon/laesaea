/*===================================================================================

  File: Tt.hh

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

#ifndef _TT_HH
#define _TT_HH

#include <vector>
#include <string>
#include <fstream>
#include <strstream>
#include <iostream>
#include <cstdlib>
#include "matrix.h"
#include "Samples.hh"

using namespace std;

#define MAX_LINE_SZ_TT 200000
#define MAX_VALUE_SZ_TT 100

class Tt {
public:
    matrix<float> *distances;
    ~Tt();
    void read(string &filename);
    void read(string &filename, int size);
    unsigned size();
    float get(unsigned i, unsigned j);
};


#endif
