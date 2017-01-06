/*===================================================================================

  File: Bubu.cc

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

#include "Bubu.hh"

void succs(coord &u, unsigned limx, unsigned limy, vector<coord> &s) {
    s.clear();
    unsigned ui = u.first;
    unsigned uj = u.second;
    if (ui < limx-1 and uj < limy-1) {
        s.push_back(coord(ui+1, uj+1));
        s.push_back(coord(ui+1, uj));
        s.push_back(coord(ui, uj+1));
    }
    else if (ui == limx-1)
        s.push_back(coord(ui, uj+1));
    else if (uj == limy-1)
        s.push_back(coord(ui+1, uj));
}

// coordinates to key
unsigned coord2key(unsigned i, unsigned j, unsigned m) {
    return j*m+i;
}

// key to coordinates
void key2coord(unsigned k, unsigned m, 
               unsigned &i, unsigned &j) {
    j = k/m;
    i = k-(j*m);
}

unsigned setContains(set<coord> &s, coord &item) {
    return s.find(item) != s.end();
}

