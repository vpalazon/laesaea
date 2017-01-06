/*===================================================================================

  File: Elem.hh

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

#ifndef _ELEM_HH
#define _ELEM_HH

#include <vector>

class Elem {
public:
    int e;
    float cost;

    Elem(int e, float cost):
        e(e), cost(cost) {}
};

bool fMaxHeap(Elem *x, Elem *y);
bool fMinHeap(Elem *x, Elem *y);

bool compve(int i1, int i2);

#endif
