/*===================================================================================

  File: Elem.cc

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

#include "Elem.hh"

bool fMaxHeap(Elem *x, Elem *y) {
    if (x->cost == y->cost)
        return x->e < y->e;
    return x->cost < y->cost;
}

bool fMinHeap(Elem *x, Elem *y) {
    if (x->cost == y->cost)
        return x->e > y->e;
    return x->cost > y->cost;
}

bool compve(int i1, int i2) {
    return i1>i2;
}

bool comp_Mm(Elem *x, Elem *y) {
    return x->cost > y->cost;
}

bool comp_mM(Elem *x, Elem *y) {
    return x->cost < y->cost;
}

