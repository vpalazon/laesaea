/*===================================================================================

  File: LocalDistances.cc

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

#include "LocalDistances.hh"
#include "Samples.hh"
#include <cmath>
#include "Defs.hh"

float levenshtein(t_value a, t_value b) {
    if ( a == LAMBDA and b == LAMBDA)
        return 0;
    else if ( a == LAMBDA or b == LAMBDA )
        return 1;
    else if ( a == b )
        return 0;
    else
        return 1;
}

