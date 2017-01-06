/*===================================================================================

  File: Ed.hh

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

#ifndef _ED_HH
#define _ED_HH

#include "Samples.hh"
#include <vector>
#include <list>
#include <utility>
#include <cmath>
#include "matrix.h"

float ed(Sample<t_value> &sx, Sample<t_value> &sy, float (*d) (t_value, t_value));
float edExt(Sample<t_value> &sx, Sample<t_value> &sy, float (*d) (t_value, t_value), float best_so_far_distance, unsigned long long &ldistances_count);
float edExt2(Sample<t_value> &sx, Sample<t_value> &sy, float (*d) (t_value, t_value), float best_so_far_distance, unsigned long long &ldistances_count);
float edAndPath(Sample<t_value> &sx, Sample<t_value> &sy, float (*d) (t_value, t_value), Sample<t_pvalue> &path);
float cedBruteForce(Sample<t_value> &sx, Sample<t_value> &sy, float (*d) (t_value, t_value));
float maesEd(Sample<t_value> &sx, Sample<t_value> &sy, float (*distance) (t_value, t_value));
float maesEdExt1bubu(Sample<t_value> &sx, Sample<t_value> &sy, float (*distance) (t_value, t_value), float best_so_far_distance, unsigned long long &ldistances_count);
int maesEdShift(Sample<t_value> &sx, Sample<t_value> &sy, float (*distance) (t_value, t_value));
float BBEd(Sample<t_value> &sx, Sample<t_value> &sy, float (*distance) (t_value, t_value));
float BBEdExt(Sample<t_value> &sx, Sample<t_value> &sy, float (*distance) (t_value, t_value), float best_so_far_distance, unsigned long long &ldistances_count);
float BBEdExt1bubu(Sample<t_value> &sx, Sample<t_value> &sy, float (*distance) (t_value, t_value), float best_so_far_distance, unsigned long long &ldistances_count);

#endif
