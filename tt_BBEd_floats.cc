/*===================================================================================

  File: tt_BBEd_floats.cc

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

#include <iostream>
#include <string>
#include <cstdlib>
#include <set>

#include "Ed.hh"
#include "Samples.hh"
#include "LocalDistances.hh"
#include "Defs.hh"
#include "Useful.hh"

#include "debug.hh"

using namespace std;

int main(int argc, char **argv) {

    if (argc != 2) {
        cerr << "use: %s file.floats \n" << endl;
        exit(1);
    }

    float (*dist) (t_value, t_value) = levenshtein;

    Samples<t_value> seqs;
    seqs.read(argv[1]);
    
    for (unsigned i=0; i<seqs.samples->size(); ++i) {
        for (unsigned j=0; j<seqs.samples->size(); ++j) {
            float cost = BBEd((*seqs.samples)[i],
                             (*seqs.samples)[j], dist);
            cout << i << " " << j << " " << cost << endl;
        }
    }
    
    return 0;
}

