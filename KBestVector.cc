/*===================================================================================

  File: KBestVector.cc

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

#include "KBestVector.hh"
#include "Defs.hh"

#include <algorithm>

KBestVector::KBestVector(unsigned n) {
    for (unsigned i=0; i<n+1; i++)
        v.push_back(i_cost(-1, INF));
}

float KBestVector::top(void) {
    return (*v.begin()).second;
}

void KBestVector::insert(unsigned i, float cost) {
    v[0].first = i;
    v[0].second = cost;

    for (unsigned k=1; k<v.size(); ++k) {
        if (v[k-1].second == v[k].second) {
            if (v[k-1].first < v[k].first) {
                swap(v[k].first, v[k-1].first);
                swap(v[k].second, v[k-1].second);
            }
        }
        else if (v[k-1].second < v[k].second) {
            swap(v[k].first, v[k-1].first);
            swap(v[k].second, v[k-1].second);
        }
    }
}

void KBestVector::insertExists(unsigned i, float cost) {

    unsigned found = 0;
    
    unsigned k;
    for (k=0; k<v.size(); ++k)
        if (v[k].first == i) {
            found = 1;
            break;
        }
    
    unsigned k_start;
    if (found)
        k_start = k;
    else
        k_start = 0;

    v[k_start].first = i;
    v[k_start].second = cost;

    for (unsigned i=k_start+1; i<v.size(); ++i) {
        if (v[i-1].second == v[i].second) {
            if (v[i-1].first < v[i].first) {
                swap(v[i].first, v[i-1].first);
                swap(v[i].second, v[i-1].second);
            }
        }
        else if (v[i-1].second < v[i].second) {
            swap(v[i].first, v[i-1].first);
            swap(v[i].second, v[i-1].second);
        }
    }
}

void KBestVector::get_sort(vector<int> &knn) {
    unsigned i=0;
    for (vector<i_cost>::iterator it=--v.end(); it!=v.begin(); --it) {
        knn[i] = (*it).first;
        i++;
    }
}

void KBestVector::get_sort(vector< pair<unsigned, float> > &knn) {
    unsigned i=0;
    for (vector<i_cost>::iterator it=--v.end(); it!=v.begin(); --it) {
        knn[i].first = (*it).first;
        knn[i].second = (*it).second;
        i++;
    }
}
