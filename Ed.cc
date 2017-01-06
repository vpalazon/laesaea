/*===================================================================================

  File: Ed.cc

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

#include "Defs.hh"
#include "Ed.hh"
#include <fstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <vector>
#include <utility>
#include <iterator>
#include <iostream>

#include "Samples.hh"
#include "matrix.h"

using namespace std;

float ed(Sample<t_value> &sx, Sample<t_value> &sy, float (*d) (t_value, t_value)) {
    vector<t_value> &x = sx.values;
    vector<t_value> &y = sy.values;

    int lenx = x.size();
    int leny = y.size();

    matrix<float> D = matrix<float>(lenx + 1, leny + 1);
    float d0, d1, d2;

    D[0][0] = 0;
    for (int i = 1; i < lenx+1; i++) D[i][0] = D[i-1][0] + d(x[i-1], LAMBDA);
    for (int j = 1; j < leny+1; j++) D[0][j] = D[0][j-1] + d(LAMBDA, y[j-1]);
    for (int i = 1; i < lenx+1; i++)
        for (int j = 1; j < leny+1; j++) {
            d0 = D[i-1][j] + d(x[i-1], LAMBDA);
            d1 = D[i][j-1] + d(LAMBDA, y[j-1]);
            d2 = D[i-1][j-1] + d(x[i-1], y[j-1]);
            if(d1 <= d0 and d1 <= d2)
                D[i][j] = d1;
            else if (d0 <= d2)
                D[i][j] = d0;
            else
                D[i][j] = d2;
        }
    return D[lenx][leny];
}

float edExt(Sample<t_value> &sx, Sample<t_value> &sy, float (*d) (t_value, t_value),
            float best_so_far_distance, unsigned long long &ldistances_count) {

    // just levenshtein {
    if (abs( (int)(sx.values.size() - sy.values.size()) ) >= best_so_far_distance)
        return INF;
    // }

    vector<t_value> &x = sx.values;
    vector<t_value> &y = sy.values;

    int lenx = x.size();
    int leny = y.size();

    matrix<float> D = matrix<float>(lenx + 1, leny + 1);
    float d0, d1, d2;

    unsigned abort;

    D[0][0] = 0;

    for (int i = 1; i < lenx+1; i++) {
        D[i][0] = D[i-1][0] + d(x[i-1], LAMBDA);
    }

    for (int j = 1; j < leny+1; j++) {
        abort = 1;

        D[0][j] = D[0][j-1] + d(LAMBDA, y[j-1]);
        if (D[0][j] < best_so_far_distance)
            abort = 0;

        for (int i = 1; i < lenx+1; i++) {
            d0 = D[i-1][j] + d(x[i-1], LAMBDA);
            d1 = D[i][j-1] + d(LAMBDA, y[j-1]);
            d2 = D[i-1][j-1] + d(x[i-1], y[j-1]);
            if(d1 <= d0 and d1 <= d2)
                D[i][j] = d1;
            else if (d0 <= d2)
                D[i][j] = d0;
            else
                D[i][j] = d2;

            if (D[i][j] < best_so_far_distance)
                abort = 0;
        }

        if (abort)
            return INF;
    }

    return D[lenx][leny];
}

float edExt2(Sample<t_value> &sx, Sample<t_value> &sy, float (*d) (t_value, t_value),
             float best_so_far_distance, unsigned long long &ldistances_count) {


    vector<t_value> &x = sx.values;
    vector<t_value> &y = sy.values;

    int lenx = x.size();
    int leny = y.size();

    matrix<float> D = matrix<float>(lenx + 1, leny + 1);
    float d0, d1, d2;

    D[0][0] = 0;

    for (int i = 1; i < lenx+1; i++) {
        D[i][0] = D[i-1][0] + d(x[i-1], LAMBDA);
    }

    int first_column = 1;
    int column = 1;

    for (int j = 1; j < leny+1; j++) {
        int found = 0;
        
        if (first_column) {
            D[0][j] = D[0][j-1] + d(LAMBDA, y[j-1]);
            if (D[0][j] >= best_so_far_distance) {
                first_column = 0;
            }
            for (int i = 1; i < lenx+1; i++) {
                d0 = D[i-1][j] + d(x[i-1], LAMBDA);
                d1 = D[i][j-1] + d(LAMBDA, y[j-1]);
                d2 = D[i-1][j-1] + d(x[i-1], y[j-1]);
                if(d1 <= d0 and d1 <= d2)
                    D[i][j] = d1;
                else if (d0 <= d2)
                    D[i][j] = d0;
                else
                    D[i][j] = d2;
            }
        }
        else {
            int i = column;
            while ((not found) && i<lenx+1) {
                d0 = D[i-1][j] + d(x[i-1], LAMBDA);
                d1 = D[i][j-1] + d(LAMBDA, y[j-1]);
                d2 = D[i-1][j-1] + d(x[i-1], y[j-1]);
                if (d1 <= d0 and d1 <= d2)
                    D[i][j] = d1;
                else if (d0 <= d2)
                    D[i][j] = d0;
                else
                    D[i][j] = d2;
                
                if (D[i][j] < best_so_far_distance) {
                    found = 1;
                    column = i;
                }
                i++;
            }

            if (i == lenx)
                return INF;

            for (i = column+1; i < lenx+1; i++) {
                d0 = D[i-1][j] + d(x[i-1], LAMBDA);
                d1 = D[i][j-1] + d(LAMBDA, y[j-1]);
                d2 = D[i-1][j-1] + d(x[i-1], y[j-1]);
                if(d1 <= d0 and d1 <= d2)
                    D[i][j] = d1;
                else if (d0 <= d2)
                    D[i][j] = d0;
                else
                    D[i][j] = d2;
            }
            
        }

    }

    return D[lenx][leny];
}

float edAndPath(Sample<t_value> &sx, Sample<t_value> &sy, float (*d) (t_value, t_value), Sample<t_pvalue> &path) {
    vector<t_value> &x = sx.values;
    vector<t_value> &y = sy.values;

    int lenx = x.size();
    int leny = y.size();

    matrix<float> D = matrix<float>(lenx + 1, leny + 1);
    matrix< pair<int, int> > backpointer(lenx + 1, leny + 1);

    float d0, d1, d2;

    D[0][0] = 0;
    for (int i = 1; i < lenx+1; i++) {
        D[i][0] = D[i-1][0] + d(x[i-1], LAMBDA);
        backpointer[i][0].first = 1;
        backpointer[i][0].second = 0;
    }
    for (int j = 1; j < leny+1; j++) {
        D[0][j] = D[0][j-1] + d(LAMBDA, y[j-1]);
        backpointer[0][j].first = 0;
        backpointer[0][j].second = 1;
    }
    for (int i = 1; i < lenx+1; i++)
        for (int j = 1; j < leny+1; j++) {
            d0 = D[i-1][j] + d(x[i-1], LAMBDA);
            d1 = D[i][j-1] + d(LAMBDA, y[j-1]);
            d2 = D[i-1][j-1] + d(x[i-1], y[j-1]);
            if(d1 <= d0 and d1 <= d2) {
                D[i][j] = d1;
                backpointer[i][j].first = 0;
                backpointer[i][j].second = 1;
            }
            else if (d0 <= d2) {
                D[i][j] = d0;
                backpointer[i][j].first = 1;
                backpointer[i][j].second = 0;
            }
            else {
                D[i][j] = d2;
                backpointer[i][j].first = 1;
                backpointer[i][j].second = 1;
            }
        }

    unsigned i = lenx;
    unsigned j = leny;

    list<t_pvalue> l;
    pair<float, float> p;
    p.first = i;
    p.second = j;
    l.push_back(p);

    unsigned ri, rj;
    while (i != 0 or j != 0) {
        ri = backpointer[i][j].first;
        rj = backpointer[i][j].second;
        i -= ri;
        j -= rj;
        p.first = i;
        p.second = j;
        l.push_back(p);
    }

//     if (i == 0)
//         for (int k=j; k>=0; --k) {
//           p.first = 0;
//           p.second = k;
//           l.push_back(p);
//       }

//     if (j == 0)
//       for (int k=i; k>=0; --k) {
//           p.first = k;
//           p.second = 0;
//           l.push_back(p);
//       }

    l.reverse();

    for (list<t_pvalue>::iterator it=l.begin(); it!=l.end(); it++)
        path.values.push_back(*it);

    return D[lenx][leny];
}

float cedBruteForce(Sample<t_value> &sx, Sample<t_value> &sy, float (*d) (t_value, t_value)) {
    vector<t_value> &x = sx.values;
    vector<float> distances;
    int lenx = x.size();
    distances.reserve(lenx);
    for (int i=0; i<lenx; i++) {
        x.insert(x.begin(), *(x.end()-1));
        x.erase(x.end()-1);
        distances.push_back(ed(sx, sy, d));
    }
    return *min_element(distances.begin(), distances.end());
}

float xed(int left, int right, int s,
	  matrix<float> &D, matrix< pair<int, int> > &backpointer,
          const vector<t_value> &X,
	  const vector<t_value> &x, const vector<t_value> &y,
	  float (*distance) (t_value, t_value),
	  matrix<int> &mins, matrix<int> &maxs) {

    int lenx = x.size();
    int leny = y.size();

    vector<int> &lmin = mins[left], &lmax = maxs[left], &rmin = mins[right], &rmax = maxs[right];

    D[s][0] = 0;
    int inicio = max(lmin[0], s);
    D[inicio-1][0] = INFINITY;

    int sr = min(s+lenx, rmax[0])+1;

    for (int i=s+1; i < sr; i++) {
        D[i][0] = D[i-1][0] + distance(X[i-1], LAMBDA);
        backpointer[i][0].first = 1;
        backpointer[i][0].second = 0;
    }

    for (int j=1; j<leny+1; j++) {
        int ini = max(lmin[j], s), fin = min(s+lenx+1, rmax[j]+1);
        D[ini-1][j] = INFINITY;
        if ((ini > lmin[j] and lmax[j-1] == lmin[j]) or ini==s)
            D[ini-1][j-1] = INFINITY;
        if (fin > rmin[j]) {
            if (rmax[j-1] < rmin[j]) D[rmin[j]][j-1] = INFINITY;
            for (int i=rmin[j]+1; i<fin+1; i++)
                D[i][j-1] = INFINITY;
        }
        else if (fin == rmin[j] and rmax[j-1] < fin)
            D[fin][j-1] = INFINITY;
        for (int i= ini; i<fin; i++) {
            float d0 = D[i-1][j] + distance(X[i-1], LAMBDA);
            float d1 = D[i][j-1] + distance(LAMBDA, y[j-1]);
            float d2 = D[i-1][j-1] + distance(X[i-1], y[j-1]);
            if(d1 <= d0 and d1 <= d2) {
                D[i][j] = d1;
                backpointer[i][j].first = 0;
                backpointer[i][j].second = 1;
            }
            else if (d0 <= d2) {
                D[i][j] = d0;
                backpointer[i][j].first = 1;
                backpointer[i][j].second = 0;
            }
            else {
                D[i][j] = d2;
                backpointer[i][j].first = 1;
                backpointer[i][j].second = 1;
            }
        }
    }

    int i = s+lenx, j = leny;
    mins[s][j] = maxs[s][j] = i;
    int ri, rj;
    while (i != 0 and j != 0) {
        ri = backpointer[i][j].first;
        rj = backpointer[i][j].second;
        i -= ri;
        j -= rj;
        if (i < mins[s][j]) mins[s][j] = i;
        if (i > maxs[s][j]) maxs[s][j] = i;
    }
    for (int k=j; k>=0; --k) {
        mins[s][k] = s;
        if (i > maxs[s][k]) maxs[s][k] = i;
    }
    return D[s+lenx][leny];
}

float maesEdStep(int left, int right,
		 matrix<float> &D,
                 matrix< pair<int, int> > &backpointer,
		 const vector<t_value> &X,
		 const vector<t_value> &x, const vector<t_value> &y,
		 float (*distance) (t_value, t_value),
		 matrix<int> &mins, matrix<int> &maxs) {

    float r = INFINITY;
    int shift = (left+right)/2;
    float rr;

    if (left < shift and shift < right) {
        int s = (left+right)/2;
        rr = xed(left, right, s, D, backpointer, X, x, y, distance, mins, maxs);
        if (rr < r) r = rr;
    }

    if (left < shift-1) {
        rr = maesEdStep(left, shift, D, backpointer, X, x, y, distance, mins, maxs);
        if (rr < r) r = rr;
    }

    if (shift < right-1) {
        rr = maesEdStep(shift, right, D, backpointer, X, x, y, distance, mins, maxs);
        if (rr < r) r = rr;
    }

    return r;
}

float maesEd(Sample<t_value> &sx, Sample<t_value> &sy, float (*distance) (t_value, t_value)) {

    vector<t_value> &x = sx.values;
    vector<t_value> &y = sy.values;

    // if (sx.values.size() > sy.values.size()) { 
    //                                            
    //     vector<t_value> &swap = y;
    //     y = x;
    //     x = swap;
    // }

    int lenx = x.size();
    int leny = y.size();

    matrix<float> D(2*lenx+1, leny+1);
    matrix< pair<int, int> > backpointer(2*lenx+1, leny+1);
    D[0][0] = 0;
    for (int i = 1; i < lenx+1; i++) {
        D[i][0] = D[i-1][0] + distance(x[i-1], LAMBDA);
        backpointer[i][0].first = 1;
        backpointer[i][0].second = 0;
    }
    for (int j = 1; j < leny+1; j++) {
        D[0][j] = D[0][j-1] + distance(LAMBDA, y[j-1]);
        backpointer[0][j].first = 0;
        backpointer[0][j].second = 1;
    }
    for (int i = 1; i < lenx+1; i++)
        for (int j = 1; j < leny+1; j++) {
            float d0 = D[i-1][j] + distance(x[i-1], LAMBDA);
            float d1 = D[i][j-1] + distance(LAMBDA, y[j-1]);
            float d2 = D[i-1][j-1] + distance(x[i-1], y[j-1]);
            if(d1 <= d0 and d1 <= d2) {
                D[i][j] = d1;
                backpointer[i][j].first = 0;
                backpointer[i][j].second = 1;
            }
            else if (d0 <= d2) {
                D[i][j] = d0;
                backpointer[i][j].first = 1;
                backpointer[i][j].second = 0;
            }
            else {
                D[i][j] = d2;
                backpointer[i][j].first = 1;
                backpointer[i][j].second = 1;
            }
        }

    float edistance = D[lenx][leny];

    matrix<int> mins(lenx+2, leny+1);
    matrix<int> maxs(lenx+2, leny+1);
    for (int i = 0; i < lenx+2; i++)
        for (int j = 0; j < leny+1; j++) {
            mins[i][j] = 2*lenx+1;
            maxs[i][j] = 0;
        }
    int i = lenx, j = leny;
    mins[0][j] = maxs[0][j] = i;

    int ri, rj;
    while (i != 0 and j != 0) {
        ri = backpointer[i][j].first;
        rj = backpointer[i][j].second;
        i -= ri;
        j -= rj;
        if (i < mins[0][j]) mins[0][j] = i;
        if (i > maxs[0][j]) maxs[0][j] = i;
    }

    for (int k=j; k>=0; --k) {
        mins[0][k] = 0;
        if (i > maxs[0][k]) maxs[0][k] = i;
    }

    for (int k=0; k<leny+1; k++) {
        mins[lenx][k] = mins[0][k] + lenx;
        maxs[lenx][k] = maxs[0][k] + lenx;
    }

    // X = x+x
    vector<t_value> X(2*lenx);
    for (int i = 0; i < lenx; i++) X[i] = x[i];
    for (int i = 0; i < lenx; i++) X[i+lenx] = x[i];

    float r =  maesEdStep(0, lenx, D, backpointer, X, x, y, distance, mins, maxs);

    return min(edistance, r);
}

float maesEdExt1bubu(Sample<t_value> &sx, Sample<t_value> &sy, float (*distance) (t_value, t_value), float best_so_far_distance, unsigned long long &ldistances_count) {

    vector<t_value> x;
    vector<t_value> y;

    if (sx.values.size() > sy.values.size()) {
                                              
        x = sy.values;
        y = sx.values;
    }
    else {
        x = sx.values;
        y = sy.values;
    }
    
    int lenx = x.size();
    int leny = y.size();

    matrix<t_value> D(2*lenx+1, leny+1);

    // X = x+x
    vector<t_value> X(2*lenx);
    for (int i = 0; i < lenx; i++) X[i] = x[i];
    for (int i = 0; i < lenx; i++) X[i+lenx] = x[i];

    for (unsigned i = 0; i < lenx+1; i++) {
        D[i][0] = 0;
    }

    for (unsigned i = lenx+1; i < 2*lenx+1; i++) {
        D[i][0] = INF;
    }

    unsigned abort;

    for (unsigned j = 1; j < leny+1; j++) {

        abort = 1;

        D[0][j] = D[0][j-1] + distance(LAMBDA, y[j-1]);
        if (D[0][j] < best_so_far_distance)
            abort = 0;

        for (unsigned i = 1; i < 2*lenx+1; i++) {

            float d0 = D[i-1][j] + distance(X[i-1], LAMBDA);
            float d1 = D[i][j-1] + distance(LAMBDA, y[j-1]);
            float d2 = D[i-1][j-1] + distance(X[i-1], y[j-1]);


            if (d1 <= d0 and d1 <= d2) {
                D[i][j] = d1;
            }
            else if (d0 <= d2) {
                D[i][j] = d0;
            }
            else {
                D[i][j] = d2;
            }

            if (D[i][j] < best_so_far_distance)
                abort = 0;
        }

        if (abort) {
            return INF;
        }
    }

    matrix< pair<int, int> > backpointer(2*lenx+1, leny+1);
    D[0][0] = 0;
    for (int i = 1; i < lenx+1; i++) {
        D[i][0] = D[i-1][0] + distance(x[i-1], LAMBDA);
        backpointer[i][0].first = 1;
        backpointer[i][0].second = 0;
    }
    for (int j = 1; j < leny+1; j++) {
        D[0][j] = D[0][j-1] + distance(LAMBDA, y[j-1]);
        backpointer[0][j].first = 0;
        backpointer[0][j].second = 1;
    }
    for (int i = 1; i < lenx+1; i++)
        for (int j = 1; j < leny+1; j++) {
            float d0 = D[i-1][j] + distance(x[i-1], LAMBDA);
            float d1 = D[i][j-1] + distance(LAMBDA, y[j-1]);
            float d2 = D[i-1][j-1] + distance(x[i-1], y[j-1]);
            if(d1 <= d0 and d1 <= d2) {
                D[i][j] = d1;
                backpointer[i][j].first = 0;
                backpointer[i][j].second = 1;
            }
            else if (d0 <= d2) {
                D[i][j] = d0;
                backpointer[i][j].first = 1;
                backpointer[i][j].second = 0;
            }
            else {
                D[i][j] = d2;
                backpointer[i][j].first = 1;
                backpointer[i][j].second = 1;
            }
        }

    float edistance = D[lenx][leny];

    matrix<int> mins(lenx+2, leny+1);
    matrix<int> maxs(lenx+2, leny+1);
    for (int i = 0; i < lenx+2; i++)
        for (int j = 0; j < leny+1; j++) {
            mins[i][j] = 2*lenx+1;
            maxs[i][j] = 0;
        }
    int i = lenx, j = leny;
    mins[0][j] = maxs[0][j] = i;

    int ri, rj;
    while (i != 0 and j != 0) {
        ri = backpointer[i][j].first;
        rj = backpointer[i][j].second;
        i -= ri;
        j -= rj;
        if (i < mins[0][j]) mins[0][j] = i;
        if (i > maxs[0][j]) maxs[0][j] = i;
    }

    for (int k=j; k>=0; --k) {
        mins[0][k] = 0;
        if (i > maxs[0][k]) maxs[0][k] = i;
    }

    for (int k=0; k<leny+1; k++) {
        mins[lenx][k] = mins[0][k] + lenx;
        maxs[lenx][k] = maxs[0][k] + lenx;
    }

    float r =  maesEdStep(0, lenx, D, backpointer, X, x, y, distance, mins, maxs);

    return min(edistance, r);
}


// Shift {

float maesEdStepShift(int left, int right,
                      matrix<float> &D,
                      matrix< pair<int, int> > &backpointer,
                      const vector<t_value> &X,
                      const vector<t_value> &x, const vector<t_value> &y,
                      float (*distance) (t_value, t_value),
                      matrix<int> &mins, matrix<int> &maxs, int &rsh) {

    float r = INFINITY;
    int shift = (left+right)/2;
    float rr;
    int sh;

    sh = shift;

    if (left < shift and shift < right) {
        int s = (left+right)/2;
        rr = xed(left, right, s, D, backpointer, X, x, y, distance, mins, maxs);
        if (rr < r) {
            r = rr;
            rsh = sh;
        }
    }

    if (left < shift-1) {
        rr = maesEdStepShift(left, shift, D, backpointer, X, x, y, distance, mins, maxs, sh);
        if (rr < r) {
            r = rr;
            rsh = sh;
        }
    }

    if (shift < right-1) {
        rr = maesEdStepShift(shift, right, D, backpointer, X, x, y, distance, mins, maxs, sh);
        if (rr < r) {
            r = rr;
            rsh = sh;
        }
    }

    return r;
}

int maesEdShift(Sample<t_value> &sx, Sample<t_value> &sy, float (*distance) (t_value, t_value)) {

    vector<t_value> &x = sx.values;
    vector<t_value> &y = sy.values;
    int lenx = x.size();
    int leny = y.size();
    int rsh, sh;

    matrix<float> D(2*lenx+1, leny+1);
    matrix< pair<int, int> > backpointer(2*lenx+1, leny+1);
    D[0][0] = 0;
    for (int i = 1; i < lenx+1; i++) {
        D[i][0] = D[i-1][0] + distance(x[i-1], LAMBDA);
        backpointer[i][0].first = 1;
        backpointer[i][0].second = 0;
    }
    for (int j = 1; j < leny+1; j++) {
        D[0][j] = D[0][j-1] + distance(LAMBDA, y[j-1]);
        backpointer[0][j].first = 0;
        backpointer[0][j].second = 1;
    }
    for (int i = 1; i < lenx+1; i++)
        for (int j = 1; j < leny+1; j++) {
            float d0 = D[i-1][j] + distance(x[i-1], LAMBDA);
            float d1 = D[i][j-1] + distance(LAMBDA, y[j-1]);
            float d2 = D[i-1][j-1] + distance(x[i-1], y[j-1]);
            if(d1 <= d0 and d1 <= d2) {
                D[i][j] = d1;
                backpointer[i][j].first = 0;
                backpointer[i][j].second = 1;
            }
            else if (d0 <= d2) {
                D[i][j] = d0;
                backpointer[i][j].first = 1;
                backpointer[i][j].second = 0;
            }
            else {
                D[i][j] = d2;
                backpointer[i][j].first = 1;
                backpointer[i][j].second = 1;
            }
        }

    float edistance = D[lenx][leny];

    matrix<int> mins(lenx+2, leny+1);
    matrix<int> maxs(lenx+2, leny+1);
    for (int i = 0; i < lenx+2; i++)
        for (int j = 0; j < leny+1; j++) {
            mins[i][j] = 2*lenx+1;
            maxs[i][j] = 0;
        }
    int i = lenx, j = leny;
    mins[0][j] = maxs[0][j] = i;

    int ri, rj;
    while (i != 0 and j != 0) {
        ri = backpointer[i][j].first;
        rj = backpointer[i][j].second;
        i -= ri;
        j -= rj;
        if (i < mins[0][j]) mins[0][j] = i;
        if (i > maxs[0][j]) maxs[0][j] = i;
    }

    for (int k=j; k>=0; --k) {
        mins[0][k] = 0;
        if (i > maxs[0][k]) maxs[0][k] = i;
    }

    for (int k=0; k<leny+1; k++) {
        mins[lenx][k] = mins[0][k] + lenx;
        maxs[lenx][k] = maxs[0][k] + lenx;
    }

    // X = x+x
    vector<t_value> X(2*lenx);
    for (int i = 0; i < lenx; i++) X[i] = x[i];
    for (int i = 0; i < lenx; i++) X[i+lenx] = x[i];

    float r =  maesEdStepShift(0, lenx, D, backpointer, X, x, y, distance, mins, maxs, sh);

    if (edistance < r)
        rsh = 0;
    else
        rsh = sh;

    return rsh;
}

// } Shift

// Branch and Bound Cyclic Edit Distance {

float lowerBound(float left_dist, float right_dist, int left, int right) {
    int delete_cost = 1;
    int insertion_cost = 1;
    return max(0., (left_dist + right_dist)/2.0 + (left - right)*(delete_cost + insertion_cost)/2.0);
}

class Element {
public:
    float lB;
    int l;
    int r;
    float l_dist;
    float r_dist;

    Element(float lB, int l, int r, float l_dist, float r_dist):
        lB(lB), l(l), r(r), l_dist(l_dist), r_dist(r_dist) {}
};

bool f(Element &x, Element &y) {
    return x.lB > y.lB;
}

float BBEd(Sample<t_value> &sx, Sample<t_value> &sy, float (*distance) (t_value, t_value)) {

    vector<t_value> x;
    vector<t_value> y;

    if (sx.values.size() > sy.values.size()) { // Cambiamos para poner la
                                               // cadena más pequeña abajo
        x = sy.values;
        y = sx.values;
    }
    else {
        x = sx.values;
        y = sy.values;
    }

    int lenx = x.size();
    int leny = y.size();

    matrix<t_value> D(2*lenx+1, leny+1);
    matrix< pair<int, int> > backpointer(2*lenx+1, leny+1);
    D[0][0] = 0;
    for (int i = 1; i < lenx+1; i++) {
        D[i][0] = D[i-1][0] + distance(x[i-1], LAMBDA);
        backpointer[i][0].first = 1;
        backpointer[i][0].second = 0;
    }
    for (int j = 1; j < leny+1; j++) {
        D[0][j] = D[0][j-1] + distance(LAMBDA, y[j-1]);
        backpointer[0][j].first = 0;
        backpointer[0][j].second = 1;
    }
    for (int i = 1; i < lenx+1; i++)
        for (int j = 1; j < leny+1; j++) {
            float d0 = D[i-1][j] + distance(x[i-1], LAMBDA);
            float d1 = D[i][j-1] + distance(LAMBDA, y[j-1]);
            float d2 = D[i-1][j-1] + distance(x[i-1], y[j-1]);
            if(d1 <= d0 and d1 <= d2) {
                D[i][j] = d1;
                backpointer[i][j].first = 0;
                backpointer[i][j].second = 1;
            }
            else if (d0 <= d2) {
                D[i][j] = d0;
                backpointer[i][j].first = 1;
                backpointer[i][j].second = 0;
            }
            else {
                D[i][j] = d2;
                backpointer[i][j].first = 1;
                backpointer[i][j].second = 1;
            }
        }

    if (x.size() == 1)
        return D[lenx][leny];

    matrix<int> mins(lenx+2, leny+1);
    matrix<int> maxs(lenx+2, leny+1);
    for (int i = 0; i < lenx+2; i++)
        for (int j = 0; j < leny+1; j++) {
            mins[i][j] = 2*lenx+1;
            maxs[i][j] = 0;
        }
    int i = lenx, j = leny;
    mins[0][j] = maxs[0][j] = i;

    int ri, rj;
    while (i != 0 and j != 0) {
        ri = backpointer[i][j].first;
        rj = backpointer[i][j].second;
        i -= ri;
        j -= rj;
        if (i < mins[0][j]) mins[0][j] = i;
        if (i > maxs[0][j]) maxs[0][j] = i;
    }

    for (int k=j; k>=0; --k) {
        mins[0][k] = 0;
        if (i > maxs[0][k]) maxs[0][k] = i;
    }

    for (int k=0; k<leny+1; k++) {
        mins[lenx][k] = mins[0][k] + lenx;
        maxs[lenx][k] = maxs[0][k] + lenx;
    }

    // X = x+x
    vector<t_value> X(2*lenx);
    for (int i = 0; i < lenx; i++) X[i] = x[i];
    for (int i = 0; i < lenx; i++) X[i+lenx] = x[i];

    float d_star = D[lenx][leny];
    int d_k = 0;

    vector<Element> S;
    S.push_back(Element(lowerBound(d_star, d_star, 0, lenx),
                        0, lenx,
                        d_star, d_star));
    make_heap(S.begin(), S.end(), f);

    while (S.size() != 0 and d_star > S[0].lB) {
        pop_heap(S.begin(), S.end(), f);
        Element e = *(S.end()-1);
        int l = e.l;
        int r = e.r;
        float l_dist = e.l_dist;
        float r_dist = e.r_dist;
        S.pop_back();
        int k = (l + r) / 2;

        float k_dist = xed(l, r, k, D, backpointer, X, x, y, distance, mins, maxs);
        if (d_star > k_dist) {
            d_star = k_dist;
            d_k = k;
        }
        float lk_lower_bound = lowerBound(l_dist, k_dist, l, k);
        float kr_lower_bound = lowerBound(k_dist, r_dist, k, r);

        if (k > l+1 and d_star > lk_lower_bound) {
            S.push_back(Element(lk_lower_bound, l, k, l_dist, k_dist));
            push_heap(S.begin(), S.end(), f);
        }
        if (r > k+1 and d_star > kr_lower_bound) {
            S.push_back(Element(kr_lower_bound, k, r, k_dist, r_dist));
            push_heap(S.begin(), S.end(), f);
        }
    }

    return d_star;
}

float BBEdExt(Sample<t_value> &sx, Sample<t_value> &sy, float (*distance) (t_value, t_value),
              float best_so_far_distance, unsigned long long &ldistances_count) {

    // just levenshtein {
    if (abs( (int)(sx.values.size() - sy.values.size()) ) >= best_so_far_distance)
        return INF;
    // }

    vector<t_value> x;
    vector<t_value> y;

    if (sx.values.size() > sy.values.size()) {
        x = sy.values;
        y = sx.values;
    }
    else {
        x = sx.values;
        y = sy.values;
    }

    int lenx = x.size();
    int leny = y.size();

    matrix<t_value> D(2*lenx+1, leny+1);
    matrix< pair<int, int> > backpointer(2*lenx+1, leny+1);
    D[0][0] = 0;
    for (int i = 1; i < lenx+1; i++) {
        D[i][0] = D[i-1][0] + distance(x[i-1], LAMBDA);
        backpointer[i][0].first = 1;
        backpointer[i][0].second = 0;
    }
    for (int j = 1; j < leny+1; j++) {
        D[0][j] = D[0][j-1] + distance(LAMBDA, y[j-1]);
        backpointer[0][j].first = 0;
        backpointer[0][j].second = 1;
    }
    for (int i = 1; i < lenx+1; i++) {
        for (int j = 1; j < leny+1; j++) {
            float d0 = D[i-1][j] + distance(x[i-1], LAMBDA);
            float d1 = D[i][j-1] + distance(LAMBDA, y[j-1]);
            float d2 = D[i-1][j-1] + distance(x[i-1], y[j-1]);
            if(d1 <= d0 and d1 <= d2) {
                D[i][j] = d1;
                backpointer[i][j].first = 0;
                backpointer[i][j].second = 1;
            }
            else if (d0 <= d2) {
                D[i][j] = d0;
                backpointer[i][j].first = 1;
                backpointer[i][j].second = 0;
            }
            else {
                D[i][j] = d2;
                backpointer[i][j].first = 1;
                backpointer[i][j].second = 1;
            }
        }
    }

    if (x.size() == 1)
        return D[lenx][leny];

    matrix<int> mins(lenx+2, leny+1);
    matrix<int> maxs(lenx+2, leny+1);
    for (int i = 0; i < lenx+2; i++)
        for (int j = 0; j < leny+1; j++) {
            mins[i][j] = 2*lenx+1;
            maxs[i][j] = 0;
        }
    int i = lenx, j = leny;
    mins[0][j] = maxs[0][j] = i;

    int ri, rj;
    while (i != 0 and j != 0) {
        ri = backpointer[i][j].first;
        rj = backpointer[i][j].second;
        i -= ri;
        j -= rj;
        if (i < mins[0][j]) mins[0][j] = i;
        if (i > maxs[0][j]) maxs[0][j] = i;
    }

    for (int k=j; k>=0; --k) {
        mins[0][k] = 0;
        if (i > maxs[0][k]) maxs[0][k] = i;
    }

    for (int k=0; k<leny+1; k++) {
        mins[lenx][k] = mins[0][k] + lenx;
        maxs[lenx][k] = maxs[0][k] + lenx;
    }

    vector<t_value> X(2*lenx);
    for (int i = 0; i < lenx; i++) X[i] = x[i];
    for (int i = 0; i < lenx; i++) X[i+lenx] = x[i];

    float d_star = D[lenx][leny];
    int d_k = 0;

    vector<Element> S;
    S.push_back(Element(lowerBound(d_star, d_star, 0, lenx),
                        0, lenx,
                        d_star, d_star));
    make_heap(S.begin(), S.end(), f);

    if (best_so_far_distance < D[lenx][leny])
        d_star = best_so_far_distance;
    else
        d_star = D[lenx][leny];

    while (S.size() != 0 and d_star > S[0].lB) {
        pop_heap(S.begin(), S.end(), f);
        Element e = *(S.end()-1);
        int l = e.l;
        int r = e.r;
        float l_dist = e.l_dist;
        float r_dist = e.r_dist;
        S.pop_back();
        int k = (l + r) / 2;

        float k_dist = xed(l, r, k, D, backpointer, X, x, y, distance, mins, maxs);
        if (d_star > k_dist) {
            d_star = k_dist;
            d_k = k;
        }
        float lk_lower_bound = lowerBound(l_dist, k_dist, l, k);
        float kr_lower_bound = lowerBound(k_dist, r_dist, k, r);

        if (k > l+1 and d_star > lk_lower_bound) {
            S.push_back(Element(lk_lower_bound, l, k, l_dist, k_dist));
            push_heap(S.begin(), S.end(), f);
        }
        if (r > k+1 and d_star > kr_lower_bound) {
            S.push_back(Element(kr_lower_bound, k, r, k_dist, r_dist));
            push_heap(S.begin(), S.end(), f);
        }
    }

    return d_star;
}

float BBEdExt1bubu(Sample<t_value> &sx, Sample<t_value> &sy, float (*distance) (t_value, t_value), float best_so_far_distance, unsigned long long &ldistances_count) {

    // just levenshtein {
    if (abs( (int)(sx.values.size() - sy.values.size()) ) >= best_so_far_distance)
        return INF;
    // }

    vector<t_value> x;
    vector<t_value> y;

    if (sx.values.size() > sy.values.size()) {
        x = sy.values;
        y = sx.values;
    }
    else {
        x = sx.values;
        y = sy.values;
    }

    int lenx = x.size();
    int leny = y.size();

    matrix<t_value> D(2*lenx+1, leny+1);

    vector<t_value> X(2*lenx);
    for (int i = 0; i < lenx; i++) X[i] = x[i];
    for (int i = 0; i < lenx; i++) X[i+lenx] = x[i];

    for (unsigned i = 0; i < lenx+1; i++) {
        D[i][0] = 0;
    }

    for (unsigned i = lenx+1; i < 2*lenx+1; i++) {
        D[i][0] = INF;
    }

    unsigned abort;

    for (unsigned j = 1; j < leny+1; j++) {

        abort = 1;

        D[0][j] = D[0][j-1] + distance(LAMBDA, y[j-1]);
        if (D[0][j] < best_so_far_distance)
            abort = 0;

        for (unsigned i = 1; i < 2*lenx+1; i++) {

            float d0 = D[i-1][j] + distance(X[i-1], LAMBDA);
            float d1 = D[i][j-1] + distance(LAMBDA, y[j-1]);
            float d2 = D[i-1][j-1] + distance(X[i-1], y[j-1]);

            if (d1 <= d0 and d1 <= d2) {
                D[i][j] = d1;
            }
            else if (d0 <= d2) {
                D[i][j] = d0;
            }
            else {
                D[i][j] = d2;
            }

            if (D[i][j] < best_so_far_distance)
                abort = 0;
        }

        if (abort) {
            return INF;
        }
    }

    matrix< pair<int, int> > backpointer(2*lenx+1, leny+1);
    D[0][0] = 0;
    for (int i = 1; i < lenx+1; i++) {
        D[i][0] = D[i-1][0] + distance(x[i-1], LAMBDA);
        backpointer[i][0].first = 1;
        backpointer[i][0].second = 0;
    }
    for (int j = 1; j < leny+1; j++) {
        D[0][j] = D[0][j-1] + distance(LAMBDA, y[j-1]);
        backpointer[0][j].first = 0;
        backpointer[0][j].second = 1;
    }
    for (int i = 1; i < lenx+1; i++) {
        for (int j = 1; j < leny+1; j++) {
            float d0 = D[i-1][j] + distance(x[i-1], LAMBDA);
            float d1 = D[i][j-1] + distance(LAMBDA, y[j-1]);
            float d2 = D[i-1][j-1] + distance(x[i-1], y[j-1]);
            if(d1 <= d0 and d1 <= d2) {
                D[i][j] = d1;
                backpointer[i][j].first = 0;
                backpointer[i][j].second = 1;
            }
            else if (d0 <= d2) {
                D[i][j] = d0;
                backpointer[i][j].first = 1;
                backpointer[i][j].second = 0;
            }
            else {
                D[i][j] = d2;
                backpointer[i][j].first = 1;
                backpointer[i][j].second = 1;
            }
        }
    }

    if (x.size() == 1)
        return D[lenx][leny];

    matrix<int> mins(lenx+2, leny+1);
    matrix<int> maxs(lenx+2, leny+1);
    for (int i = 0; i < lenx+2; i++)
        for (int j = 0; j < leny+1; j++) {
            mins[i][j] = 2*lenx+1;
            maxs[i][j] = 0;
        }
    int i = lenx, j = leny;
    mins[0][j] = maxs[0][j] = i;

    int ri, rj;
    while (i != 0 and j != 0) {
        ri = backpointer[i][j].first;
        rj = backpointer[i][j].second;
        i -= ri;
        j -= rj;
        if (i < mins[0][j]) mins[0][j] = i;
        if (i > maxs[0][j]) maxs[0][j] = i;
    }

    for (int k=j; k>=0; --k) {
        mins[0][k] = 0;
        if (i > maxs[0][k]) maxs[0][k] = i;
    }

    for (int k=0; k<leny+1; k++) {
        mins[lenx][k] = mins[0][k] + lenx;
        maxs[lenx][k] = maxs[0][k] + lenx;
    }


    float d_star = D[lenx][leny];
    int d_k = 0;

    vector<Element> S;
    S.push_back(Element(lowerBound(d_star, d_star, 0, lenx),
                        0, lenx,
                        d_star, d_star));
    make_heap(S.begin(), S.end(), f);

    if (best_so_far_distance < D[lenx][leny])
        d_star = best_so_far_distance;
    else
        d_star = D[lenx][leny];

    while (S.size() != 0 and d_star > S[0].lB) {
        pop_heap(S.begin(), S.end(), f);
        Element e = *(S.end()-1);
        int l = e.l;
        int r = e.r;
        float l_dist = e.l_dist;
        float r_dist = e.r_dist;
        S.pop_back();
        int k = (l + r) / 2;

        float k_dist = xed(l, r, k, D, backpointer, X, x, y, distance, mins, maxs);
        if (d_star > k_dist) {
            d_star = k_dist;
            d_k = k;
        }
        float lk_lower_bound = lowerBound(l_dist, k_dist, l, k);
        float kr_lower_bound = lowerBound(k_dist, r_dist, k, r);

        if (k > l+1 and d_star > lk_lower_bound) {
            S.push_back(Element(lk_lower_bound, l, k, l_dist, k_dist));
            push_heap(S.begin(), S.end(), f);
        }
        if (r > k+1 and d_star > kr_lower_bound) {
            S.push_back(Element(kr_lower_bound, k, r, k_dist, r_dist));
            push_heap(S.begin(), S.end(), f);
        }
    }

    return d_star;
}


// Branch and Bound Cyclic Edit Distance }


