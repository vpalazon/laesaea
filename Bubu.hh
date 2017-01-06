/*===================================================================================

  File: Bubu.hh

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

#ifndef _BUBU_HH
#define _BUBU_HH

#include "Samples.hh"
#include <vector>
#include <set>
#include "matrix.h"
#include "Heap.hh"
#include "Defs.hh"

template <typename T> float bb(Sample<T> &sx, Sample<T> &sy, float (*distance) (T, T));

template <typename T> float bubu(Sample<T> &sx, Sample<T> &sy, float (*distance) (T, T), 
                                 float best_so_far_distance=INF);

template <typename T>
float bububs(Sample<T> &sx, Sample<T> &sy, float (*distance) (T, T), float best_so_far_distance=INF);

template <typename T>
float bb(Sample<T> &sx, Sample<T> &sy, float (*distance) (T, T)) {
  vector<T> &x = sx.values;
  vector<T> &y = sy.values;
  
  matrix<float> D = matrix<float>(2*x.size(), y.size());
  float dis;
  float min_dis = INF;

  matrix<float> d(2*x.size(), y.size());
  for (unsigned i = 0; i < x.size(); i++)
    for (unsigned j = 0; j < y.size(); j++)
      d[i+x.size()][j] = d[i][j] = distance(x[i], y[j]);

  for (unsigned i = 0; i < x.size(); i++) D[i][0] = d[i][0];
  for (unsigned i = x.size(); i < 2*x.size(); i++) D[i][0] = D[i-1][0] + d[i][0];
  for (unsigned j = 1; j < y.size(); j++) D[0][j] = D[0][j-1] + d[0][j];
  for (unsigned i = 1; i < 2*x.size(); i++)
    for (unsigned j = 1; j < y.size(); j++) {
      float a=D[i-1][j-1], b = D[i-1][j], c = D[i][j-1];
      dis = d[i][j];
      if (a <= b)
	if (a <= c) D[i][j] = a+dis;
	else D[i][j] = c+dis;
      else
	if (b <= c) D[i][j] = b+dis;
	else D[i][j] = c+dis;
    }
  for (unsigned i = x.size()-1; i < 2*x.size(); i++) {
    if (D[i][y.size()-1] < min_dis)
      min_dis = D[i][y.size()-1];
  }
  return min_dis;
}

template <typename T>
float bubu(Sample<T> &sx, Sample<T> &sy, float (*distance) (T, T), float best_so_far_distance=INF) {
    vector<T> &x = sx.values;
    vector<T> &y = sy.values;

    unsigned lenx = x.size();
    unsigned leny = y.size();

    matrix<float> d(2*lenx+1, leny);
    matrix<float> D(2*lenx+1, leny);

    for (unsigned i = 0; i < lenx; i++)
        D[i][0] = d[i][0] = distance(x[i], y[0]);

    for (unsigned i = lenx; i < 2*lenx; i++) {
        d[i][0] = d[i-lenx][0];
        D[i][0] = D[i-1][0] + d[i][0];
    }

    for (unsigned j = 1; j < leny; j++) {
        d[0][j] = distance(x[0], y[j]);
        D[0][j] = D[0][j-1] + d[0][j];
    }
    
    for (unsigned j = 1; j < leny; j++) {
        unsigned abort = 1;

        for (unsigned i = 1; i < lenx; i++) {

            float a = D[i-1][j-1];
            float b = D[i-1][j];
            float c = D[i][j-1];
            float dis = d[i][j] = distance(x[i], y[j]);

            if (a <= b) 
                if (a <= c)
                    D[i][j] = a+dis;
                else
                    D[i][j] = c+dis;
            else
                if (b <= c)
                    D[i][j] = b+dis;
                else
                    D[i][j] = c+dis;
            if (D[i][j] < best_so_far_distance)
                abort = 0;
        }

        for (unsigned i = lenx; i < 2*lenx; i++) {

            float a = D[i-1][j-1];
            float b = D[i-1][j];
            float c = D[i][j-1];
            float dis = d[i][j] = d[i-lenx][j];

            if (a <= b) 
                if (a <= c)
                    D[i][j] = a+dis;
                else
                    D[i][j] = c+dis;
            else
                if (b <= c)
                    D[i][j] = b+dis;
                else
                    D[i][j] = c+dis;
            if (D[i][j] < best_so_far_distance)
                abort = 0;
        }

        if (abort)
            return INF;
    }
    
    float dist = INF;
    for (unsigned i=lenx-1; i<2*lenx; i++) {
        if (D[i][leny-1] < dist)
            dist = D[i][leny-1];
    }
    return dist;
}

typedef pair<unsigned, unsigned> coord;
typedef pair<unsigned, float> keyscore;

void succs(coord &u, unsigned limx, unsigned limy, vector<coord> &s);

unsigned coord2key(unsigned i, unsigned j, unsigned m);
void key2coord(unsigned k, unsigned m, 
               unsigned &i, unsigned &j);

unsigned setContains(set<coord> &s, coord &item);

template <typename T>
float bububs(Sample<T> &sx, Sample<T> &sy, float (*distance) (T, T), float best_so_far_distance=INF) {
    vector<T> &x = sx.values;
    vector<T> &y = sy.values;

    unsigned lenx = x.size();
    unsigned leny = y.size();

    unsigned m = 2*lenx+1;

    matrix<float> d(m, leny);
    matrix<float> D(m, leny);

    set<coord> S;
    for (unsigned i=0; i<m; i++)
        for (unsigned j=0; j<leny; j++)
            D[i][j] = INF;
    set<coord> E;
    
    unsigned key;
    float score;
    
    coord u, v;
    vector<coord> lv;
    
    for (unsigned i=lenx-1; i<2*lenx; i++)
        E.insert(coord(i, leny-1));
    
    for (unsigned i=0; i<lenx; i++)
        D[i][0] = d[i][0] = distance(x[i], y[0]);

    vector<keyscore> v_key_score;
    for (unsigned i=0; i<lenx; i++)
        v_key_score.push_back(keyscore(coord2key(i, 0, m), D[i][0]));

    IndexedMinHeap Q = IndexedMinHeap(leny*(m), v_key_score);
    // IndexedMinHeap Q;
    // Q.init(leny*m);
    // for (unsigned i=0; i<lenx; i++)
    //     Q.insert(coord2key(i, 0, m), D[i][0]);
    
    while (not Q.is_empty()) {
        Q.extract_min(key, score);
        key2coord(key, m, u.first, u.second);
        if (setContains(E, u))
            break;
        S.insert(u);
        succs(u, 2*lenx, leny, lv);
        for (unsigned i=0; i<lv.size(); i++)
            if (not setContains(S, lv[i])) {
                unsigned di;
                unsigned dj = lv[i].second;
                if (lv[i].first >= lenx)
                    di = lv[i].first - lenx;
                else
                    di = lv[i].first;
                unsigned lv_key = coord2key(lv[i].first, lv[i].second, m);
                if (not Q.contains(lv_key)) {
                    d[di][dj] = distance(x[di], y[dj]);
                    D[lv[i].first][lv[i].second] = D[u.first][u.second] + d[di][dj];
                    Q.insert(lv_key, D[lv[i].first][lv[i].second]);
                }
                else {
                    float aux;
                    aux = D[u.first][u.second] + d[di][dj];
                    if (aux < D[lv[i].first][lv[i].second]) {
                        D[lv[i].first][lv[i].second] = aux;
                        Q.decrease(lv_key, D[lv[i].first][lv[i].second]);
                    }
                }
            }
    }
    return D[u.first][u.second];
}

#endif

