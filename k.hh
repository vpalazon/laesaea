/*===================================================================================

  File: k.hh

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

#ifndef _K_HH
#define _K_HH

#include <list>
#include <vector>
#include <algorithm>
#include <iostream>

#include "Samples.hh"
#include "matrix.h"
#include "Defs.hh"
#include "Elem.hh"
#include "KBestVector.hh"
#include "debug.hh"
#include "LocalDistances.hh"

using namespace std;

// K vecinos más próximos

template <class T>
void k(list<int> &P, unsigned x, Samples<T> &seqs,
       float (*f)(Sample<T> &sx, 
                  Sample<T> &sy, 
                  float (*dist) (T, T)),
       float (*dist) (T, T),
       // Salida:
       vector< pair<unsigned, float> > &knn,
       // Optional:
       float omega=INF);

template <class T>
void k_heap(list<int> &P, unsigned x, Samples<T> &seqs,
            float (*f)(Sample<T> &sx, 
                       Sample<T> &sy, 
                       float (*dist) (T, T)),
            float (*dist) (T, T),
            // Salida:
            vector< pair<unsigned, float> > &knn);

template <class T>
void k_sort(list<int> &P, unsigned x, Samples<T> &seqs,
            float (*f)(Sample<T> &sx, 
                       Sample<T> &sy, 
                       float (*dist) (T, T)),
            float (*dist) (T, T),
            // Salida:
            vector< pair<unsigned, float> > &knn);

bool compk(int i1, int i2);

template <class T>
void k(list<int> &P, unsigned x, Samples<T> &seqs,
       float (*f)(Sample<T> &sx, 
                  Sample<T> &sy, 
                  float (*dist) (T, T)),
       float (*dist) (T, T),
       // Salida:
       vector< pair<unsigned, float> > &knn,
       // Output:
       float omega=INF) {

    unsigned K = knn.size();

    KBestVector KB(K);
    
    int s = 0;                  // Elemento arbitrario de P
    list<int>::iterator s_it;

    float (*fomega)(Sample<T> &sx, 
                    Sample<T> &sy, 
                    float (*dist) (T, T), float omega) = NULL;

    // if (omega != INF) {
    //     fomega = maesSakoe;
    // }

    while (P.size() > 0) {
        s_it = P.begin();
        advance(s_it, s);
        
        float d_s;
        if (omega == INF)
            d_s = (*f)((*seqs.samples)[x], (*seqs.samples)[*s_it], dist);
        else
            d_s = (*fomega)((*seqs.samples)[x], (*seqs.samples)[*s_it], dist, omega);
        
        if (d_s <= KB.top()) {
            KB.insert(*s_it, d_s);
        }

        P.erase(s_it);
    }

    KB.get_sort(knn);
}

template <class T>
void k_heap(list<int> &P, unsigned x, Samples<T> &seqs,
            float (*f)(Sample<T> &sx, 
                       Sample<T> &sy, 
                       float (*dist) (T, T)),
            float (*dist) (T, T),
            // Salida:
            vector< pair<unsigned, float> > &knn) {

    unsigned K = knn.size();

    vector<Elem *> NN;
    for (unsigned k=0; k<K; k++)
        NN.push_back(new Elem(-1, INF));

    int s = 0;                  // Elemento arbitrario de P
    list<int>::iterator s_it;

    while (P.size() > 0) {
        s_it = P.begin();
        advance(s_it, s);
        float d_s = (*f)((*seqs.samples)[x], (*seqs.samples)[*s_it], dist);
        
        if (d_s < NN[0]->cost) {       // Actualización del heap de los nns
            NN[0]->e = *s_it;
            NN[0]->cost = d_s;
            make_heap(NN.begin(), NN.end(), fMaxHeap);
        }

        P.erase(s_it);
    }

    sort_heap(NN.begin(), NN.end(), fMaxHeap);
    for (unsigned i=0; i<NN.size(); i++)
        knn[i] = NN[i]->e;
}

template <class T>
void k_sort(list<int> &P, unsigned x, Samples<T> &seqs,
            float (*f)(Sample<T> &sx, 
                       Sample<T> &sy, 
                       float (*dist) (T, T)),
            float (*dist) (T, T),
            // Salida:
            vector< pair<unsigned, float> > &knn) {
    
    unsigned K = knn.size();

    vector<Elem *> NN;

    int s = 0;                  // Elemento arbitrario de P
    list<int>::iterator s_it;

    while (P.size() > 0) {
        s_it = P.begin();
        advance(s_it, s);
        float d_s = (*f)((*seqs.samples)[x], (*seqs.samples)[*s_it], dist);
        
        NN.push_back(new Elem(*s_it, d_s));
        
        P.erase(s_it);
    }

    sort(NN.begin(), NN.end(), fMaxHeap);

    for (unsigned i=0; i<K; i++)
        knn[i] = NN[i]->e;
}


#endif

