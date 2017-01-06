/*===================================================================================

  File: Aesa.hh

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

#ifndef _AESA_HH
#define _AESA_HH

#include <list>
#include <vector>
#include <algorithm>
#include <iostream>
#include <ext/hash_set>
#include <set>

#include "Samples.hh"
#include "matrix.h"
#include "Defs.hh"
#include "Elem.hh"
#include "KBestVector.hh"
#include "debug.hh"
#include "Ed.hh"


using namespace std;
using namespace __gnu_cxx;

int isin(int v, set<int> &P);

template <class T>
void kAesa(set<int> &P, unsigned x, matrix<float> &D, Samples<T> &seqs,
           float (*f)(Sample<T> &sx, 
                      Sample<T> &sy, 
                      float (*dist) (T, T)),
           float (*dist) (T, T),
           vector< pair<unsigned, float> > &knn, unsigned long &bound_count,
           float H=0, float omega=INF);

template <class T>
void klAesa1(set<int> &P, set<int> &B, unsigned x, matrix<float> &D, Samples<T> &seqs,
             float (*f)(Sample<T> &sx, 
                        Sample<T> &sy, 
                        float (*dist) (T, T)),
             float (*dist) (T, T),
             vector< pair<unsigned, float> > &knn, 
             unsigned long &bound_count, unsigned long long &ldistances_count,
             float H=0, float omega=INF);

template <class T>
void kblAesa1(set<int> &P, set<int> &B, unsigned x, matrix<float> &D, Samples<T> &seqs,
              float (*f)(Sample<T> &sx, 
                         Sample<T> &sy, 
                         float (*dist) (T, T)),
              float (*dist) (T, T),
              vector< pair<unsigned, float> > &knn, 
              unsigned long &bound_count, unsigned long long &ldistances_count,
              float H=0, float omega=INF);

template <class T>
void kblAesa1BBEdExt(set<int> &P, set<int> &B, unsigned x, matrix<float> &D, Samples<T> &seqs,
                     float (*f)(Sample<T> &sx, 
                                Sample<T> &sy, 
                                float (*dist) (T, T)),
                     float (*dist) (T, T),
                     vector< pair<unsigned, float> > &knn, 
                     unsigned long &bound_count, unsigned long long &ldistances_count,
                     float H=0, float omega=INF);

template <class T>
void kblAesa1BBEdExt1bubu(set<int> &P, set<int> &B, unsigned x, matrix<float> &D, Samples<T> &seqs,
                          float (*f)(Sample<T> &sx, 
                                     Sample<T> &sy, 
                                     float (*dist) (T, T)),
                          float (*dist) (T, T),
                          vector< pair<unsigned, float> > &knn, 
                          unsigned long &bound_count, unsigned long long &ldistances_count,
                          float H=0, float omega=INF);

template <class T>
void kblAesa1maesEdExt1bubu(set<int> &P, set<int> &B, unsigned x, matrix<float> &D, Samples<T> &seqs,
                            float (*f)(Sample<T> &sx, 
                                       Sample<T> &sy, 
                                       float (*dist) (T, T)),
                            float (*dist) (T, T),
                            vector< pair<unsigned, float> > &knn, 
                            unsigned long &bound_count, unsigned long long &ldistances_count,
                            float H=0, float omega=INF);

template <class T>
void kAesa(set<int> &P, unsigned x, matrix<float> &D, Samples<T> &seqs,
           float (*f)(Sample<T> &sx, 
                      Sample<T> &sy, 
                      float (*dist) (T, T)),
           float (*dist) (T, T),
           vector< pair<unsigned, float> > &knn, unsigned long &bound_count,
           float H=0, float omega=INF) {

    unsigned K = knn.size();

    vector<float> G(P.size()+1);
    for (unsigned p=0; p<G.size(); p++)
        G[p] = 0;
    
    KBestVector KB(K);
    
    int s = *(P.begin());

    bound_count = 0;
    while (P.size() > 0) {

        float d_s = (*f)((*seqs.samples)[x], (*seqs.samples)[s], dist);
        P.erase(s);

        if (d_s <= KB.top()) {
            KB.insert(s, d_s);
        }
        
        int siguiente = -1;
        float gmin = INF;
        for (set<int>::iterator p_it=P.begin(); p_it!=P.end(); p_it++) {
            int p = *p_it;
            G[p] = max(G[p], abs(D[s][p] - d_s) + H);
            if (G[p] > KB.top()) {
                P.erase(p);
            }
            else
                if (G[p] < gmin) {
                    gmin = G[p];
                    siguiente = p;
                }
        }
        
        s = siguiente;
        
        bound_count++;
    }
    
    KB.get_sort(knn);
}

template <class T>
void klAesa1(set<int> &P, set<int> &B, unsigned x, matrix<float> &D, Samples<T> &seqs,
             float (*f)(Sample<T> &sx, 
                        Sample<T> &sy, 
                        float (*dist) (T, T)),
             float (*dist) (T, T),
             vector< pair<unsigned, float> > &knn, 
             unsigned long &bound_count, unsigned long long &ldistances_count,
             float H=0, float omega=INF) {

    unsigned K = knn.size();

    vector<float> G(P.size()+1);
    for (unsigned p=0; p<G.size(); p++)
        G[p] = 0;

    KBestVector KB(K);

    int s = *(B.begin());

    bound_count = 0;
    while (P.size() > 0) {
        float d_s = (*f)((*seqs.samples)[x], (*seqs.samples)[s], dist);
        P.erase(s);

        if (d_s <= KB.top()) {
            KB.insert(s, d_s);
        }
        
        int siguiente_B = -1;
        float gmin_B = INF;
        int siguiente = -1;
        float gmin = INF;
        for (set<int>::iterator p_it=P.begin(); p_it!=P.end(); p_it++) {
            int p = *p_it;
            if (isin(s, B)) {
                G[p] = max(G[p], abs(D[s][p] - d_s) + H);
            }

            int p_is_in_B = isin(p, B);

            if ((not p_is_in_B) and G[p] > KB.top()) {
                P.erase(p);                           
            }
            else {
                if (p_is_in_B) {
                    if (G[p] < gmin_B) {
                        gmin_B = G[p];
                        siguiente_B = p;
                    }
                }
                else
                    if (G[p] < gmin) {
                        gmin = G[p];
                        siguiente = p;
                    }
            }
        }
        
        if (siguiente_B != -1) 
            s = siguiente_B;   
        else
            s = siguiente;     
        
        bound_count++;
    }

    ldistances_count = bound_count * (*seqs.samples)[0].size() * (*seqs.samples)[0].size();
    
    KB.get_sort(knn);
}

template <class T>
void kblAesa1(set<int> &P, set<int> &B, unsigned x, matrix<float> &D, Samples<T> &seqs,
              float (*f)(Sample<T> &sx, 
                         Sample<T> &sy, 
                         float (*dist) (T, T)),
              float (*dist) (T, T),
              vector< pair<unsigned, float> > &knn, 
              unsigned long &bound_count, unsigned long long &ldistances_count,
              float H=0, float omega=INF) {

    unsigned K = knn.size();

    vector<float> G(P.size()+1);
    for (unsigned p=0; p<G.size(); p++)
        G[p] = 0;

    KBestVector KB(K);

    int s = *(B.begin());

    bound_count = 0;
    ldistances_count = 0;
    while (P.size() > 0) {
        float d_s;
        if (isin(s, B)) {
            d_s = (*f)((*seqs.samples)[x], (*seqs.samples)[s], dist);
            bound_count++;
        }
        else {
            d_s = cdtw1bubu((*seqs.samples)[x], (*seqs.samples)[s], dist, KB.top(), ldistances_count);
            if (d_s != INF)
                bound_count++;
        }
        
        P.erase(s);

        if (d_s <= KB.top()) {
            KB.insert(s, d_s);
        }
        
        int siguiente_B = -1;
        float gmin_B = INF;
        int siguiente = -1;
        float gmin = INF;
        for (set<int>::iterator p_it=P.begin(); p_it!=P.end(); p_it++) {
            int p = *p_it;
            if (isin(s, B)) {
                G[p] = max(G[p], abs(D[s][p] - d_s) + H);
            }

            int p_is_in_B = isin(p, B);

            if ((not p_is_in_B) and G[p] > KB.top()) {
                P.erase(p);                           
            }
            else {
                if (p_is_in_B) {
                    if (G[p] < gmin_B) {
                        gmin_B = G[p];
                        siguiente_B = p;
                    }
                }
                else
                    if (G[p] < gmin) {
                        gmin = G[p];
                        siguiente = p;
                    }
            }
        }
        
        if (siguiente_B != -1) 
            s = siguiente_B;
        else
            s = siguiente;
        
    }
    
    KB.get_sort(knn);
}

template <class T>
void kblAesa1BBEdExt(set<int> &P, set<int> &B, unsigned x, matrix<float> &D, Samples<T> &seqs,
                     float (*f)(Sample<T> &sx, 
                                Sample<T> &sy, 
                                float (*dist) (T, T)),
                     float (*dist) (T, T),
                     vector< pair<unsigned, float> > &knn, 
                     unsigned long &bound_count, unsigned long long &ldistances_count,
                     float H=0, float omega=INF) {

    unsigned K = knn.size();

    vector<float> G(P.size()+1);
    for (unsigned p=0; p<G.size(); p++)
        G[p] = 0;

    KBestVector KB(K);

    int s = *(B.begin());

    bound_count = 0;
    ldistances_count = 0;
    while (P.size() > 0) {
        float d_s;
        if (isin(s, B)) {
            d_s = (*f)((*seqs.samples)[x], (*seqs.samples)[s], dist);
            bound_count++;
        }
        else {
            d_s = BBEdExt((*seqs.samples)[x], (*seqs.samples)[s], dist, KB.top(), ldistances_count);
            if (d_s == KB.top()) {
                d_s = INF;
            }
            else
                bound_count++;
        }
        
        P.erase(s);

        if (d_s <= KB.top()) {
            KB.insert(s, d_s);
        }
        
        int siguiente_B = -1;
        float gmin_B = INF;
        int siguiente = -1;
        float gmin = INF;
        for (set<int>::iterator p_it=P.begin(); p_it!=P.end(); p_it++) {
            int p = *p_it;
            if (isin(s, B)) {
                G[p] = max(G[p], abs(D[s][p] - d_s) + H);
            }

            int p_is_in_B = isin(p, B);

            if ((not p_is_in_B) and G[p] > KB.top()) {
                P.erase(p);                           
            }
            else {
                if (p_is_in_B) {
                    if (G[p] < gmin_B) {
                        gmin_B = G[p];
                        siguiente_B = p;
                    }
                }
                else
                    if (G[p] < gmin) {
                        gmin = G[p];
                        siguiente = p;
                    }
            }
        }
        
        if (siguiente_B != -1) 
            s = siguiente_B;
        else
            s = siguiente;
        
    }
    
    KB.get_sort(knn);
}

template <class T>
void kblAesa1BBEdExt1bubu(set<int> &P, set<int> &B, unsigned x, matrix<float> &D, Samples<T> &seqs,
                          float (*f)(Sample<T> &sx, 
                                     Sample<T> &sy, 
                                     float (*dist) (T, T)),
                          float (*dist) (T, T),
                          vector< pair<unsigned, float> > &knn, 
                          unsigned long &bound_count, unsigned long long &ldistances_count,
                          float H=0, float omega=INF) {

    unsigned K = knn.size();

    vector<float> G(P.size()+1);
    for (unsigned p=0; p<G.size(); p++)
        G[p] = 0;

    KBestVector KB(K);

    int s = *(B.begin());

    bound_count = 0;
    ldistances_count = 0;
    while (P.size() > 0) {
        float d_s;
        if (isin(s, B)) {
            d_s = (*f)((*seqs.samples)[x], (*seqs.samples)[s], dist);
            bound_count++;
        }
        else {
            d_s = BBEdExt1bubu((*seqs.samples)[x], (*seqs.samples)[s], dist, KB.top(), ldistances_count);
            if (d_s == KB.top() or d_s == INF) {
                d_s = INF;
            }
            else
                bound_count++;
        }
        
        P.erase(s);

        if (d_s <= KB.top()) {
            KB.insert(s, d_s);
        }
        
        int siguiente_B = -1;
        float gmin_B = INF;
        int siguiente = -1;
        float gmin = INF;
        for (set<int>::iterator p_it=P.begin(); p_it!=P.end(); p_it++) {
            int p = *p_it;
            if (isin(s, B)) {
                G[p] = max(G[p], abs(D[s][p] - d_s) + H);
            }

            int p_is_in_B = isin(p, B);

            if ((not p_is_in_B) and G[p] > KB.top()) {
                P.erase(p);                           
            }
            else {
                if (p_is_in_B) {
                    if (G[p] < gmin_B) {
                        gmin_B = G[p];
                        siguiente_B = p;
                    }
                }
                else
                    if (G[p] < gmin) {
                        gmin = G[p];
                        siguiente = p;
                    }
            }
        }
        
        if (siguiente_B != -1) 
            s = siguiente_B;
        else
            s = siguiente;
        
    }
    
    KB.get_sort(knn);
}

template <class T>
void kblAesa1maesEdExt1bubu(set<int> &P, set<int> &B, unsigned x, matrix<float> &D, Samples<T> &seqs,
                            float (*f)(Sample<T> &sx, 
                                       Sample<T> &sy, 
                                       float (*dist) (T, T)),
                            float (*dist) (T, T),
                            vector< pair<unsigned, float> > &knn, 
                            unsigned long &bound_count, unsigned long long &ldistances_count,
                            float H=0, float omega=INF) {

    unsigned K = knn.size();

    vector<float> G(P.size()+1);
    for (unsigned p=0; p<G.size(); p++)
        G[p] = 0;

    KBestVector KB(K);

    int s = *(B.begin());

    bound_count = 0;
    ldistances_count = 0;
    while (P.size() > 0) {
        float d_s;
        if (isin(s, B)) {
            d_s = (*f)((*seqs.samples)[x], (*seqs.samples)[s], dist);
            bound_count++;
        }
        else {
            d_s = maesEdExt1bubu((*seqs.samples)[x], (*seqs.samples)[s], dist, KB.top(), ldistances_count);
            if (d_s == KB.top() or d_s == INF) {
                d_s = INF;
            }
            else
                bound_count++;
        }
        
        P.erase(s);

        if (d_s <= KB.top()) {
            KB.insert(s, d_s);
        }
        
        int siguiente_B = -1;
        float gmin_B = INF;
        int siguiente = -1;
        float gmin = INF;
        for (set<int>::iterator p_it=P.begin(); p_it!=P.end(); p_it++) {
            int p = *p_it;
            if (isin(s, B)) {
                G[p] = max(G[p], abs(D[s][p] - d_s) + H);
            }

            int p_is_in_B = isin(p, B);

            if ((not p_is_in_B) and G[p] > KB.top()) {
                P.erase(p);                          
            }
            else {
                if (p_is_in_B) {
                    if (G[p] < gmin_B) {
                        gmin_B = G[p];
                        siguiente_B = p;
                    }
                }
                else
                    if (G[p] < gmin) {
                        gmin = G[p];
                        siguiente = p;
                    }
            }
        }
        
        if (siguiente_B != -1) 
            s = siguiente_B;
        else
            s = siguiente;
        
    }
    
    KB.get_sort(knn);
}

template <class T>
void kblAesa1EdExt(set<int> &P, set<int> &B, unsigned x, matrix<float> &D, Samples<T> &seqs,
                   float (*f)(Sample<T> &sx, 
                              Sample<T> &sy, 
                              float (*dist) (T, T)),
                   float (*dist) (T, T),
                   vector< pair<unsigned, float> > &knn, 
                   unsigned long &bound_count, unsigned long long &ldistances_count,
                   float H=0, float omega=INF) {

    unsigned K = knn.size();

    vector<float> G(P.size()+1);
    for (unsigned p=0; p<G.size(); p++)
        G[p] = 0;

    KBestVector KB(K);

    int s = *(B.begin());

    bound_count = 0;
    ldistances_count = 0;
    while (P.size() > 0) {
        float d_s;
        if (isin(s, B)) {
            d_s = (*f)((*seqs.samples)[x], (*seqs.samples)[s], dist);
            bound_count++;
        }
        else {
            d_s = edExt((*seqs.samples)[x], (*seqs.samples)[s], dist, KB.top(), ldistances_count);
            if (d_s != INF)
                bound_count++;
        }
        
        P.erase(s);

        if (d_s <= KB.top()) {
            KB.insert(s, d_s);
        }
        
        int siguiente_B = -1;
        float gmin_B = INF;
        int siguiente = -1;
        float gmin = INF;
        for (set<int>::iterator p_it=P.begin(); p_it!=P.end(); p_it++) {
            int p = *p_it;
            if (isin(s, B)) {
                G[p] = max(G[p], abs(D[s][p] - d_s) + H);
            }

            int p_is_in_B = isin(p, B);

            if ((not p_is_in_B) and G[p] > KB.top()) {
                P.erase(p);                 
            }
            else {
                if (p_is_in_B) {
                    if (G[p] < gmin_B) {
                        gmin_B = G[p];
                        siguiente_B = p;
                    }
                }
                else
                    if (G[p] < gmin) {
                        gmin = G[p];
                        siguiente = p;
                    }
            }
        }
        
        if (siguiente_B != -1) 
            s = siguiente_B;
        else
            s = siguiente;
        
    }
    
    KB.get_sort(knn);
}

template <class T>
void kblAesa1EdExt2(set<int> &P, set<int> &B, unsigned x, matrix<float> &D, Samples<T> &seqs,
                   float (*f)(Sample<T> &sx, 
                              Sample<T> &sy, 
                              float (*dist) (T, T)),
                   float (*dist) (T, T),
                   vector< pair<unsigned, float> > &knn, 
                   unsigned long &bound_count, unsigned long long &ldistances_count,
                   float H=0, float omega=INF) {

    unsigned K = knn.size();

    vector<float> G(P.size()+1);
    for (unsigned p=0; p<G.size(); p++)
        G[p] = 0;

    KBestVector KB(K);

    int s = *(B.begin());

    bound_count = 0;
    ldistances_count = 0;
    while (P.size() > 0) {
        float d_s;
        if (isin(s, B)) {
            d_s = (*f)((*seqs.samples)[x], (*seqs.samples)[s], dist);
            bound_count++;
        }
        else {
            d_s = edExt2((*seqs.samples)[x], (*seqs.samples)[s], dist, KB.top(), ldistances_count);
            if (d_s != INF)
                bound_count++;
        }
        
        P.erase(s);

        if (d_s <= KB.top()) {
            KB.insert(s, d_s);
        }
        
        int siguiente_B = -1;
        float gmin_B = INF;
        int siguiente = -1;
        float gmin = INF;
        for (set<int>::iterator p_it=P.begin(); p_it!=P.end(); p_it++) {
            int p = *p_it;
            if (isin(s, B)) {
                G[p] = max(G[p], abs(D[s][p] - d_s) + H);
            }

            int p_is_in_B = isin(p, B);

            if ((not p_is_in_B) and G[p] > KB.top()) {
                P.erase(p);                           
            }
            else {
                if (p_is_in_B) {
                    if (G[p] < gmin_B) {
                        gmin_B = G[p];
                        siguiente_B = p;
                    }
                }
                else
                    if (G[p] < gmin) {
                        gmin = G[p];
                        siguiente = p;
                    }
            }
        }
        
        if (siguiente_B != -1) 
            s = siguiente_B;
        else
            s = siguiente;
        
    }
    
    KB.get_sort(knn);
}

#endif

