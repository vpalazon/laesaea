/*===================================================================================

  File: tc.cc

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
#include <ext/hash_set>
#include <set>

#include "Ed.hh"
#include "Samples.hh"
#include "LocalDistances.hh"
#include "Defs.hh"
#include "Tt.hh"
#include "Aesa.hh"
#include "k.hh"
#include "Useful.hh"
#include <getopt.h>

#include <algorithm>        // random_shuffle

extern "C" {
#include "chronometer.h"
}

#include "debug.hh"

using namespace std;
using namespace __gnu_cxx;

template <class T> void experiment_kReal(Samples<T> &seqs,
                                         float (*f)(Sample<T> &sx, 
                                                    Sample<T> &sy, 
                                                    float (*dist) (T, T)), 
                                         float (*dist) (T, T), 
                                         unsigned K, float omega=INF, unsigned step=1) {
    float t1;
    list<int>::iterator it;

    vector< pair<unsigned, float> > knn(K);

    ClockReset();

    for (unsigned i=0; i<seqs.samples->size(); i+=step) {

        ClockStop();
        list<int> P(seqs.samples->size()-1);
        it = P.begin();
        for (unsigned j=0; j<seqs.samples->size(); j++) {
            if (i == j) continue;
            *it = j;
            it++;
        }

        float tbefore = ClockTotal();
        ClockContinue();
        
        k(P, i, seqs, f, dist, knn, omega);
        
        ClockStop();
        float tafter = ClockTotal();
        cout << i << " ";
        for (unsigned k=0; k<K; k++)
            cout << knn[k].first << " ";
        cout << "bcsample:" << (unsigned long)(seqs.samples->size());
        cout << " " << "timesample:" << tafter-tbefore;
        cout << " cost:" << knn[0].second;
        for (unsigned k=1; k<K; k++)
            cout << "," << knn[k].second;
        cout << endl;
        ClockContinue();
    }

    t1 = ClockTotal();

    cout << "time: " << t1 << endl;

    cout << "bound_count: " 
         << (unsigned long)(seqs.samples->size()/step)*seqs.samples->size() 
         << "/" 
         << (unsigned long)(seqs.samples->size()/step)*seqs.samples->size() << endl;

    cout << "local_distances: "
         << ( (unsigned long long)((seqs.samples->size()/step) * seqs.samples->size()) * 
              (unsigned long long)((*seqs.samples)[0].size() * (*seqs.samples)[0].size()) )
         << "/" 
         << ( (unsigned long long)((seqs.samples->size()/step) * seqs.samples->size()) * 
              (unsigned long long)((*seqs.samples)[0].size() * (*seqs.samples)[0].size()) ) << endl;
}

template <class T> void experiment_kReal_dats(Samples<t_value> &dats,
                                              Samples<T> &seqs,
                                              float (*f)(Sample<T> &sx,
                                                         Sample<T> &sy,
                                                         float (*dist) (T, T)),
                                              float (*dist) (T, T),
                                              unsigned K, float omega=INF) {
    int step = 1;

    float t1;
    list<int>::iterator it;

    vector< pair<unsigned, float> > knn(K);

    list<int> seqs_sub;

    for (unsigned i=0; i<dats.samples->size()-1; i++)
        for (unsigned j=0; j<(*dats.samples)[i].size(); j++)
            seqs_sub.push_back((*dats.samples)[i][j]);

    seqs_sub.sort();
    
    ClockReset();

    for (list<int>::iterator it_sub=seqs_sub.begin(); it_sub!=seqs_sub.end(); it_sub++) {

        unsigned i = *it_sub;

        ClockStop();
        list<int> P(seqs.samples->size()-1);
        it = P.begin();
        for (unsigned j=0; j<seqs.samples->size(); j++) {
            if (i == j) continue;
            *it = j;
            it++;
        }

        float tbefore = ClockTotal();
        ClockContinue();
        
        k(P, i, seqs, f, dist, knn, omega);
        
        ClockStop();
        float tafter = ClockTotal();
        cout << i << " ";
        for (unsigned k=0; k<K; k++)
            cout << knn[k].first << " ";
        cout << "bcsample:" << (unsigned long)(seqs.samples->size());
        cout << " " << "timesample:" << tafter-tbefore;
        cout << " cost:" << knn[0].second;
        for (unsigned k=1; k<K; k++)
            cout << "," << knn[k].second;
        cout << endl;
        ClockContinue();
    }

    t1 = ClockTotal();

    cout << "time: " << t1 << endl;

    cout << "bound_count: " 
         << (unsigned long)(seqs_sub.size()/step)*seqs.samples->size() 
         << "/" 
         << (unsigned long)(seqs_sub.size()/step)*seqs.samples->size() << endl;

    cout << "local_distances: "
         << ( (unsigned long long)((seqs_sub.size()/step) * seqs.samples->size()) * 
              (unsigned long long)((*seqs.samples)[0].size() * (*seqs.samples)[0].size()) )
         << "/" 
         << ( (unsigned long long)((seqs_sub.size()/step) * seqs.samples->size()) * 
              (unsigned long long)((*seqs.samples)[0].size() * (*seqs.samples)[0].size()) ) << endl;
}


template <class T> void experiment_kAesa(Samples<T> &seqs,
                                         float (*f)(Sample<T> &sx, 
                                                    Sample<T> &sy, 
                                                    float (*dist) (T, T)),
                                         float (*dist) (T, T),
                                         string &filename_tt,
                                         unsigned K,
                                         void (*fkaesa) (set<int> &P, unsigned x, 
                                                         matrix<float> &D, Samples<T> &seqs, 
                                                         float (*f)(Sample<T> &sx, Sample<T> &sy, 
                                                                    float (*dist) (T, T)), 
                                                         float (*dist) (T, T),
                                                         vector< pair<unsigned, float> > &knn, 
                                                         unsigned long &bound_count,
                                                         float H, float omega),
                                         float H=0, float omega=INF, unsigned step=1) {
    float t1 = 0;
    list<int>::iterator it;

    Tt tt;
    tt.read(filename_tt);
    matrix<float> &D = (*tt.distances);
    // unsigned lenD = seqs.samples->size();
    // matrix<float> D(lenD, lenD);
    // for (unsigned i=0; i<lenD; i++)
    //     for (unsigned j=0; j<lenD; j++)
    //         D[i][j] = 1.;

    vector< pair<unsigned, float> > knn(K);

    ClockReset();

    // unsigned num = 0;
    
    unsigned long sum_bound_count = 0;
    for (unsigned i=0; i<seqs.samples->size(); i+=step) {

        ClockStop();
        set<int> P;
        for (unsigned j=0; j<seqs.samples->size(); j++) {
            if (i == j) continue;
            P.insert(j);
        }

        float tbefore = ClockTotal();
        ClockContinue();
        
        unsigned long bound_count = 0;
        fkaesa(P, i, D, seqs, f, dist, knn, bound_count, H, omega);
        sum_bound_count += bound_count;
        
        ClockStop();
        float tafter = ClockTotal();
        cout << i << " ";
        for (unsigned k=0; k<K; k++)
            cout << knn[k].first << " ";
        cout << "bcsample:" << bound_count;
        cout << " " << "timesample:" << tafter-tbefore;
        cout << " cost:" << knn[0].second;
        for (unsigned k=1; k<K; k++)
            cout << "," << knn[k].second;
        cout << endl;
        ClockContinue();

        // if (num > 10) break;
        // num++;
    }

    t1 = ClockTotal();
    cout << "time: " << t1 << endl;
    cout << "bound_count: " << sum_bound_count << "/" 
         << (unsigned long)(seqs.samples->size()/step)*seqs.samples->size() << endl;
    cout << "local_distances: " 
         << ( (unsigned long long)(sum_bound_count) * 
              (unsigned long long)((*seqs.samples)[0].size() * (*seqs.samples)[0].size()) )
         << "/" 
         << ( (unsigned long long)((seqs.samples->size()/step) * seqs.samples->size()) * 
              (unsigned long long)((*seqs.samples)[0].size() * (*seqs.samples)[0].size()) ) << endl;
}

template <class T> void experiment_kAesa_dats(Samples<t_value> &dats,
                                              Samples<T> &seqs,
                                              float (*f)(Sample<T> &sx, 
                                                         Sample<T> &sy, 
                                                         float (*dist) (T, T)),
                                              float (*dist) (T, T),
                                              string &filename_tt,
                                              unsigned K,
                                              void (*fkaesa) (set<int> &P, unsigned x, 
                                                              matrix<float> &D, Samples<T> &seqs, 
                                                              float (*f)(Sample<T> &sx, Sample<T> &sy, 
                                                                         float (*dist) (T, T)), 
                                                              float (*dist) (T, T),
                                                              vector< pair<unsigned, float> > &knn, 
                                                              unsigned long &bound_count,
                                                              float H, float omega),
                                              float H=0, float omega=INF) {
    int step = 1;
    
    float t1 = 0;
    list<int>::iterator it;

    Tt tt;
    tt.read(filename_tt);
    matrix<float> &D = (*tt.distances);
    // unsigned lenD = seqs.samples->size();
    // matrix<float> D(lenD, lenD);
    // for (unsigned i=0; i<lenD; i++)
    //     for (unsigned j=0; j<lenD; j++)
    //         D[i][j] = 1.;

    vector< pair<unsigned, float> > knn(K);

    list<int> seqs_sub;

    for (unsigned i=0; i<dats.samples->size()-1; i++) // don't consider default class
        for (unsigned j=0; j<(*dats.samples)[i].size(); j++)
            seqs_sub.push_back((*dats.samples)[i][j]);

    seqs_sub.sort();

    ClockReset();

    // unsigned num = 0;
    
    unsigned long sum_bound_count = 0;
    for (list<int>::iterator it_sub=seqs_sub.begin(); it_sub!=seqs_sub.end(); it_sub++) {
        unsigned i = *it_sub;

        ClockStop();
        set<int> P;
        for (unsigned j=0; j<seqs.samples->size(); j++) {
            if (i == j) continue;
            P.insert(j);
        }
        
        float tbefore = ClockTotal();
        ClockContinue();
        
        unsigned long bound_count = 0;
        fkaesa(P, i, D, seqs, f, dist, knn, bound_count, H, omega);
        sum_bound_count += bound_count;
        
        ClockStop();
        float tafter = ClockTotal();
        cout << i << " ";
        for (unsigned k=0; k<K; k++)
            cout << knn[k].first << " ";
        cout << "bcsample:" << bound_count;
        cout << " " << "timesample:" << tafter-tbefore;
        cout << " cost:" << knn[0].second;
        for (unsigned k=1; k<K; k++)
            cout << "," << knn[k].second;
        cout << endl;
        ClockContinue();

        // if (num > 10) break;
        // num++;
    }

    t1 = ClockTotal();
    cout << "time: " << t1 << endl;
    cout << "bound_count: " << sum_bound_count << "/"
         << (unsigned long)(seqs_sub.size()/step)*seqs.samples->size() << endl;
    cout << "local_distances: " 
         << ( (unsigned long long)(sum_bound_count) * 
              (unsigned long long)((*seqs.samples)[0].size() * (*seqs.samples)[0].size()) )
         << "/" 
         << ( (unsigned long long)((seqs_sub.size()/step) * seqs.samples->size()) * 
              (unsigned long long)((*seqs.samples)[0].size() * (*seqs.samples)[0].size()) ) << endl;
}

template <class T> void experiment_klAesa(Samples<T> &seqs,
                                          float (*f)(Sample<T> &sx, 
                                                     Sample<T> &sy, 
                                                     float (*dist) (T, T)),
                                          float (*dist) (T, T),
                                          string &filename_tt,
                                          unsigned K,
                                          void (*fklaesa) (set<int> &P, set<int> &B, unsigned x, 
                                                           matrix<float> &D, Samples<T> &seqs, 
                                                           float (*f)(Sample<T> &sx, Sample<T> &sy, 
                                                                      float (*dist) (T, T)),
                                                           float (*dist) (T, T),
                                                           vector< pair<unsigned, float> > &knn, 
                                                           unsigned long &bound_count, 
                                                           unsigned long long &ldistances_count,
                                                           float H, float omega),
                                          unsigned pivots, float H=0, float omega=INF, unsigned step=1) {

    float t1 = 0;
    list<int>::iterator it;

    Tt tt;
    tt.read(filename_tt);
    matrix<float> &D = (*tt.distances);

    vector< pair<unsigned, float> > knn(K);

    ClockReset();

    if (pivots > seqs.samples->size()-1) {
        pivots = seqs.samples->size()-100;
        cerr << "Warning: number of pivots corrected to " << pivots << "." << endl;
    }

    unsigned long sum_bound_count = 0;
    unsigned long long sum_ldistances_count = 0;
    for (unsigned i=0; i<seqs.samples->size(); i+=step) {

        ClockStop();

        set<int> P;
        for (unsigned j=0; j<seqs.samples->size(); j++) {
            if (i == j) continue;
            P.insert(j);
        }

        // SMD {
        set<int> B;
        set<int> PB;
        for (set<int>::iterator it=P.begin(); it!=P.end(); ++it)
            PB.insert(*it);
        unsigned b_prime = *(PB.begin());
        B.insert(b_prime);
        vector<float> A(PB.size()+1);
        for (unsigned k=0; k<A.size(); k++) 
            A[k] = 0;
        PB.erase(PB.begin());
        while (B.size() < pivots) {
            float max_b = 0;
            unsigned b = b_prime;
            for (set<int>::iterator it=PB.begin(); it!=PB.end(); ++it) {
                int p = *it;
                A[p] += D[b][p];
                if (A[p] > max_b) {
                    b_prime = p;
                    max_b = A[p];
                }
            }
            B.insert(b_prime);
            PB.erase(b_prime);
        }
        // }

        float tbefore = ClockTotal();
        ClockContinue();
        
        unsigned long bound_count = 0;
        unsigned long long ldistances_count = 0;
        fklaesa(P, B, i, D, seqs, f, dist, knn, bound_count, ldistances_count, H, omega);
        sum_bound_count += bound_count;
        sum_ldistances_count += ldistances_count;
        
        ClockStop();
        float tafter = ClockTotal();
        cout << i << " ";
        for (unsigned k=0; k<K; k++)
            cout << knn[k].first << " ";
        cout << "bcsample:" << bound_count;
        cout << " " << "timesample:" << tafter-tbefore;
        cout << " cost:" << knn[0].second;
        for (unsigned k=1; k<K; k++)
            cout << "," << knn[k].second;
        cout << endl;
        ClockContinue();

    }

    t1 = ClockTotal();
    cout << "time: " << t1 << endl;
    cout << "bound_count: " << sum_bound_count << "/" 
         << (unsigned long)(seqs.samples->size()/step)*seqs.samples->size() << endl;
    cout << "local_distances: " << sum_ldistances_count << "/" 
         << ( (unsigned long long)((seqs.samples->size()/step) * seqs.samples->size()) * 
              (unsigned long long)((*seqs.samples)[0].size() * (*seqs.samples)[0].size()) ) << endl;
}

template <class T> void experiment_klAesa_dats(Samples<t_value> &dats,
                                               Samples<T> &seqs,
                                               float (*f)(Sample<T> &sx, 
                                                          Sample<T> &sy, 
                                                          float (*dist) (T, T)),
                                               float (*dist) (T, T),
                                               string &filename_tt,
                                               unsigned K,
                                               void (*fklaesa) (set<int> &P, set<int> &B, unsigned x, 
                                                                matrix<float> &D, Samples<T> &seqs, 
                                                                float (*f)(Sample<T> &sx, Sample<T> &sy, 
                                                                           float (*dist) (T, T)),
                                                                float (*dist) (T, T),
                                                                vector< pair<unsigned, float> > &knn, 
                                                                unsigned long &bound_count, 
                                                                unsigned long long &ldistances_count,
                                                                float H, float omega),
                                               unsigned pivots, float H=0, float omega=INF) {

    int step = 1;

    float t1 = 0;
    list<int>::iterator it;

    Tt tt;
    tt.read(filename_tt);
    matrix<float> &D = (*tt.distances);

    vector< pair<unsigned, float> > knn(K);

    list<int> seqs_sub;

    for (unsigned i=0; i<dats.samples->size()-1; i++) // don't consider default class
        for (unsigned j=0; j<(*dats.samples)[i].size(); j++)
            seqs_sub.push_back((*dats.samples)[i][j]);

    seqs_sub.sort();

    ClockReset();

    if (pivots > seqs.samples->size()-1) {
        pivots = seqs.samples->size()-100;
        cerr << "Warning: number of pivots corrected to " << pivots << "." << endl;
    }

    unsigned long sum_bound_count = 0;
    unsigned long long sum_ldistances_count = 0;
    for (list<int>::iterator it_sub=seqs_sub.begin(); it_sub!=seqs_sub.end(); it_sub++) {
        unsigned i = *it_sub;

        ClockStop();

        set<int> P;
        for (unsigned j=0; j<seqs.samples->size(); j++) {
            if (i == j) continue;
            P.insert(j);
        }

        // SMD {
        set<int> B;
        set<int> PB;
        for (set<int>::iterator it=P.begin(); it!=P.end(); ++it)
            PB.insert(*it);
        unsigned b_prime = *(PB.begin());
        B.insert(b_prime);
        vector<float> A(PB.size()+1);
        for (unsigned k=0; k<A.size(); k++) 
            A[k] = 0;
        PB.erase(PB.begin());
        while (B.size() < pivots) {
            float max_b = 0;
            unsigned b = b_prime;
            for (set<int>::iterator it=PB.begin(); it!=PB.end(); ++it) {
                int p = *it;
                A[p] += D[b][p];
                if (A[p] > max_b) {
                    b_prime = p;
                    max_b = A[p];
                }
            }
            B.insert(b_prime);
            PB.erase(b_prime);
        }
        // }

        float tbefore = ClockTotal();
        ClockContinue();
        
        unsigned long bound_count = 0;
        unsigned long long ldistances_count = 0;
        fklaesa(P, B, i, D, seqs, f, dist, knn, bound_count, ldistances_count, H, omega);
        sum_bound_count += bound_count;
        sum_ldistances_count += ldistances_count;
        
        ClockStop();
        float tafter = ClockTotal();
        cout << i << " ";
        for (unsigned k=0; k<K; k++)
            cout << knn[k].first << " ";
        cout << "bcsample:" << bound_count;
        cout << " " << "timesample:" << tafter-tbefore;
        cout << " cost:" << knn[0].second;
        for (unsigned k=1; k<K; k++)
            cout << "," << knn[k].second;
        cout << endl;
        ClockContinue();

    }

    t1 = ClockTotal();
    cout << "time: " << t1 << endl;
    cout << "bound_count: " << sum_bound_count << "/" 
         << (unsigned long)(seqs_sub.size()/step)*seqs.samples->size() << endl;
    cout << "local_distances: " << sum_ldistances_count << "/" 
         << ( (unsigned long long)((seqs_sub.size()/step) * seqs.samples->size()) * 
              (unsigned long long)((*seqs.samples)[0].size() * (*seqs.samples)[0].size()) ) << endl;
}

void usage(char *program) {
    cout << "Usage:\n";
    cout << program << " [Options] (-f filename)\n";
    cout << "Options:\n";
    cout << "-e experiment (kreal, kbounded, kaesa -filename_tt required-)\n";
    cout << "-u subexperiment  (kreal: kreal;\n";
    cout << "                   kaesa: kaesa; \n";
    cout << "                   klaesa: klaesa1, kblaesa1, kblaesa1bbed, kblaesa1bbed1bubu; \n";
    cout << "                   )\n";
    cout << "-m method \n";
    cout << "-l local_distance \n";
    cout << "-d type \n";
    cout << "-s size\n";
    cout << "-k k-nearest-neighbours\n";
    cout << "-t filename_tt\n";
    cout << "-v metadata.dats\n";
    cout << "-p step\n";
    cout << "-h looseness\n";
    cout << "-o omega (Sakoe)\n";
    cout << "-b filename_tt_ub\n";
    cout << "-i number of pivots (klaesa)\n";

    cout << endl;
}

int main(int argc, char **argv) {

    extern char *optarg;
    int c;

    string experiment = "";
    string subexperiment = "";
    string method = "";
    string local_distance = "";
    string type = "";
    unsigned vector_size = 0;
    string datatype = "";
    string filename = "";
    string filename_tt = "";
    string filename_dats = "";
    unsigned K = 0;
    unsigned step = 1;
    float H = 0;
    float omega = INF;
    string filename_tt_ub = "";
    unsigned pivots = 0;

    while ((c = getopt(argc, argv, "e:u:m:l:d:f:t:v:s:k:p:h:o:b:i:")) != EOF)
        switch (c) {
          case 'e': {
              experiment = optarg;
              break;
          }
          case 'u': {
              subexperiment = optarg;
              break;
          }
          case 'm': {
              method = optarg;
              break;
          }
          case 'l': {
              local_distance = optarg;
              break;
          }
          case 'd': {
              datatype = optarg;
              break;
          }
          case 'f': {
              filename = optarg;
              break;
          }
          case 't': {
              filename_tt = optarg;
              break;
          }
          case 'v': {
              filename_dats = optarg;
              break;
          }
          case 's': {
              vector_size = atoi(optarg);
              break;
          }
          case 'k': {
              K = atoi(optarg);
              break;
          }
          case 'p': {
              step = atoi(optarg);
              break;
          }
          case 'h': {
              H = myatof(optarg);
              break;
          }
          case 'o': {
              omega = myatof(optarg);
              break;
          }
          case 'b': {
              filename_tt_ub = optarg;
              break;
          }
          case 'i': {
              pivots = atoi(optarg);
              break;
          }
        }

    if (filename == "") {
        usage(argv[0]);
        exit(1);
    }

    if (experiment == "") {
        usage(argv[0]);
        exit(1);
    }

    if (subexperiment == "") {
        usage(argv[0]);
        exit(1);
    }

    if (datatype == "") {
        usage(argv[0]);
        exit(1);
    }

    if (experiment == "kaesa" and filename_tt == "") {
        usage(argv[0]);
        exit(1);
    }

    if (experiment == "klaesa" and (filename_tt == "" or pivots == 0)) {
        usage(argv[0]);
        exit(1);
    }

    if (vector_size == 0) {
        usage(argv[0]);
        exit(1);
    }

    if (K == 0)
        K = 1;

    Samples<t_value> dats;
    if (filename_dats != "") {
        if (step != 1)
            cerr << "Warning: there's no step option with dats file." << endl;
        dats.read(filename_dats.c_str());
    }

    if (datatype == "floats") {
        Samples<t_value> seqs;
    
        seqs.read(filename.c_str());

        if (experiment == "kreal")
            if (method == "ed")
                experiment_kReal(seqs, ed, levenshtein, K, omega, step);
            else if (method == "BBEd")
                experiment_kReal(seqs, BBEd, levenshtein, K, omega, step);
            else {
                    cerr << "Error: method is not suitable for this experiment" << endl;
                    exit(1);
            }

        else if (experiment == "kaesa") {
            if (subexperiment == "kaesa")

                if (method == "ed")
                    experiment_kAesa(seqs, ed, levenshtein, filename_tt, K, kAesa, 
                                     H, omega, step);
                else if (method == "BBEd")
                    experiment_kAesa(seqs, BBEd, levenshtein, filename_tt, K, kAesa, 
                                 H, omega, step);
                else {
                    cerr << "Error: method is not suitable for this experiment" << endl;
                    exit(1);
                }

            else {
                cerr << "Error: subexperiment is not suitable for this experiment" << endl;
                exit(1);
            }
        }
        else if (experiment == "klaesa") {
            if (subexperiment == "klaesa1") {
                if (method == "ed")
                    experiment_klAesa(seqs, ed, levenshtein, filename_tt, K, klAesa1, pivots,
                                      H, omega, step);
                else if (method == "BBEd")
                    experiment_klAesa(seqs, BBEd, levenshtein, filename_tt, K, klAesa1, pivots, 
                                      H, omega, step);
                else {
                    cerr << "Error: method is not suitable for this experiment" << endl;
                    exit(1);
                }
            }
            else if (subexperiment == "kblaesa1BBEdExt")
                experiment_klAesa(seqs, BBEd, levenshtein, filename_tt, K, kblAesa1BBEdExt, pivots, 
                                  H, omega, step);
            else if (subexperiment == "kblaesa1BBEdExt1bubu")
                experiment_klAesa(seqs, BBEd, levenshtein, filename_tt, K, kblAesa1BBEdExt1bubu, pivots, 
                                  H, omega, step);
            else if (subexperiment == "kblaesa1maesEdExt1bubu")
                experiment_klAesa(seqs, maesEd, levenshtein, filename_tt, K, kblAesa1maesEdExt1bubu, pivots, 
                                  H, omega, step);
            else if (subexperiment == "kblaesa1EdExt")
                experiment_klAesa(seqs, ed, levenshtein, filename_tt, K, kblAesa1EdExt, pivots, 
                                  H, omega, step);
            else if (subexperiment == "kblaesa1EdExt2")
                experiment_klAesa(seqs, ed, levenshtein, filename_tt, K, kblAesa1EdExt2, pivots, 
                                  H, omega, step);
            else {
                cerr << "Error: subexperiment is not suitable for this experiment" << endl;
                exit(1);
            }
        }
    }
    else {
        cerr << "Error: don't know datatype" << endl;
        exit(1);
    }
    return 0;
}
