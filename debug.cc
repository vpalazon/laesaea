/*===================================================================================

  File: debug.cc

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

#include <vector>
#include <list>
#include <iostream>

#include "debug.hh"

using namespace std;

void showVector(vector<float> &v) {
    vector<float>::iterator it;
    cout << "[";
    for (it = v.begin(); it != v.end(); it++) {
        cout << *it;
        if (it != v.end()-1) cout << ",";
    }
    cout << "]" << endl;
}

void showVector(vector<int> &v) {
    vector<int>::iterator it;
    cout << "[";
    for (it = v.begin(); it != v.end(); it++) {
        cout << *it;
        if (it != v.end()-1) cout << ",";
    }
    cout << "]" << endl;
}

void showVectorSpaced(vector<float> &v) {
    vector<float>::iterator it;
    cout << "[ ";
    for (it = v.begin(); it != v.end(); it++) {
        cout << *it;
        if (it != v.end()-1) cout << ", ";
    }
    cout << " ]" << endl;
}

void showMatrix(vector< vector<int> > &m) {
    vector< vector<int> >::iterator it1;
    vector<int>::iterator it2;

    for (it1 = m.begin(); it1 != m.end(); it1++) {
        for (it2 = (*it1).begin(); it2 != (*it1).end(); it2++) {
            cout << *it2;
            if (it2 != (*it1).end()-1) cout << ", ";
        }
        cout << endl;
    }
}

void showMatrix(vector< vector<float> > &m) {
    vector< vector<float> >::iterator it1;
    vector<float>::iterator it2;

    for (it1 = m.begin(); it1 != m.end(); it1++) {
        for (it2 = (*it1).begin(); it2 != (*it1).end(); it2++) {
            cout << *it2;
            if (it2 != (*it1).end()-1) cout << ", ";
        }
        cout << endl;
    }
}

void showList(list<float> &l) {
    list<float>::iterator it;
    cout << "[";
    for (it = l.begin(); it != l.end(); it++) {
        cout << *it;
        if (it != --l.end()) cout << ",";
    }
    cout << "]" << endl;
}

void showList(list<int> &l) {
    list<int>::iterator it;
    cout << "[";
    for (it = l.begin(); it != l.end(); it++) {
        cout << *it;
        if (it != --l.end()) cout << ",";
    }
    cout << "]" << endl;
}

void showSet(set<int> &l) {
    set<int>::iterator it;
    cout << "[";
    for (it = l.begin(); it != l.end(); it++) {
        cout << *it;
        if (it != --l.end()) cout << ",";
    }
    cout << "]" << endl;
}
