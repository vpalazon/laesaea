/*===================================================================================

  File: Tt.cc

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

#include <string>
#include "Tt.hh"
#include "Useful.hh"

using namespace std;

Tt::~Tt() {
    delete this->distances;
}

void Tt::read(string &filename) {

    char buf[MAX_LINE_SZ_TT];
    string aux;
    string svalue;
    Sample<float> s;
    char value[MAX_VALUE_SZ_TT];
    unsigned current_i=0, last_i=0;
    unsigned i, j;

    string fn(filename);
    if (not FileExists(fn)) {
        cerr << "Error: File does not exist.\n" << endl;
        exit(1);
    }

    fstream file_op0(filename.c_str(), ios::in);

    while(!file_op0.eof()) {
        file_op0.getline(buf, MAX_LINE_SZ_TT);
        svalue = buf;
        if (svalue == "") continue;
        istrstream isvalue(svalue.c_str(), svalue.size());
        while(isvalue.get(value, MAX_VALUE_SZ_TT, ' ')) {
            isvalue.get();
            float aux;
            aux = convertToFloat(std::string(value));
            s.values.push_back(aux);
        }
        current_i = s.values[0];
        s.values.erase(s.values.begin(), s.values.end());
    }

    file_op0.close();

    last_i = current_i;

    fstream file_op(filename.c_str(), ios::in);

    matrix<float> *distances = new matrix<float>(last_i+1, last_i+1);

    while(!file_op.eof()) {
        file_op.getline(buf, MAX_LINE_SZ_TT);
        svalue = buf;
        if (svalue == "") continue;
        istrstream isvalue(svalue.c_str(), svalue.size());
        while(isvalue.get(value, MAX_VALUE_SZ_TT, ' ')) {
            isvalue.get();
            float aux;
            aux = convertToFloat(std::string(value));
            s.values.push_back(aux);
        }
        i = int(s.values[0]);
        j = int(s.values[1]);
        (*distances)[i][j] = (*distances)[j][i] = s.values[2];
        s.values.erase(s.values.begin(), s.values.end());
    }
    file_op.close();
    this->distances = distances;
}

void Tt::read(string &filename, int size) {

    char buf[MAX_LINE_SZ_TT];
    string aux;
    string svalue;
    Sample<float> s;
    char value[MAX_VALUE_SZ_TT];
    unsigned i, j;
    unsigned current_i = -1;

    string fn(filename);
    if (not FileExists(fn)) {
        cerr << "Error: File does not exist.\n" << endl;
        exit(1);
    }

    fstream file_op(filename.c_str(), ios::in);

    matrix<float> *distances = new matrix<float>(size, size);

    while(!file_op.eof()) {
        file_op.getline(buf, MAX_LINE_SZ_TT);
        svalue = buf;
        if (svalue == "") continue;
        istrstream isvalue(svalue.c_str(), svalue.size());
        while(isvalue.get(value, MAX_VALUE_SZ_TT, ' ')) {
            isvalue.get();
            float aux;
            aux = convertToFloat(std::string(value));
            s.values.push_back(aux);
        }
        i = int(s.values[0]);
        j = int(s.values[1]);
        if (i != current_i) {
            current_i = i;
            cout << current_i << endl;
        }
        (*distances)[i][j] = (*distances)[j][i] = s.values[2];
        s.values.erase(s.values.begin(), s.values.end());
    }
    file_op.close();
    this->distances = distances;
}

unsigned Tt::size() {
    return this->distances->size_x();
}

float Tt::get(unsigned i, unsigned j) {
    return (*(this->distances))[i][j];
}
