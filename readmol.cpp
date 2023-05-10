#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include "readmol.h"


ostream &operator<<(ostream & os, const Atom & atom) {
    os << atom.name << " " << atom.type << " " << atom.x << " " << atom.y << " " << atom.z;
    return os;
}

// output bond
ostream &operator<<(ostream & os, const Bond & bond) {
    os << bond.atom1 << " " << bond.atom2 << " " << bond.type;
    return os;
}

// output molecule
ostream &operator<<(ostream & os, const Molecule & mol) {
    os << mol.name << endl;
    for (auto atom : mol.atoms) {
        os << atom << endl;
    }
    return os;
}

//

Molecule readxyz(const string filename){
    Molecule mol;
    ifstream file(filename);
    string line;
    getline(file, line);
    int natom = stoi(line);
    getline(file, line);
    mol.name = line;
    for (int i = 0; i < natom; i++) {
        getline(file, line);
        istringstream iss(line);
        Atom atom;
        iss >> atom.name >> atom.x >> atom.y >> atom.z;
        mol.atoms.push_back(atom);
    }

    return mol;
}

vector<vector<double> > calculated_distance_matrix(Molecule mol){
    vector<vector<double> > distance_matrix;
    for (int i = 0; i < mol.atoms.size(); i++) {
        vector<double> row;
        for (int j = 0; j < mol.atoms.size(); j++) {
            double x1 = mol.atoms[i].x;
            double y1 = mol.atoms[i].y;
            double z1 = mol.atoms[i].z;
            double x2 = mol.atoms[j].x;
            double y2 = mol.atoms[j].y;
            double z2 = mol.atoms[j].z;
            double distance = sqrt(pow(x1-x2, 2) + pow(y1-y2, 2) + pow(z1-z2, 2));
            row.push_back(distance);
        }
        distance_matrix.push_back(row);
    }
    return distance_matrix;
}

void gen_bond(Molecule& mol){
    vector<vector<double> > distance_matrix = calculated_distance_matrix(mol);
    for (int i = 0; i < distance_matrix.size(); i++) {
        for (int j = 0; j < distance_matrix.size(); j++) {
            if (distance_matrix[i][j] < 1.7 && distance_matrix[i][j] > 0.8) {
                Bond bond;
                bond.atom1 = i;
                bond.atom2 = j;
                bond.type = "1";
                mol.bonds.push_back(bond);
            }
            if (distance_matrix[i][j] < 1.0 && distance_matrix[i][j] > 0.8) {
                Bond bond;
                bond.atom1 = i;
                bond.atom2 = j;
                bond.type = "2";
                mol.bonds.push_back(bond);
            }
        }
    }
}

#ifdef TEST
int main(int argc, const char * argv[]) {
    Molecule mol = readxyz("/Users/zhouoh/Downloads/scenes/CH4.xyz");
    cout << mol << endl;
    vector<vector<double> > distance_matrix = calculated_distance_matrix(mol);
    gen_bond(mol);
    for (auto bond : mol.bonds) {
        cout << bond << endl;
    }
    return 0;
}
#endif