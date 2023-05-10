using namespace std;
struct Atom {
    string name;
    string type;
    double x;
    double y;
    double z;
};

struct Bond {
    int atom1;
    int atom2;
    string type;
};

struct Molecule {
    string name;
    vector<Atom> atoms;
    vector<Bond> bonds;
};
ostream &operator<<(ostream & os, const Atom & atom);
ostream &operator<<(ostream & os, const Bond & bond);
ostream &operator<<(ostream & os, const Molecule & mol);
Molecule readxyz(const string filename);
vector<vector<double> > calculated_distance_matrix(Molecule mol);
void gen_bond(Molecule& mol);


