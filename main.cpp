#include <iostream>
#include "tinyxml2.h"
#include <string> 
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include "readmol.h"
#include <map>

using namespace std;
using namespace tinyxml2; 

const float bond_radius = 0.16;

const float double_bond_radius = 0.08;
const float double_bond_shift = 0.16;

const float triple_bond_radius = 0.04;
const float triple_bond_shift = 0.1;

const float dash_bond_radius = 0.08;

map<string,float> atom_radius = {
    {"H", 0.3},
    {"C", 0.5},
    {"N", 0.5},
    {"O", 0.5},
    {"F", 0.5},
    {"P", 1.06},
    {"S", 1.02},
    {"Cl", 0.99},
    {"Br", 1.14},
    {"I", 1.33}
};

map<string,vector<float> > atom_color = {
    {"H", {0.7, 0.7, 0.7}},
    {"C", {0.2, 0.2, 0.2}},
    {"N", {0.188, 0.313, 0.973}},
    {"O", {0.9, 0.05, 0.05}},
    {"F", {0.565, 0.878, 0.314}},
    {"P", {1.0, 0.502, 0.0}},
    {"S", {1.0, 1.0, 0.188}},
    {"Cl", {0.121, 0.941, 0.121}},
    {"Br", {0.651, 0.161, 0.161}},
    {"I", {0.580, 0.0, 0.580}}
};

float vector_angle(vector<float> v1, vector<float> v2){
    float dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
    float v1_norm = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
    float v2_norm = sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);
    float cos = dot/(v1_norm*v2_norm);
    float angle = acos(cos);
    return angle;
}

vector<float> vector_add(vector<float> v1, vector<float> v2){
    float x = v1[0] + v2[0];
    float y = v1[1] + v2[1];
    float z = v1[2] + v2[2];
    vector<float> v;
    v.push_back(x);
    v.push_back(y);
    v.push_back(z);
    return v;
}

vector<float> cross_product(vector<float> v1, vector<float> v2){
    vector<float> v;
    v.push_back(v1[1]*v2[2] - v1[2]*v2[1]);
    v.push_back(v1[2]*v2[0] - v1[0]*v2[2]);
    v.push_back(v1[0]*v2[1] - v1[1]*v2[0]);
    return v;
}

vector<float> normalize(vector<float> v, float length){
    float norm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    vector<float> v_norm;
    v_norm.push_back(length * v[0]/norm);
    v_norm.push_back(length * v[1]/norm);
    v_norm.push_back(length * v[2]/norm);
    return v_norm;
}

vector<float> calculate_bond_vector(Molecule mol, int atom1, Atom atom2){
    vector<float> v;
    v.push_back(atom2.x - mol.atoms[atom1].x);
    v.push_back(atom2.y - mol.atoms[atom1].y);
    v.push_back(atom2.z - mol.atoms[atom1].z);
    return v;
}

vector<float> calculate_bond_vector(Molecule mol, int atom1, int atom2){
    vector<float> v;
    v.push_back(mol.atoms[atom2].x - mol.atoms[atom1].x);
    v.push_back(mol.atoms[atom2].y - mol.atoms[atom1].y);
    v.push_back(mol.atoms[atom2].z - mol.atoms[atom1].z);
    return v;
}

vector<float> calculate_shift_vector(Molecule mol, int atom1, int atom2, float bond_shift){
    vector<float> v;
    Atom ref_atom;
    int ref_atom_index;
    vector<int> atom1_bonded_atoms;
    vector<int> atom2_bonded_atoms;
    vector<Atom> atom1_bonded_atoms_vector;
    vector<Atom> atom2_bonded_atoms_vector;
    bool found = false;
    //check all the bonds of atom1
    for (auto bond : mol.bonds){
        if (bond.atom1 == atom1){
            if (bond.atom2 == atom2){
                continue;
            }
            else{
                atom1_bonded_atoms.push_back(bond.atom2);
                atom1_bonded_atoms_vector.push_back(mol.atoms[bond.atom2]);
                found = true;
            }
        }
        else if (bond.atom2 == atom1){
            if (bond.atom1 == atom2){
                continue;
            }
            else{
                atom1_bonded_atoms.push_back(bond.atom1);
                atom1_bonded_atoms_vector.push_back(mol.atoms[bond.atom1]);
                found = true;
            }
        }
    }
    // if atom1 is not bonded to atom2, check all the bonds of atom2
    if (!found){
        for (auto bond : mol.bonds){
            if (bond.atom1 == atom2){
                if (bond.atom2 == atom1){
                    continue;
                }
                else{
                    atom2_bonded_atoms.push_back(bond.atom2);
                    atom2_bonded_atoms_vector.push_back(mol.atoms[bond.atom2]);
                    found = true;
                }
            }
            else if (bond.atom2 == atom2){
                if (bond.atom1 == atom1){
                    continue;
                }
                else{
                    atom2_bonded_atoms.push_back(bond.atom1);
                    atom2_bonded_atoms_vector.push_back(mol.atoms[bond.atom1]);
                    found = true;
                }
            }
        }
    }

    v = calculate_bond_vector(mol, atom1, atom2);
    if (!found){
        vector<float> shift = cross_product(v, {0.0f, 0.0f, -10.0f});
        vector<float> shift_norm = normalize(shift, 0.1);
        return shift_norm;
    }

    if (atom1_bonded_atoms.size() <= atom2_bonded_atoms.size()){
        vector<float> ref_v = {0.0f, 0.0f, 0.0f};
        vector<float> bv;
        for (auto atom : atom1_bonded_atoms_vector){
            if (atom.x == mol.atoms[atom2].x && atom.y == mol.atoms[atom2].y && atom.z == mol.atoms[atom2].z){
                continue;
            }
            bv = calculate_bond_vector(mol, atom1, atom);
            ref_v = vector_add(ref_v, bv);
        }
    }
    else{
        vector<float> ref_v = {0.0f, 0.0f, 0.0f};
        vector<float> bv;
        for (auto atom : atom2_bonded_atoms_vector){
            if (atom.x == mol.atoms[atom1].x && atom.y == mol.atoms[atom1].y && atom.z == mol.atoms[atom1].z){
                continue;
            }
            bv = calculate_bond_vector(mol, atom2, atom);
            ref_v = vector_add(ref_v, {atom.x, atom.y, atom.z});
        }
    }

    vector<float> ref_v = calculate_bond_vector(mol, atom1, ref_atom_index);
    vector<float> shift = cross_product(v, ref_v);
    vector<float> real_shift = cross_product(v, shift);
    vector<float> shift_norm = normalize(real_shift, bond_shift);
    return shift_norm;

}



void draw_atoms(Atom atom, XMLDocument &doc, XMLElement *scence){
    
        auto shape = doc.NewElement("shape");
        scence->InsertEndChild(shape);
        shape->SetAttribute("type", "sphere");
        auto float2 = doc.NewElement("float");
        shape->InsertEndChild(float2);
        float2->SetAttribute("name", "radius");
        float radius = atom_radius.at(atom.name);
        float2->SetAttribute("value", radius);
        auto transform = doc.NewElement("transform");
        shape->InsertEndChild(transform);
        transform->SetAttribute("name", "toWorld");
        auto translate = doc.NewElement("translate");
        transform->InsertEndChild(translate);
        translate->SetAttribute("x", atom.x);
        translate->SetAttribute("y", atom.y);
        translate->SetAttribute("z", atom.z);
        auto ref = doc.NewElement("ref");
        shape->InsertEndChild(ref);
        ref->SetAttribute("id", atom.name.c_str());
    
}

void draw_bond(Atom atom1, Atom atom2, XMLDocument &doc, XMLElement *scence){

        auto shape = doc.NewElement("shape");
        scence->InsertEndChild(shape);
        shape->SetAttribute("type", "cylinder");
        auto float2 = doc.NewElement("float");
        shape->InsertEndChild(float2);
        float2->SetAttribute("name", "radius");
        float2->SetAttribute("value", bond_radius);

        auto point = doc.NewElement("point");
        shape->InsertEndChild(point);
        point->SetAttribute("name", "p0");
        string p0 = to_string(atom1.x) + ", " + to_string(atom1.y) + ", " + to_string(atom1.z);
        point->SetAttribute("value", p0.c_str());

        auto point1 = doc.NewElement("point");
        shape->InsertEndChild(point1);
        point1->SetAttribute("name", "p1");
        string p1 = to_string(atom2.x) + ", " + to_string(atom2.y) + ", " + to_string(atom2.z);
        point1->SetAttribute("value", p1.c_str());

        auto transform = doc.NewElement("transform");

        shape->InsertEndChild(transform);
        transform->SetAttribute("name", "toWorld");
        auto ref = doc.NewElement("ref");
        shape->InsertEndChild(ref);
        ref->SetAttribute("id", "ball");
        
}

void draw_double_bond(Atom atom1, Atom atom2, vector<float> shift_vector,  XMLDocument &doc, XMLElement *scence){
    for (int i = 0; i < 2; i++){
        auto shape = doc.NewElement("shape");
        scence->InsertEndChild(shape);
        shape->SetAttribute("type", "cylinder");
        auto float2 = doc.NewElement("float");
        shape->InsertEndChild(float2);
        float2->SetAttribute("name", "radius");
        float2->SetAttribute("value", double_bond_radius);

        auto point = doc.NewElement("point");
        shape->InsertEndChild(point);
        point->SetAttribute("name", "p0");
        float x = atom1.x + (pow(-1,i)) * shift_vector[0];
        float y = atom1.y + (pow(-1,i)) * shift_vector[1];
        float z = atom1.z + (pow(-1,i)) * shift_vector[2];
        string p0 = to_string(x) + ", " + to_string(y) + ", " + to_string(z);
        point->SetAttribute("value", p0.c_str());

        auto point1 = doc.NewElement("point");
        shape->InsertEndChild(point1);
        point1->SetAttribute("name", "p1");
        x = atom2.x + (pow(-1,i)) * shift_vector[0];
        y = atom2.y + (pow(-1,i)) * shift_vector[1];
        z = atom2.z + (pow(-1,i)) * shift_vector[2];
        string p1 = to_string(x) + ", " + to_string(y) + ", " + to_string(z);
        point1->SetAttribute("value", p1.c_str());

        auto transform = doc.NewElement("transform");

        shape->InsertEndChild(transform);
        transform->SetAttribute("name", "toWorld");
        auto ref = doc.NewElement("ref");
        shape->InsertEndChild(ref);
        ref->SetAttribute("id", "ball");
    };

}

void draw_triple_bond(Atom atom1, Atom atom2, vector<float> shift_vector,  XMLDocument &doc, XMLElement *scence){
    for (int i = -1; i < 2; i++){
        auto shape = doc.NewElement("shape");
        scence->InsertEndChild(shape);
        shape->SetAttribute("type", "cylinder");
        auto float2 = doc.NewElement("float");
        shape->InsertEndChild(float2);
        float2->SetAttribute("name", "radius");
        float2->SetAttribute("value", triple_bond_radius);

        auto point = doc.NewElement("point");
        shape->InsertEndChild(point);
        point->SetAttribute("name", "p0");
        float x = atom1.x + i * shift_vector[0];
        float y = atom1.y + i * shift_vector[1];
        float z = atom1.z + i * shift_vector[2];
        string p0 = to_string(x) + ", " + to_string(y) + ", " + to_string(z);
        point->SetAttribute("value", p0.c_str());

        auto point1 = doc.NewElement("point");
        shape->InsertEndChild(point1);
        point1->SetAttribute("name", "p1");
        x = atom2.x + i * shift_vector[0];
        y = atom2.y + i * shift_vector[1];
        z = atom2.z + i * shift_vector[2];
        string p1 = to_string(x) + ", " + to_string(y) + ", " + to_string(z);
        point1->SetAttribute("value", p1.c_str());

        auto transform = doc.NewElement("transform");

        shape->InsertEndChild(transform);
        transform->SetAttribute("name", "toWorld");
        auto ref = doc.NewElement("ref");
        shape->InsertEndChild(ref);
        ref->SetAttribute("id", "ball");
    };

}

void draw_dash_bond(Atom atom1, Atom atom2, float dash_length, float interval_length ,  XMLDocument &doc, XMLElement *scence){
    float distance = sqrt(pow(atom1.x - atom2.x, 2) + pow(atom1.y - atom2.y, 2) + pow(atom1.z - atom2.z, 2));
    int num = distance / (dash_length + interval_length);
    vector<float> bond_vector = {atom2.x - atom1.x, atom2.y - atom1.y, atom2.z - atom1.z};

    for (int i=0; i<num; i++){
        auto shape = doc.NewElement("shape");
        scence->InsertEndChild(shape);
        shape->SetAttribute("type", "cylinder");
        auto float2 = doc.NewElement("float");
        shape->InsertEndChild(float2);
        float2->SetAttribute("name", "radius");
        float2->SetAttribute("value", dash_bond_radius);

        auto point = doc.NewElement("point");
        shape->InsertEndChild(point);
        point->SetAttribute("name", "p0");
        float x = atom1.x + i * (dash_length + interval_length) * bond_vector[0] / distance;
        float y = atom1.y + i * (dash_length + interval_length) * bond_vector[1] / distance;
        float z = atom1.z + i * (dash_length + interval_length) * bond_vector[2] / distance;
        string p0 = to_string(x) + ", " + to_string(y) + ", " + to_string(z);
        point->SetAttribute("value", p0.c_str());

        auto point1 = doc.NewElement("point");
        shape->InsertEndChild(point1);
        point1->SetAttribute("name", "p1");
        x = atom1.x + (i + 1) * dash_length * bond_vector[0] / distance;
        y = atom1.y + (i + 1) * dash_length * bond_vector[1] / distance;
        z = atom1.z + (i + 1) * dash_length * bond_vector[2] / distance;
        string p1 = to_string(x) + ", " + to_string(y) + ", " + to_string(z);
        point1->SetAttribute("value", p1.c_str());

        auto transform = doc.NewElement("transform");

        shape->InsertEndChild(transform);
        transform->SetAttribute("name", "toWorld");
        auto ref = doc.NewElement("ref");
        shape->InsertEndChild(ref);
        ref->SetAttribute("id", "ball");
    };
}
    




void write_mol_xml(Molecule mol){
    XMLDocument doc;
    auto decl_elec = doc.NewDeclaration();    //不传参会添加默认XML头信息
    doc.InsertFirstChild(decl_elec);

    auto scence = doc.NewElement("scene");
    doc.InsertAfterChild(decl_elec, scence);
    scence->SetAttribute("version", "0.5.0");
    
    auto integrator = doc.NewElement("integrator");
    scence->InsertEndChild(integrator);
    integrator->SetAttribute("type", "path");
    auto integer = doc.NewElement("integer");
    integrator->InsertEndChild(integer);
    integer->SetAttribute("name", "maxDepth");
    integer->SetAttribute("value", 5);

    auto sensor = doc.NewElement("sensor");
    scence->InsertEndChild(sensor);
    sensor->SetAttribute("type", "perspective");
    auto float1 = doc.NewElement("float");
    sensor->InsertEndChild(float1);
    float1->SetAttribute("name", "fov");
    float1->SetAttribute("value", 60);
    auto transform = doc.NewElement("transform");
    sensor->InsertEndChild(transform);
    transform->SetAttribute("name", "toWorld");
    auto lookat = doc.NewElement("lookat");
    transform->InsertEndChild(lookat);
    lookat->SetAttribute("origin", "0, 0, 15");
    lookat->SetAttribute("target", "0, 0, 0");

    auto bsdf = doc.NewElement("bsdf");
    scence->InsertEndChild(bsdf);
    bsdf->SetAttribute("type", "diffuse");
    bsdf->SetAttribute("id", "ball");

    for (auto atom_type: atom_color){
        auto bsdf = doc.NewElement("bsdf");
        scence->InsertEndChild(bsdf);
        bsdf->SetAttribute("type", "diffuse");
        bsdf->SetAttribute("id", atom_type.first.c_str());
        auto rgb = doc.NewElement("rgb");
        bsdf->InsertEndChild(rgb);
        rgb->SetAttribute("name", "reflectance");
        string color = to_string(atom_type.second[0]) + ", " + to_string(atom_type.second[1]) + ", " + to_string(atom_type.second[2]);
        rgb->SetAttribute("value", color.c_str());
    }

    auto rgb = doc.NewElement("rgb");
    bsdf->InsertEndChild(rgb);
    rgb->SetAttribute("name", "reflectance");
    rgb->SetAttribute("value", "0.5, 0.5, 0.5");

    for (auto atom : mol.atoms){
    draw_atoms(atom, doc, scence);
    }
    //draw bond
    for (auto bond : mol.bonds){
        if (bond.type == "1"){
            draw_bond(mol.atoms[bond.atom1], mol.atoms[bond.atom2], doc, scence);
        }
        if (bond.type == "2"){
            cout << "double bond" << endl;
            cout << bond.atom1 << " " << bond.atom2 << endl;
            vector<float> shift_vector = calculate_shift_vector(mol, bond.atom1, bond.atom2,double_bond_shift);
            cout << shift_vector[0] << " " << shift_vector[1] << " " << shift_vector[2] << endl;

            draw_double_bond(mol.atoms[bond.atom1], mol.atoms[bond.atom2], shift_vector, doc, scence);
        }


    }

    auto emitter = doc.NewElement("emitter");
    scence->InsertEndChild(emitter);
    emitter->SetAttribute("type", "constant");
    auto e_rgb = doc.NewElement("rgb");
    emitter->InsertEndChild(e_rgb);
    e_rgb->SetAttribute("name", "radiance");
    e_rgb->SetAttribute("value", "1.0");

    #ifdef DEBUG
    auto shape = doc.NewElement("shape");
    scence->InsertEndChild(shape);
    shape->SetAttribute("type", "rectangle");
    auto transform1 = doc.NewElement("transform");
    shape->InsertEndChild(transform1);
    transform1->SetAttribute("name", "toWorld");
    auto translate1 = doc.NewElement("translate");
    transform1->InsertEndChild(translate1);
    translate1->SetAttribute("x", 0);
    translate1->SetAttribute("y", 0);
    translate1->SetAttribute("z", -2);
    auto scale = doc.NewElement("scale");
    transform1->InsertEndChild(scale);
    scale->SetAttribute("x", 10);
    scale->SetAttribute("y", 10);
    scale->SetAttribute("z", 10);
    auto ref1 = doc.NewElement("ref");
    shape->InsertEndChild(ref1);
    ref1->SetAttribute("id", "ball");
    auto emitter = doc.NewElement("emitter");
    shape->InsertEndChild(emitter);
    emitter->SetAttribute("type", "area");
    auto rgb1 = doc.NewElement("rgb");
    emitter->InsertEndChild(rgb1);
    rgb1->SetAttribute("name", "radiance");
    rgb1->SetAttribute("value", "10, 10, 10");
    #endif

    doc.SaveFile("/Users/zhouoh/Downloads/scenes/mol.xml");


}


int main(int argc, const char * argv[]) {
    

    Molecule mol = readxyz("/Users/zhouoh/Downloads/scenes/CH4.xyz");
    gen_bond(mol);
    
    write_mol_xml(mol);
    return EXIT_SUCCESS;
}

