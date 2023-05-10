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

const map<string,float> atom_radius = {
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

const map<string,vector<float>> atom_color = {
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
    //draw bond
    for (auto bond : mol.bonds){
        auto shape = doc.NewElement("shape");
        scence->InsertEndChild(shape);
        shape->SetAttribute("type", "cylinder");
        auto float2 = doc.NewElement("float");
        shape->InsertEndChild(float2);
        float2->SetAttribute("name", "radius");
        float2->SetAttribute("value", 0.1);

        auto point = doc.NewElement("point");
        shape->InsertEndChild(point);
        point->SetAttribute("name", "p0");
        string p0 = to_string(mol.atoms[bond.atom1].x) + ", " + to_string(mol.atoms[bond.atom1].y) + ", " + to_string(mol.atoms[bond.atom1].z);
        point->SetAttribute("value", p0.c_str());

        auto point1 = doc.NewElement("point");
        shape->InsertEndChild(point1);
        point1->SetAttribute("name", "p1");
        string p1 = to_string(mol.atoms[bond.atom2].x) + ", " + to_string(mol.atoms[bond.atom2].y) + ", " + to_string(mol.atoms[bond.atom2].z);
        point1->SetAttribute("value", p1.c_str());


        auto transform = doc.NewElement("transform");

        shape->InsertEndChild(transform);
        transform->SetAttribute("name", "toWorld");
        auto ref = doc.NewElement("ref");
        shape->InsertEndChild(ref);
        ref->SetAttribute("id", "ball");
    }
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

    doc.SaveFile("/mnt/e/mol.xml");


}


int main(int argc, const char * argv[]) {
    

    Molecule mol = readxyz("/home/zhouoh/mol2xml/CH4.xyz");
    gen_bond(mol);
    write_mol_xml(mol);
    return EXIT_SUCCESS;
}


