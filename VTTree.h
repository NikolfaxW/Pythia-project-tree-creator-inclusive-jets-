//
// Created by nikol on 11/28/2024.
//

#ifndef PYTHIA_PROJECT_TREE_CREATOR_INCLUSIVE_JETS_VTTREE_H
#define PYTHIA_PROJECT_TREE_CREATOR_INCLUSIVE_JETS_VTTREE_H


#include "TTree.h"

class VTTree {
public:
    TTree * T;
    Float_t _l11, _l105, _l115, _l12, _l13, _l20, Jet_PT, rapidity;
    std::string _seed;
    unsigned int _id;
VTTree(const std::string seed, const unsigned int id);
bool saveTTree(std::string path);
~VTTree();


};


#endif //PYTHIA_PROJECT_TREE_CREATOR_INCLUSIVE_JETS_VTTREE_H
