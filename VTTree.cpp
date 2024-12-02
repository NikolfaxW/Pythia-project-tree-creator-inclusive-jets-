//
// Created by nikol on 11/28/2024.
//
#include "VTTree.h"

#include <mutex>
#include <iostream>

#include "TFile.h"


VTTree::VTTree(std::string seed, const unsigned int id) {
    _id = id;
    _seed = seed;
    T = new TTree("T", "Stores information about inclusive jets");
    T->Branch("l11", &_l11, "l11");
    T->Branch("l105", &_l105, "l105");
    T->Branch("l115", &_l115, "l115");
    T->Branch("l12", &_l12, "l12");
    T->Branch("l13", &_l13, "l13");
    T->Branch("l20", &_l20, "l20");
    T->Branch("jet_pT", &Jet_PT, "jet_pT");
    T->Branch("eta", &rapidity, "eta");

}

bool VTTree::saveTTree(std::string path) {
    std::lock_guard<std::mutex> lock(std::mutex);
    TFile *file = new TFile((path + std::to_string(_id) + " : " + _seed + ".root").c_str(),
                            "RECREATE"); //create a file to store the tree as an output
    if (file->IsZombie())
        return false;

    T->Print();
    T->Write();

    file->Close();
    delete file;
    return true;
}

VTTree::~VTTree() {
    saveTTree("../results/");
    delete T;
}