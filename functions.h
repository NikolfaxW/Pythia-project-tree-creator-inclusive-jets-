//
// Created by nikol on 11/14/2024.
//

#ifndef PYTHIA_PROJECT_TREE_CREATOR_FUNCTIONS_H
#define PYTHIA_PROJECT_TREE_CREATOR_FUNCTIONS_H

#include "Rtypes.h"
#include "TTree.h"
#include "VTTree.h"



void showProgressBar(int progress, int total);
Double_t delta_R(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2);
int getId(); //Ok
std::string getStatus(); //Ok
bool increaseIdOrChageStatus(int id, std::string status); //ok status can be "true" or "false" !!!
void mainSec(const int numThreads, std::string  seed, TTree *&T, Float_t &Jet_Pt, Float_t &rapidity,
             const unsigned int & requiredNumberOfJets, unsigned int & numberOfJetsFound,
             Float_t &l11, Float_t &l105, Float_t &l115, Float_t &l12, Float_t &l13, Float_t &l20);
void mainSecTest(const int numThreads, const std::string seed, const unsigned int &requiredNumberOfJets,
                 unsigned int &numberOfJetsFound, const unsigned int id);


#endif //PYTHIA_PROJECT_TREE_CREATOR_FUNCTIONS_H
