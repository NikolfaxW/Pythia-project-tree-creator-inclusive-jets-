//
// Created by nikol on 11/14/2024.
//

#include "functions.h"
#include "particleUnit.h"
#include "MyInfo.h"

#include <iostream>
#include <string>
#include <fstream>
#include <deque>
#include <vector>
#include <cmath>

#include "Pythia8/Pythia.h"
#include "TVector2.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "TH2F.h"


void showProgressBar(int progress, int total) {
    double ratio = static_cast<double>(progress) / total;
    std::cout << "(required number of jets found: " << progress << " | " << total << ") ";
    std::cout << int(ratio * 100.0) << "%\r";
    std::cout.flush();
}

Double_t delta_R(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2) {
    Double_t delta_eta = eta1 - eta2;
    Double_t delta_phi = TVector2::Phi_mpi_pi(phi1 - phi2);
    return TMath::Sqrt(pow(delta_eta, 2) + pow(delta_phi, 2));
}

int getId() {
    std::string fileName = "../status.txt";
    std::ifstream file(fileName);
    if (!file.is_open()) {
        std::cerr << "Error: file not opened" << std::endl;
        return -1;
    }
    int id;
    std::string temp;
    file >> temp >> id;
    file.close();
    return id;
}

std::string getStatus() {
    std::string fileName = "../status.txt";
    std::ifstream file(fileName);
    if (!file.is_open()) {
        std::cerr << "Error: file not opened" << std::endl;
        return "Error";
    }
    std::string temp;
    file >> temp;
    file.close();
    return temp;
}

bool increaseIdOrChageStatus(int id, std::string status) {
    std::string fileName = "../status.txt";
    std::fstream file(fileName, std::ios::out);
    if (!file.is_open()) {
        std::cerr << "Error: file not opened" << std::endl;
        return false;
    }
    file << status << std::endl << (id) << std::endl;
    file.close();
    return true;
}

void mainSec(const int numThreads, std::string  seed, TTree *&T, Float_t &Jet_Pt, Float_t &rapidity,
             const unsigned int & requiredNumberOfJets, unsigned int & numberOfJetsFound,
             Float_t &l11, Float_t &l105, Float_t &l115, Float_t &l12, Float_t &l13, Float_t &l20) {


    Pythia8::Pythia pythia; //create pythia object


    {
        std::lock_guard<std::mutex> lock(std::mutex);  // Lock the mutex
        pythia.readFile("../config1.cmnd"); //read config file and intialize pythia
        pythia.readString("Random:setSeed = on");  // Enable setting of the seed
        pythia.readString("Random:seed = " + seed);
        pythia.init();

    }


    std::list<particleUnit> j_constituents;
    std::map<TString, fastjet::JetDefinition> jetDefs; //map to store jet definitions

    //parameters for jet finding and to make program less hardcoded
    double R = 0.4; //jet radius
    double pTmin_jet = 5, pThadron = 0.2; //minimum pT for jets and hadrons

    double pTMinTrig = 0; //minimum pT for the particle to be found
    double pTMaxTrig = 5.0;
    double mTemp; //This variable are needed to recount momentum after particle mass resets
    Pythia8::Vec4 pTemp; //This variable are needed to recount momentum after particle mass resets


    Double_t R_frac, pT_frac, temp_l11, temp_l105, temp_l115, temp_l12, temp_l13, temp_l20;

    //define jet finding algorithms here:
    jetDefs["#it{k_{t}} jets, #it{R} = " + std::to_string(R)] = fastjet::JetDefinition(
            fastjet::kt_algorithm, R, fastjet::E_scheme, fastjet::Best);


    auto &event = pythia.event; //create a reference to the Pythia8 collision event
    std::vector<fastjet::PseudoJet> fjInputs; //to store particles to be clustered
    std::vector<fastjet::PseudoJet> selectedJets; //to store jets after all cuts



    for (int iEvent = 0; numberOfJetsFound < requiredNumberOfJets; ++iEvent) { //loop over needed number of events
        if (!pythia.next()) continue; //generate next event, if it is not possible, continue

        fjInputs.clear(); //clears the vector of particles to be clustered
        for (int i = 0; i < event.size(); ++i) { // saves particles in order to make jets
            auto &p = event[i];
            if (event[i].isFinal()) { //checks if particle is final state
                pTemp = p.p();
                if (p.idAbs() == 22) //changes mass for any other particle except photons to pion mass
                    mTemp = 0;
                else
                    mTemp = 0.13957;
                fastjet::PseudoJet particleTemp = fastjet::PseudoJet(p.px(), p.py(), p.pz(),
                                                                     sqrt(pTemp.pAbs2() + mTemp *
                                                                                          mTemp)); //recounts 4 momentum and energy after mass reset
                particleTemp.set_user_info(
                        new MyInfo(p.idAbs(), i, event[i].isCharged())); //adds additional info to the particle
                fjInputs.push_back(particleTemp); //saves the particle to the vector to be clustered

            }
        }


        if (fjInputs.empty()) { // Abort to avoid running algorithms through empty vectors
            std::cout << "Error: event with no final state particles" << std::endl;
            continue;
        }


        for (auto jetDef: jetDefs) { //for each jet definition runs jet clustering sequence, then saves it to ROOT tree
            selectedJets.clear(); //empties the vector of jets - prevents from saving jets from previous event
            fastjet::ClusterSequence clusterSequence(fjInputs, jetDef.second); //sets up the cluster sequence
            auto jets = sorted_by_pt(clusterSequence.inclusive_jets(0)); //runs the clustering and sorts jets by pT
            fastjet::Selector eta_selector = fastjet::SelectorEtaRange(-1 + R, 1 - R); // Selects jets with |eta| < 1-R
            selectedJets = eta_selector(jets); //applies the eta selector to the jets

            for (const auto &jet: selectedJets) { //loop through all jets

                temp_l11 = 0;
                temp_l105 = 0;
                temp_l115 = 0;
                temp_l12 = 0;
                temp_l13 = 0;
                temp_l20 = 0;


                for (const auto &c: jet.constituents()) { //loop through all jet constituents to calculate l11, l105, l115, l12, l13, l20
                    if (!c.user_info<MyInfo>().isCharged())
                        continue; // neutral particles
                    pT_frac = c.pt() / jet.pt();
                    R_frac = delta_R(jet.eta(), jet.phi(), c.eta(), c.phi()) / R;
                    temp_l11 += pT_frac * R_frac;
                    temp_l105 += pT_frac * pow(R_frac, 0.5);
                    temp_l115 = pT_frac * pow(R_frac, 1.5);
                    temp_l12 += pT_frac * pow(R_frac, 2);
                    temp_l13 += pT_frac * pow(R_frac, 3);
                    temp_l20 += pow(pT_frac, 2);
                }


                {
                    //! Not safe to fill the tree from multiple threads
                    std::lock_guard<std::mutex> lock(std::mutex);  // Lock the mutex
                    l11 = temp_l11;//Often is more than 1
                    l105 = temp_l105;//Often is more than 1
                    l115 = temp_l115; //Rarely is more than one
                    l12 = temp_l12; //Rarely is more than one
                    l13 = temp_l13;//Often is more than 1
                    l20 = temp_l20;//Rarely is more than 1
                    // Jet_Pt = jet.pt();
                    rapidity = jet.rapidity();
                    T->Fill();  // Fill the data vector safely
                    ++numberOfJetsFound;
                    if (numberOfJetsFound % 100 == 0) showProgressBar(numberOfJetsFound, requiredNumberOfJets);
                }
            }

        }
    }
}

void mainSecTest(const int numThreads, const std::string seed, const unsigned int &requiredNumberOfJets, const unsigned int id) {
    unsigned int numberOfJetsFound = 0;
    Pythia8::Pythia pythia; //create pythia object
    {
        std::lock_guard<std::mutex> lock(std::mutex);  // Lock the mutex
        pythia.readFile("../config1.cmnd"); //read config file and intialize pythia
        pythia.readString("Random:setSeed = on");  // Enable setting of the seed
        pythia.readString("Random:seed = " + seed);
        pythia.init();
    }

    std::list<particleUnit> j_constituents;
    std::map<TString, fastjet::JetDefinition> jetDefs; //map to store jet definitions

    //parameters for jet finding and to make program less hardcoded
    double R = 0.4; //jet radius
    double pTmin_jet = 5, pThadron = 0.2; //minimum pT for jets and hadrons
    double pTMinTrig = 0; //minimum pT for the particle to be found
    double pTMaxTrig = 5.0;
    double mTemp; //This variable are needed to recount momentum after particle mass resets
    Pythia8::Vec4 pTemp; //This variable are needed to recount momentum after particle mass resets
    VTTree T(seed, id);

    //define jet finding algorithms here:
    jetDefs["#it{k_{t}} jets, #it{R} = " + std::to_string(R)] = fastjet::JetDefinition(
            fastjet::kt_algorithm, R, fastjet::E_scheme, fastjet::Best);

    auto &event = pythia.event; //create a reference to the Pythia8 collision event
    std::vector<fastjet::PseudoJet> fjInputs; //to store particles to be clustered
    std::vector<fastjet::PseudoJet> selectedJets; //to store jets after all cuts
    Float_t pT_frac, R_frac;
    for (int iEvent = 0; numberOfJetsFound < requiredNumberOfJets; ++iEvent) { //loop over needed number of events
        if (!pythia.next()) continue; //generate next event, if it is not possible, continue
        fjInputs.clear(); //clears the vector of particles to be clustered
        for (int i = 0; i < event.size(); ++i) { // saves particles in order to make jets
            auto &p = event[i];
            if (event[i].isFinal()) { //checks if particle is final state
                pTemp = p.p();
                if (p.idAbs() == 22) //changes mass for any other particle except photons to pion mass
                    mTemp = 0;
                else
                    mTemp = 0.13957;
                fastjet::PseudoJet particleTemp = fastjet::PseudoJet(p.px(), p.py(), p.pz(),
                                                                     sqrt(pTemp.pAbs2() + mTemp *
                                                                                          mTemp)); //recounts 4 momentum and energy after mass reset
                particleTemp.set_user_info(
                        new MyInfo(p.idAbs(), i, event[i].isCharged())); //adds additional info to the particle
                fjInputs.push_back(particleTemp); //saves the particle to the vector to be clustered

            }
        }
        if (fjInputs.empty()) { // Abort to avoid running algorithms through empty vectors
            std::cout << "Error: event with no final state particles" << std::endl;
            continue;
        }

        for (auto jetDef: jetDefs) { //for each jet definition runs jet clustering sequence, then saves it to ROOT tree
            selectedJets.clear(); //empties the vector of jets - prevents from saving jets from previous event
            fastjet::ClusterSequence clusterSequence(fjInputs, jetDef.second); //sets up the cluster sequence
            auto jets = sorted_by_pt(clusterSequence.inclusive_jets(0)); //runs the clustering and sorts jets by pT
            fastjet::Selector eta_selector = fastjet::SelectorEtaRange(-1 + R, 1 - R); // Selects jets with |eta| < 1-R
            selectedJets = eta_selector(jets); //applies the eta selector to the jets

            for (const auto &jet: selectedJets) { //loop through all jets
                T._l11 = 0;
                T._l105 = 0;
                T._l115 = 0;
                T._l12 = 0;
                T._l13 = 0;
                T._l20 = 0;
                for (const auto &c: jet.constituents()) { //loop through all jet constituents to calculate l11, l105, l115, l12, l13, l20
                    if (!c.user_info<MyInfo>().isCharged())
                        continue; // neutral particles
                    pT_frac = c.pt() / jet.pt();
                    R_frac = delta_R(jet.eta(), jet.phi(), c.eta(), c.phi()) / R;
                    T._l11 += pT_frac * R_frac;
                    T._l105 += pT_frac * pow(R_frac, 0.5);
                    T._l115 = pT_frac * pow(R_frac, 1.5);
                    T._l12 += pT_frac * pow(R_frac, 2);
                    T._l13 += pT_frac * pow(R_frac, 3);
                    T._l20 += pow(pT_frac, 2);
                }
                T.T->Fill();
                numberOfJetsFound++;
                std::cout << numberOfJetsFound << std::endl;
//                    std::lock_guard<std::mutex> lock(std::mutex);
//                    ++numberOfJetsFound;
//                    if (numberOfJetsFound % 100 == 0) showProgressBar(numberOfJetsFound, requiredNumberOfJets);
//                }

            }
        }

    }
}
