#include <iostream>
#include <fstream>
#include <string>
// ROOT headers
#include "TROOT.h"
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TObject.h>
#include <vector>
// Delphes headers
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"
// fastjet
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/PseudoJet.hh"

 float jetpt1;
 float jetpt2;
 
 float jeteta1;
 float jeteta2;
 
 float mjj;
 float mbb;
 
 //---- h>bb
 float hbb_pt;
 float hbb_eta;
 float hbb_phi;
 float hbb_e;
 float hbb_mass;

 float gen_hbb_pt;
 float gen_hbb_phi;
 float gen_hbb_eta;
 float gen_hbb_e;
 float gen_hbb_mass;

 float bjetpt1;
 float bjetpt2;

 float bjeteta1;
 float bjeteta2;

 float bjetphi1;
 float bjetphi2;

 float bjete1;
 float bjete2;
 
 //---- h>WW
 float hww_mt;
 float hww_pt;
 float hww_phi;
 float hww_etap; //---- ambiguity on the sign
 float hww_etam;
 
 float gen_hww_mt;
 float gen_hww_pt;
 float gen_hww_phi;
 float gen_hww_eta;

 float hw1_pt;
 float hw2_pt;
 float hw1_eta;
 float hw2_eta;
 float hw1_phi;
 float hw2_phi;
 float hw1_e;
 float hw2_e;
 
 float pt1;
 float pt2;
 float nlep;
 float channel;
 float mll;
 float ptll;
 float pzll;
 float dphill;
 float pfmet;
 
 float gen_pt1;
 float gen_pt2;
 float gen_nlep;
 float gen_channel;
 float gen_mll;
 float gen_dphill;
 float gen_ptll;
 float gen_pzll;
 float gen_pfmet;
 float gen_pfmez;
 float gen_mvv;


 //---- x>hh_m (ww)
 float xhh_ww_mt;
 
 float xhh_m_ww_pt;
 float xhh_m_ww_eta;
 float xhh_m_ww_phi;
 float xhh_m_ww_m;
 float xhh_p_ww_pt;
 float xhh_p_ww_eta;
 float xhh_p_ww_phi;
 float xhh_p_ww_m;
 
 GenParticle *particle;
 Electron *electron;
 Photon *photon;
 Muon *muon;
 MissingET *met;
 Track *track;
 Tower *tower;
 TObject *object;
 TLorentzVector momentum;
 Float_t Eem, Ehad;
 Bool_t skip;
 Long64_t entry;
 Int_t i, j, pdgCode;
