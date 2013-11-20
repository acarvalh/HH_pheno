/*******************************************************************************/
/** development:
    Andrea Massironi 
    Olivier Boundu
    Alexandra Carvalho                                                        **/
/*******************************************************************************/
/*  to run
 make ntupleProducerVbfHHbbXX.exe
 ./ntupleProducerVbfHHbbXX.exe -t test -i GluGluToHHTo2B2G_M-125_8TeV_madgraph_v2_DEL_v03.root -o testout.root
******************************************************************************/
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
// fastjet headers
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/MassDropTagger.hh"
// basic parameters
#include "baseline.h"
// analysis cuts
#include "cuts.dat"
// Verbosity
#define DEBUG 0
#include <boost/program_options.hpp>
// ntuple variables
#include "ntupleProducerVbfHHbbXX.h"
using namespace std;
//using namespace fastjet;
typedef fastjet::JetDefinition::DefaultRecombiner DefRecomb;

struct myclassMin {
 bool operator() (std::pair<float, std::pair <int,Int_t> > i, std::pair<float, std::pair <int,Int_t> > j) { 
  return (i.first < j.first);
 }
} myObjMin;


struct myclassMax {
 bool operator() (std::pair<float, std::pair <int,Int_t> > i, std::pair<float, std::pair <int,Int_t> > j) { 
  return (i.first > j.first);
 }
} myObjMax;


///////////////////////////////////////////////////
int fourb();
int dobranches(TTree* outtree);
bool findleptons(TClonesArray *branchMissingET ,TClonesArray *branchElectron, TClonesArray *branchMuon,
		ExRootTreeReader* treeReader, bool doHwwselection); 
int findphotons(TClonesArray *branchPhoton,ExRootTreeReader* treeReader, bool doHwwselection);
bool findjets(TClonesArray *branchJet,ExRootTreeReader* treeReader, 
		bool doHbbselection, int & countJets, int & counttags, std::vector<int> & tagentry);
bool istagged(Jet *jet);
bool isThisJetALepton(TLorentzVector* jet, TLorentzVector* l1, TLorentzVector* l2);
/////////////////////////////////////////////////////
bool analyse_4b(int countJets, int counttags, std::vector<int> tagentry);
bool jets_semi_hadronic(int countJets, int counttags, std::vector<int> tagentry);
bool analyse_2b2w(bool findlepton);
bool VBFjes(TLorentzVector* jet1,TLorentzVector* jet2);
////////////////////////////////////////////////////

namespace po = boost::program_options;

TLorentzVector pho1, pho2; // save to compare to jets
TLorentzVector l1, l2; // save to compare to jets
std::vector<TLorentzVector> Jets;// jet1, jet2, jet3, jet4, jet5, jet6; // all the jets

int main (int argc, char **argv) {
 int r = fourb();
 // declare arguments
 std::string inputfile;
 std::string outputfile;
 std::string outputtree;
 bool doHwwselection;
 bool doHggselection;
 bool doHbbselection;
 // print out passed arguments
 std::copy(argv, argv + argc, std::ostream_iterator<char*>(std::cout, " "));
 // argument parsing
 try {
  po::options_description desc("Allowed options");
  desc.add_options()
  ("help,h", "produce help message")
  ("inputfile,i", po::value<std::string>(&inputfile)->default_value("../GluGluToHHTo2B2G_M-125_8TeV_madgraph_v2_DEL_v03.root"), "input file")
  ("outputfile,o", po::value<std::string>(&outputfile)->default_value("output.root"), "output file")
  ("outputtree,t", po::value<std::string>(&outputtree)->default_value("GluGluToHHTo2B2G_8TeV"), "output tree")
  ("doHwwselection", po::value<bool>(&doHwwselection)->default_value(false),  "apply Hww selection")  
  ("doHggselection", po::value<bool>(&doHggselection)->default_value(false), "apply Hgg selection")  
  ("doHbbselection", po::value<bool>(&doHbbselection)->default_value(true), "apply Hbb selection")  
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  if (vm.count("help")) {std::cout << desc << "\n";return 1;}
  } catch(std::exception& e) {std::cerr << "error: " << e.what() << "\n";return 1;} 
    catch(...) {std::cerr << "Exception of unknown type!\n";}
 // end of argument parsing
 //################################################

 // check if at least one selection is turned on
 if( (!doHwwselection) && (!doHggselection) && (!doHbbselection) ) 
   {std::cerr << std::endl << "ERROR: Exactly one Higgs selection must be turned on, exiting" << std::endl; return 1;} 
 // check if at most one selection is turned on
 if( ((doHwwselection) && (doHggselection)) || ((doHggselection) && (doHbbselection)) || ((doHwwselection) && (doHbbselection)) ) 
  {std::cerr << std::endl << "ERROR: Exactly one Higgs selection must be turned on, exiting" << std::endl; return 1;}
 
 //---- analysis channel ----
// float  KindSelection = -1;  //---- 0 = WWbb,  1=ggbb,   2=bbbb
// if (doHwwselection) KindSelection = 0;
// if (doHggselection) KindSelection = 1;
// if (doHbbselection) KindSelection = 2;
 
 TChain *chain = new TChain("Delphes");
 chain->Add(inputfile.c_str());      
 ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
 TFile *outfile = new TFile(outputfile.c_str(), "RECREATE");
 TTree *outtree = new TTree(outputtree.c_str(), "reduced");
 
 int dobranch = dobranches(outtree);

 //---- objects in Delphes format 
// TClonesArray *branchParticle = treeReader->UseBranch("Particle");
// TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
// TClonesArray *branchEFlowTower = treeReader->UseBranch("EFlowTower");
// TClonesArray *branchEFlowMuon = treeReader->UseBranch("EFlowMuon");
 
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");

 //---- events
 Long64_t allEntries = treeReader->GetEntries();
 std::cout << "** Chain contains " << allEntries << " events" << std::endl;
 
 // Loop over all events
 for(entry = 0; entry < 10; entry++) {
  // Load selected branches with data from specified event
  treeReader->ReadEntry(entry);
  // analysis itself -- it alredy fill the variables
  // jets
  //if(!doHbbselection) {
      // EW objects first
      if(doHwwselection) { // save 2 leptons, and MET
	bool findlepton = findleptons(branchMissingET,branchElectron,branchMuon,treeReader, doHwwselection); 
        bool bbww = analyse_2b2w(findlepton); 
        std::cout<<"leptons"<<std::endl;
      } // ww sel
      if(doHggselection) { // save 2 photons
			//bool aabb = analyse_2b2a();  
	int findphoton = findphotons(branchPhoton,treeReader, doHggselection); 
        std::cout<<"photons"<<std::endl;
      } // gg sel
  // jets
  // find, tag, return all jets
  int countJets = 0, counttags=0; std::vector<int> tagentry;
  bool findjet = findjets(branchJet, treeReader, doHbbselection,countJets, counttags,tagentry); 
  //std::cout<<"jets"<<std::endl;
  // analyse
  if(!doHbbselection) {bool semi = jets_semi_hadronic(countJets, counttags,tagentry);}
  if(doHbbselection) {bool fourb = analyse_4b(countJets, counttags,tagentry);}
  //  
  outtree->Fill();
 } // close loop entry
 outfile->cd();
 outtree->Write();
 outfile->Close();
 return 0;
} // close main
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int fourb(){
  std::cout<<"hi!!!!"<<std::endl;
  //int r = 0;
  return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool analyse_4b(int countJets, int counttags, std::vector<int> tagentry){ // to Alexandra
return false;
} // close 4b analysis
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool jets_semi_hadronic(int countJets, int counttags, std::vector<int> tagentry){ // to us
  // 
/*  float mjj12 = (Jets[0]+Jets[1]).M();
  float mjj13 = (Jets[0]+Jets[2]).M();
  float mjj14 = (Jets[0]+Jets[3]).M();
  float mjj23 = (Jets[1]+Jets[2]).M();
  float mjj24 = (Jets[1]+Jets[3]).M();
  float mjj34 = (Jets[2]+Jets[3]).M();
  //
  std::map<float, int> m_maxmjj;
  m_maxmjj[-mjj12] = 1;
  m_maxmjj[-mjj13] = 2;
  m_maxmjj[-mjj14] = 3;
  m_maxmjj[-mjj23] = 4;
  m_maxmjj[-mjj24] = 5;
  m_maxmjj[-mjj34] = 6;
  //
  // write here the resolved case -- things are already TLorentzVectors 
  std::vector<PseudoJet> bJets, vbfJets; 
  std::map<float, int>::iterator it_type_m_maxmjj = m_maxmjj.begin();
  if (it_type_m_maxmjj->second == 1) { 
		vbfJets.push_back(Jets[0]); vbfJets.push_back(Jets[1]); 
		bJet1 = jet3; bJet2 = jet4; }; // and so on
*/
/*  if (it_type_m_maxmjj->second == 2) { Jet1 = jet1; Jet2 = jet3; bJet1 = jet2; bJet2 = jet4; };
  if (it_type_m_maxmjj->second == 3) { Jet1 = jet1; Jet2 = jet4; bJet1 = jet2; bJet2 = jet3; };
  if (it_type_m_maxmjj->second == 4) { Jet1 = jet2; Jet2 = jet3; bJet1 = jet1; bJet2 = jet4; };
  if (it_type_m_maxmjj->second == 5) { Jet1 = jet2; Jet2 = jet4; bJet1 = jet1; bJet2 = jet3; };
  if (it_type_m_maxmjj->second == 6) { Jet1 = jet3; Jet2 = jet4; bJet1 = jet1; bJet2 = jet2; };
*/  
  /*
  //---- sub-order in pt: jetpt1 > jetpt2
  if (Jet1.Pt() < Jet2.Pt()) {
   TLorentzVector tempjet = Jet1;
   Jet1 = Jet2;
   Jet2 = tempjet;
  }
  
  //---- sub-order in pt: bjetpt1 > bjetpt2
  if (bJet1.Pt() < bJet2.Pt()) {
   TLorentzVector tempjet = bJet1;
   bJet1 = bJet2;
   bJet2 = tempjet;
  }  
  
  
  //---- save information
  
  jetpt1 = Jet1.Pt();
  jetpt2 = Jet2.Pt();
  bjetpt1 = bJet1.Pt();
  bjetpt2 = bJet2.Pt();
  
  jeteta1 = Jet1.Eta();
  jeteta2 = Jet2.Eta();
  bjeteta1 = bJet1.Eta();
  bjeteta2 = bJet2.Eta();
  */
  
  /*
  TLorentzVector hbb;
  hbb = bJet1 + bJet2;
  
  mjj = (Jet1 +  Jet2 ).M();
  mbb = (bJet1 + bJet2).M();
  
  //---- h>bb
  hbb_pt = (bJet1 +  bJet2 ).Pt();
  hbb_eta = (bJet1 +  bJet2 ).Eta();
  hbb_phi = (bJet1 +  bJet2 ).Phi();
  */  

return false;
} // close semi hadronic
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool analyse_2b2w(bool findlepton){ // to Andrea
  //-------------
/*
  //-----------------
  //---- leptons ----
  hww_mt = -99;
  xhh_ww_mt = -99.;
  
  xhh_m_ww_pt  = -99;
  xhh_m_ww_eta = -99;
  xhh_m_ww_phi = -99;
  xhh_m_ww_m   = -99;
  
  xhh_p_ww_pt  = -99;
  xhh_p_ww_eta = -99;
  xhh_p_ww_phi = -99;
  xhh_p_ww_m   = -99;    
  
  //-------------------
  //---- hh > WWbb ----
  //---- at least 2 leptons ----
   
   TLorentzVector gen_met_vector;
   int nH = 0;
   for(int iPart = 0; iPart < branchParticle->GetEntriesFast(); iPart++) {
    GenParticle* particle = (GenParticle*) branchParticle->At(iPart);
    //---- neutrinos
    int pdgCode = TMath::Abs(particle->PID);
    int IsPU = particle->IsPU;
    int status = particle->Status;
    
    if (IsPU == 0 && status == 3 && (pdgCode == 12 || pdgCode == 14 || pdgCode == 16) ) {
     gen_met_vector = gen_met_vector + particle->P4();
     //     if (particle->M1 != -1) std::cout << " particle->M1 = " << particle->M1 << std::endl;
     //     if (particle->M2 != -1) std::cout << " particle->M2 = " << particle->M2 << std::endl;
     //     if (particle->D1 != -1) std::cout << " particle->D1 = " << particle->D1 << std::endl;
     //     if (particle->D2 != -1) std::cout << " particle->D2 = " << particle->D2 << std::endl;
     
     //    h ->  W W -> lvlv
     //     if (particle->M1 != -1) {
     //      GenParticle* possibleW = (GenParticle*) (branchParticle->At(particle->M1));    
     //      if (possibleW  && TMath::Abs(possibleW->PID) == 24 ) {
     //       GenParticle* possibleH = (GenParticle*) (branchParticle->At(possibleW->M1));
     //       if (possibleH && possibleH->PID == 25  ) {
     //        gen_met_vector = gen_met_vector + particle->P4();
     //       }
     //      }
     //     }
     
    }
    
    
    if (IsPU == 0  &&  pdgCode == 35) { //--- 35 = "modified Higgs" (the "25" one is the one decaying into 2b")
     gen_hww_pt  = particle->P4().Pt(); 
     gen_hww_phi = particle->P4().Phi(); 
     gen_hww_eta = particle->P4().Eta(); 
     nH++;
    }
    
    if (IsPU == 0  &&  pdgCode == 25) { //--- the "25" higgs is the one decaying into 2b"
     gen_hbb_pt  = particle->P4().Pt(); 
     gen_hbb_phi = particle->P4().Phi(); 
     gen_hbb_eta = particle->P4().Eta(); 
     nH++;
    }
    
    
   }
   
   //   std::cout << " nH = " << nH << std::endl;
   
   gen_pfmet = gen_met_vector.Pt();
   gen_pfmez = gen_met_vector.Pz();
   gen_mvv = gen_met_vector.M();
   
   
   
   TLorentzVector hww;
   TLorentzVector hwwp;
   TLorentzVector hwwm;
   
   if (pfmet != -99) {
    //   HiggsMass
    //---- h>ww
    TLorentzVector vmet;
    //--- IMPORTANT: h>ww, mll ~ mvv, otherwise something missing in higgs kinematic reconstruction
    vmet.SetPtEtaPhiM(met->MET, 0, met->Phi, mll);
    
    hww = l1 + l2 + vmet;
    
    hww_pt =  (l1 + l2 + vmet ).Pt();
    hww_phi = (l1 + l2 + vmet ).Phi();
    
    //--- transverse mass
    hww_mt = sqrt((l1.Pt() + l2.Pt() + vmet.Pt())*(l1.Pt() + l2.Pt() + vmet.Pt()) - hww_pt*hww_pt);
    
    //---- kinematic fit for eta
    float sintheta2 = (hww_pt*hww_pt / (hww.E() * hww.E() - HiggsMass*HiggsMass ));
    float sintheta;
    if (sintheta2 > 0) sintheta = sqrt (sintheta2);
    if (sintheta2 > 0) {
     hww_etap = - log (tan ( asin ( sintheta ) / 2. )) ;
     hww_etam = + log (tan ( asin ( sintheta ) / 2. )) ;
     
     hwwp = hww;
     //     std::cout << " hww_pt = " << hww_pt << std::endl;
     hwwp.SetPtEtaPhiM(hww_pt, hww_etap, hww_phi, HiggsMass);
     hwwm.SetPtEtaPhiM(hww_pt, hww_etam, hww_phi, HiggsMass);
     
     //---- x>hh
     TLorentzVector xhh;
     xhh = hww + hbb;
     
     TLorentzVector xhh_p;
     TLorentzVector xhh_m;
     xhh_p = hwwm + hbb;
     xhh_m = hwwp + hbb;
     
     xhh_m_ww_pt  = xhh_m.Pt();
     xhh_m_ww_eta = xhh_m.Eta();
     xhh_m_ww_phi = xhh_m.Phi();
     xhh_m_ww_m   = xhh_m.M();
     
     xhh_p_ww_pt  = xhh_p.Pt();
     xhh_p_ww_eta = xhh_p.Eta();
     xhh_p_ww_phi = xhh_p.Phi();
     xhh_p_ww_m   = xhh_p.M();
     
     
     //--- transverse mass
     xhh_ww_mt = sqrt((l1.Pt() + l2.Pt() + vmet.Pt() + bJet1.Pt() + bJet2.Pt())*(l1.Pt() + l2.Pt() + vmet.Pt() + bJet1.Pt() + bJet2.Pt()) - xhh.Pt()*xhh.Pt());
     
    }
    else {
     hww_etap = -99;
     hww_etam = -99;
    }
    
   }
   else {
    //---- h>ww
    hww_pt = -99;
    hww_etam = -99;
    hww_etap = -99;
    hww_phi = -99;
   }
   
  */

return false;
} // close wwbb analysis
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//bool aabb analyse_2b2a(){ // to Olivier
//return false;
//} // close aabb analysis
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool findleptons( TClonesArray *branchMissingET, TClonesArray *branchElectron,TClonesArray *branchMuon,ExRootTreeReader* treeReader, bool doHwwselection){

  TLorentzVector gen_l1, gen_l2;
  ///---- take the two highest pt leptons in the event (m or e)
  //    maps are ordered in ascending order
  std::map <float, int> m_maxptleptons;
  
  // Loop over all electrons in event
  for(i = 0; i < branchElectron->GetEntriesFast(); i++) {
   electron = (Electron*) branchElectron->At(i);
   double pt = electron->PT;
   m_maxptleptons[-pt] = -(i+1);
  }
  
  // Loop over all muons in event
  for(i = 0; i < branchMuon->GetEntriesFast(); i++) {
   muon = (Muon*) branchMuon->At(i);
   double pt = muon->PT;
   m_maxptleptons[-pt] = (i+1);
  }

  if (doHwwselection && m_maxptleptons.size() > 1) {
   // kind = 0/1 if m/e
   
   std::map<float, int>::iterator it_type_m_lepton = m_maxptleptons.begin();
   int flav1 = (it_type_m_lepton->second<0);  // m>0, e<0 ---> m=0, e=1
   
   it_type_m_lepton++;
   int flav2 = (it_type_m_lepton->second<0);  // m>0, e<0 ---> m=0, e=1
   
   nlep = 0;
   for(it_type_m_lepton = m_maxptleptons.begin(); it_type_m_lepton != m_maxptleptons.end(); it_type_m_lepton++) 
		if ( -(it_type_m_lepton->first) > 10) nlep++;
  
  
   it_type_m_lepton = m_maxptleptons.begin(); 
   
   if (nlep >= 1) {
    if (it_type_m_lepton->second>0) { 
     l1     =                 ((Muon*)         branchMuon->At(  it_type_m_lepton->second - 1 ))->P4();
     gen_l1 = ((GenParticle*) ((Muon*)         branchMuon->At(  it_type_m_lepton->second - 1 ))->Particle.GetObject())->P4();
    } else                            {
     l1     =                 ((Electron*) branchElectron->At(-(it_type_m_lepton->second + 1)))->P4();
     gen_l1 = ((GenParticle*) ((Electron*) branchElectron->At(-(it_type_m_lepton->second + 1)))->Particle.GetObject())->P4();
    }
    
    if (nlep >= 2) {
     it_type_m_lepton++;
     if (it_type_m_lepton->second>0) { 
      l2     =                 ((Muon*)         branchMuon->At(  it_type_m_lepton->second - 1 ))->P4();
      gen_l2 = ((GenParticle*) ((Muon*)         branchMuon->At(  it_type_m_lepton->second - 1 ))->Particle.GetObject())->P4();
     }
     else                            { 
      l2     =                 ((Electron*) branchElectron->At(-(it_type_m_lepton->second + 1)))->P4();
      gen_l2 = ((GenParticle*) ((Electron*) branchElectron->At(-(it_type_m_lepton->second + 1)))->Particle.GetObject())->P4();
     }
    }
   }
  }  

  if (doHwwselection && m_maxptleptons.size() >= 2) {
   
   // kind = 0/1 if m/e
   
   std::map<float, int>::iterator it_type_m_lepton = m_maxptleptons.begin();
   int flav1 = (it_type_m_lepton->second<0);  // m>0, e<0 ---> m=0, e=1
   pt1 = - it_type_m_lepton->first;
   it_type_m_lepton++;
   int flav2 = (it_type_m_lepton->second<0);  // m>0, e<0 ---> m=0, e=1
   pt2 = - it_type_m_lepton->first;
   
   nlep = 0;
   for(it_type_m_lepton = m_maxptleptons.begin(); it_type_m_lepton != m_maxptleptons.end(); it_type_m_lepton++) {
    if ( -(it_type_m_lepton->first) > 10) nlep++;
   }
   
   //                       ee/mm          e   m           m    e
   channel =             flav1*flav2+2*(flav1>flav2)+3*(flav1<flav2);
   
   // # mumu #    channel == 0
   // # mue #     channel == 3
   // # emu #     channel == 2
   // # ee #      channel == 1
   
   pt1 = l1.Pt();
   pt2 = l2.Pt();
   mll = (l1+l2).M();
   ptll = (l1+l2).Pt();
   pzll = (l1+l2).Pz();
   dphill = l1.DeltaPhi(l2);
   
   gen_pt1 = gen_l1.Pt();
   gen_pt2 = gen_l2.Pt();
   gen_mll = (gen_l1+gen_l2).M();
   gen_ptll = (gen_l1+gen_l2).Pt();
   gen_pzll = (gen_l1+gen_l2).Pz();
   gen_dphill = gen_l1.DeltaPhi(gen_l2);

  }
  else {
   pt1 = -99;
   pt2 = -99;
   nlep = m_maxptleptons.size();
   channel = -1;
   mll = -99;
   dphill = -99;
   
   hww_pt = -99;
   hww_etam = -99;
   hww_etap = -99;
   hww_phi = -99;
  }

  //---- met ----  
  if(branchMissingET->GetEntriesFast() > 0) {
   met   = (MissingET*) branchMissingET->At(0);
   pfmet = met->MET;
  } else pfmet = -99;

  return false;
} // end findlepton
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool findjets(TClonesArray *branchJet,ExRootTreeReader* treeReader, 
		bool doHbbselection, int & countJets, int & counttags, std::vector<int> & tagentry){ 
  Jet *jet; // P4 returns a TLorentzVector
  //---- at least 4 jets with pt>MINPTJET GeV 
  // loop in all jets
  for(i = 0; i < branchJet->GetEntriesFast(); i++) {
   jet = (Jet*) branchJet->At(i);
   TLorentzVector jetP4 = jet->P4();
   // save jets if baseline cuts
   if ( jet->PT > MINPTJET 
	&& ((!doHbbselection ) ||  doHbbselection )){ // && !isThisJetALepton(&jetP4, &l1, &l2)
		countJets++; 
		//double Px=jetP4.Px(),Py=jetP4.Py(),Pz=jetP4.Pz(),E=jetP4.E();
		Jets.push_back(jetP4);
		//std::cout<<"jet px! "<<Px<<std::endl;
	}
  // among all jets does it have sub?
  bool tagged = istagged(jet);
  if(tagged) {counttags++, tagentry.push_back(i);}
  // if is tagged keep it
//    //cout << tagged_jet.m() << endl;
//    //if  {
//    tagged_jets.push_back(tagged_jet); // mass drop found
//    tagged_jets_index.push_back(i); // identify the fat jet
  } // close for each jet
  //if (Jets[0].Pt() < MINPTJET) std::cout << "We have a problem; countJets = " << countJets 
   //				<< "; ijet = " << ijet << " and jet4.Pt() = " << jet4.Pt() << std::endl; 
  std::cout<<"number of jets "<<countJets<<std::endl;
  //if ((!doHbbselection && countJets > 2) || (doHbbselection && countJets > 3)){
  // shrink to account substructure   
  return false;
} // end findjets
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool istagged(Jet *jet){
  TObject *object;
  //Jet *jet;
  // check constituents
  std::cout<<jet->Constituents.GetEntriesFast()<<std::endl;
  vector<fastjet::PseudoJet> particles;
      for(int j = 0; j < jet->Constituents.GetEntriesFast(); ++j){
	//jet = jetentry;
	momentum.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
	//particles[j].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
        object = jet->Constituents.At(j);
//	TLorentzVector* blabla = (TLorentzVector*) ((GenParticle*) object)->P4();
	//double px = blabla->Pt();//, py = blabla.Py(), pz = blabla.Pz(),pe = blabla.E();
	std::cout<<"Entered the loop "<<std::endl;
	//particles.push_back(fastjet::PseudoJet());
	//jetP4 = jet->Constituents.At(j)->P4();
        // Check if the constituent are accessible
        if(object == 0) continue;
        if(object->IsA() == GenParticle::Class()) 
		{momentum += ((GenParticle*) object)->P4(); }
        else if(object->IsA() == Track::Class())
		{momentum += ((Track*) object)->P4(); } 
        else if(object->IsA() == Tower::Class()) 
		{momentum += ((Tower*) object)->P4(); }
        else if(object->IsA() == Muon::Class()) 
		{ momentum += ((Muon*) object)->P4(); }
        else if(object->IsA() == Electron::Class()) 
		{ momentum += ((Electron*) object)->P4(); }
        else if(object->IsA() == Photon::Class()) 
		{ momentum += ((Photon*) object)->P4(); }
	//particles.push_back(jet->Constituents.At(j));
      particles.push_back(momentum);
      std::cout<<"the hardest core! "<<momentum.Px()<<std::endl;
      } // close for jet constituents
  //fastjet::PseudoJet teste = particles[0];
  //std::cout<<"the hardest core! "<<teste<<std::endl;
  fastjet::Selector jet_selector = fastjet::SelectorPtMin(MINPTJET) && fastjet::SelectorAbsRapMax(rapmax);
  //JetDefinition akt(antikt_algorithm, jetR);
  //ClusterSequence cs_akt(particles, akt);
  //std::vector<PseudoJet> jets_akt;
  //jets_akt = sorted_by_pt(jet_selector(cs_akt.inclusive_jets()));
  // jet definition for substructure
  fastjet::JetDefinition CA10(fastjet::cambridge_algorithm, Rsb);
  // first recluster with some large CA (needed for mass-drop)
  fastjet::ClusterSequence cs_tmp(particles, CA10);
  // next get hardest jet
  std::vector<fastjet::PseudoJet> ca_jet;
  ca_jet = fastjet::sorted_by_pt(cs_tmp.inclusive_jets()); // find the cores
  //if(ca_jet[0].pt()>0)std::cout<<"the hardest core! "<<ca_jet[0].pt()<<std::endl;
  // now run mass drop tagger / compare the hardest core with the rest of the jet
//  MassDropTagger md_tagger(mu, ycut); // define the cut on mass drop
	// mu: ratio in between mass of cores, symetric splitting
//  PseudoJet tagged_jet; // save to check if survives mass drop
//  tagged_jet = md_tagger(ca_jet);
  //if(tagged_jet.pt() > 0) {std::cout<<"tag!"<<std::endl;}
  // Filter definition to improve mass resolution // after
  //Filter filter(JetDefinition(cambridge_algorithm, Rfilt), SelectorNHardest(n_subjet));
  //JetDefinition akt(antikt_algorithm, jetR);
  //ClusterSequence cs_akt(particles, akt);
  //
  // by now we will have only resolved analysis
  return false;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int findphotons(TClonesArray *branchPhoton,ExRootTreeReader* treeReader, bool doHwwselection){

return 0;
} // end findphotons
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool isThisJetALepton(TLorentzVector* jet, TLorentzVector* l1, TLorentzVector* l2){
 bool isLep = false;
 if (l1) if (jet->DeltaR(*l1) < DRmax) isLep = true;
 if (l2) if (jet->DeltaR(*l2) < DRmax) isLep = true;
 return isLep;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 int dobranches(TTree* outtree){
 // outtree->Branch("KindSelection",  &KindSelection,  "KindSelection/F");
 // the 2 vbf jets
 outtree->Branch("jetpt1",  &jetpt1,  "jetpt1/F");
 outtree->Branch("jetpt2",  &jetpt2,  "jetpt2/F");
 outtree->Branch("jeteta1",  &jeteta1,  "jeteta1/F");
 outtree->Branch("jeteta2",  &jeteta2,  "jeteta2/F");
 // both candidates invariant mass
 outtree->Branch("mjj", &mjj, "mjj/F");
 outtree->Branch("mbb", &mbb, "mbb/F");
 // H1
 outtree->Branch("hbb_pt", &hbb_pt, "hbb_pt/F");
 outtree->Branch("hbb_eta", &hbb_eta, "hbb_eta/F");
 outtree->Branch("hbb_phi", &hbb_phi, "hbb_phi/F");
 outtree->Branch("hbb_e", &hbb_e, "hbb_e/F");
 outtree->Branch("hbb_mass", &hbb_mass, "hbb_mass/F");
 //
 outtree->Branch("gen_hbb_pt", &gen_hbb_pt, "gen_hbb_pt/F");
 outtree->Branch("gen_hbb_phi", &gen_hbb_phi, "gen_hbb_phi/F");
 outtree->Branch("gen_hbb_eta", &gen_hbb_eta, "gen_hbb_eta/F");
 outtree->Branch("gen_hbb_e", &gen_hbb_e, "gen_hbb_e/F");
 outtree->Branch("gen_hbb_mass", &gen_hbb_mass, "gen_hbb_mass/F");
 // H1 constituents
 outtree->Branch("bjeteta1", &bjeteta1, "bjeteta1/F");
 outtree->Branch("bjeteta2", &bjeteta2, "bjeteta2/F");
 outtree->Branch("bjetpt1", &bjetpt1, "bjetpt1/F");
 outtree->Branch("bjetpt2", &bjetpt2, "bjetpt2/F");
 outtree->Branch("bjetphi1", &bjetphi1, "bjetphi1/F");
 outtree->Branch("bjetphi2", &bjetphi2, "bjetphi2/F");
 outtree->Branch("bjete1", &bjetphi1, "bjetphi1/F");
 outtree->Branch("bjete2", &bjetphi2, "bjete2/F");
 // H2
 outtree->Branch("hww_mt", &hww_mt, "hww_mt/F");
 outtree->Branch("hww_pt", &hww_pt, "hww_pt/F");
 outtree->Branch("hww_etap", &hww_etap, "hww_etap/F");
 outtree->Branch("hww_etam", &hww_etam, "hww_etam/F");
 outtree->Branch("hww_phi", &hww_phi, "hww_phi/F");
 // H2 constituents
 outtree->Branch("gen_hww_mt", &gen_hww_mt, "gen_hww_mt/F");
 outtree->Branch("gen_hww_pt", &gen_hww_pt, "gen_hww_pt/F");
 outtree->Branch("gen_hww_phi", &gen_hww_phi, "gen_hww_phi/F");
 outtree->Branch("gen_hww_eta", &gen_hww_eta, "gen_hww_eta/F");
 //
 outtree->Branch("hw1_pt", &hw1_pt, "hw1_pt/F");
 outtree->Branch("hw1_eta", &hw1_eta, "hw1_eta/F");
 outtree->Branch("hw1_phi", &hw1_phi, "hw1_phi/F");
 outtree->Branch("hw1_e", &hw1_e, "hw1_e/F");
 outtree->Branch("hw2_pt", &hw2_pt, "hw2_pt/F");
 outtree->Branch("hw2_eta", &hw2_eta, "hw2_eta/F");
 outtree->Branch("hw2_phi", &hw2_phi, "hw2_phi/F");
 outtree->Branch("hw2_e", &hw2_e, "hw2_e/F");
 // for W
 outtree->Branch("pfmet", &pfmet, "pfmet/F");
 outtree->Branch("pt1", &pt1, "pt1/F");
 outtree->Branch("pt2", &pt2, "pt2/F");
 outtree->Branch("ptll", &ptll, "ptll/F");
 outtree->Branch("pzll", &pzll, "pzll/F");
 outtree->Branch("mll", &mll, "mll/F");
 outtree->Branch("dphill", &dphill, "dphill/F");
 outtree->Branch("gen_pfmet", &gen_pfmet, "gen_pfmet/F");
 outtree->Branch("gen_pfmez", &gen_pfmez, "gen_pfmez/F");
 outtree->Branch("gen_mvv", &gen_mvv, "gen_mvv/F");
 outtree->Branch("gen_pt1", &gen_pt1, "gen_pt1/F");
 outtree->Branch("gen_pt2", &gen_pt2, "gen_pt2/F");
 outtree->Branch("gen_ptll", &gen_ptll, "gen_ptll/F");
 outtree->Branch("gen_pzll", &gen_pzll, "gen_pzll/F");
 outtree->Branch("gen_mll", &gen_mll, "gen_mll/F");
 outtree->Branch("gen_dphill", &gen_dphill, "gen_dphill/F");
 outtree->Branch("nlep", &nlep, "nlep/F");
 outtree->Branch("channel", &channel, "channel/F");
 // ??
 outtree->Branch("xhh_ww_mt",  &xhh_ww_mt,  "xhh_ww_mt/F");
 outtree->Branch("xhh_p_ww_pt",  &xhh_p_ww_pt,  "xhh_p_ww_pt/F");
 outtree->Branch("xhh_p_ww_eta", &xhh_p_ww_eta, "xhh_p_ww_eta/F");
 outtree->Branch("xhh_p_ww_phi", &xhh_p_ww_phi, "xhh_p_ww_phi/F");
 outtree->Branch("xhh_p_ww_m",   &xhh_p_ww_m,   "xhh_p_ww_m/F");
 //
 outtree->Branch("xhh_m_ww_pt",  &xhh_m_ww_pt,  "xhh_m_ww_pt/F");
 outtree->Branch("xhh_m_ww_eta", &xhh_m_ww_eta, "xhh_m_ww_eta/F");
 outtree->Branch("xhh_m_ww_phi", &xhh_m_ww_phi, "xhh_m_ww_phi/F");
 outtree->Branch("xhh_m_ww_m",   &xhh_m_ww_m,   "xhh_m_ww_m/F");

  return 0;
}
