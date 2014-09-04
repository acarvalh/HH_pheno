/// main11.cc is a part of the PYTHIA event generator.
// Copyright (C) 2011 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

#include "Pythia.h"
//#include <fstream>
using namespace Pythia8;
int main() {

  Pythia pythia;

string path= "/afs/cern.ch/work/a/acarvalh/phenoHH/model_LHEfiles/nonresonantHH/Events/13tev/parton/";
for(unsigned ifile=4; ifile<6; ifile++){
//unsigned ifile=0;
string namefile_in=path + "pp_hh_vbf_";
if(ifile==0) namefile_in += "BSM_13tev_VBFcuts_CV_p0p5_C2V_p1p0_C3_p1p0.lhe";
if(ifile==1) namefile_in += "BSM_13tev_VBFcuts_CV_p1p5_C2V_p1p0_C3_p1p0.lhe";
//
if(ifile==2) namefile_in += "BSM_13tev_VBFcuts_CV_p1p0_C2V_p2p0_C3_p1p0.lhe";
if(ifile==3) namefile_in += "BSM_13tev_VBFcuts_CV_p1p0_C2V_p0p0_C3_p1p0.lhe";
//
if(ifile==4) namefile_in += "BSM_13tev_VBFcuts_CV_p1p0_C2V_p1p0_C3_p0p0.lhe";
if(ifile==5) namefile_in += "BSM_13tev_VBFcuts_CV_p1p0_C2V_p1p0_C3_p2p0.lhe";
//
if(ifile==6) namefile_in += "BSM_13tev_VBFcuts_CV_p1p0_C2V_p1p0_C3_p1p0.lhe";
if(ifile==7) namefile_in += "SM_13tev_nocuts.lhe";
if(ifile==8) namefile_in += "SM_13tev_VBFcuts.lhe";
/*string namefile_in=path + "Madgraphcg";
if(ifile==0) namefile_in += "0/MGraviton_260.lhe";
if(ifile==1) namefile_in += "0/MGraviton_500.lhe";
if(ifile==2) namefile_in += "1/MGraviton_260.lhe";
if(ifile==3) namefile_in += "1/MGraviton_500.lhe";
if(ifile==4) namefile_in += "0_0137/Graviton_Parton/MGraviton_260.lhe";
if(ifile==5) namefile_in += "0_0137/Graviton_Parton/MGraviton_500.lhe";
if(ifile==6) namefile_in += "0/MGraviton_260.lhe";
if(ifile==7) namefile_in += "600.lhe";
if(ifile==8) namefile_in += "650.lhe";
if(ifile==9) namefile_in += "700.lhe";
if(ifile==10) namefile_in += "750.lhe";
*/
string namefile_out=namefile_in + ".decayed";

    //    string namefile_in=path + "atEightTeV_events_patched.lhe";
    //string namefile_out=namefile_in + ".pythia";
    //namefile_in += "test-MR610.lhe";
    cout<<"\n namefile_in = "<<namefile_in<<endl;
    cout<<"\n namefile_out = "<<namefile_out<<endl;

   
    // output file
    // we want to store the list of all final state particles
    ofstream out_pythia;
    // Highest precision required for jet clustering
    out_pythia.precision(15);
    // Generator. We here stick with default values, but changes
    // could be inserted with readString or readFile
    // Initialize Les Houches Event File run. List initialization information.
    pythia.readString("Beams:frameType = 4");
     // the analysis program
    string sfile = "Beams:LHEF ="+namefile_in;
    pythia.readString(sfile.c_str());
    out_pythia.open(namefile_out.c_str());
  // Allow for possibility of a few faulty events.
  int nAbort = 10;
  int iAbort = 0;
 // turn of hadronization settings - for testing  
 //  pythia.readString("25:mayDecay = no");
 ////////////////////////////////////////////////////////////////////
 // read decay table
 //pythia.readString("ProcessLevel:resonanceDecays = off"); // do not decay anything
pythia.readString("SLHA:readFrom = 2");
pythia.readString("SLHA:file = Susy.txt "); // input the decay table
////////////////////////////////////////////////////////////////////
pythia.readString("PartonLevel:MI = off"); // Off multiple interactions
pythia.readString("PartonLevel:ISR = off"); // Shower on
pythia.readString("PartonLevel:FSR = off"); // Shower on
pythia.readString("PartonLevel:FSRinResonances  = off"); // Off multiple interactions
pythia.readString("HadronLevel:all = off"); // Of hadronization

  pythia.init();

   for (int iEvent = 0; ; ++iEvent) {
    cout<<"\n ievent = "<<iEvent<<"\n"<<endl;
    // Generate events, and check whether generation failed.
    if (!pythia.next()) {
      // If failure because reached end of file then exit event loop.
      if (pythia.info.atEndOfFile()) break;
      // First few failures write off as "acceptable" errors, then quit.
      if (++iAbort < nAbort) continue;
      break;
    }
    cout<<"hi"<<endl;
    // Acess event record  pythia.event.size()
    cout<<"Number of particles = "<<pythia.event.size()<<endl;
    vector<int> pID;
    vector<double> px;
    vector<double> py;
    vector<double> pz;
    vector<double> E;
    vector<int> mother;
    vector<int> code;
    // Some checks on the event record
    // Check for example that at least we have two bs and two bbars
    for (int i = 0; i < pythia.event.size(); i++){
      int particle_id = pythia.event[i].id();
      int particle_status = pythia.event[i].status();
      int particle_mother = pythia.event[i].mother1();
      // save only final state particles
      if(particle_status>0){
        cout<<i<<" "<<particle_id<<" "<<particle_mother<<" "<<particle_status<<endl;
        double ppx= pythia.event[i].px();
        double ppy= pythia.event[i].py();
        double ppz= pythia.event[i].pz();
        double EE= pythia.event[i].e();
        //cout<<px<<" "<<py<<" "<<pz<<" "<<E<<endl;
        pID.push_back(particle_id);
        px.push_back(ppx);
        py.push_back(ppy);
        pz.push_back(ppz);
        E.push_back(EE);
        mother.push_back(particle_mother);
        code.push_back(particle_id);
      }
    }
    // Save into file
    out_pythia<<"#"<<endl;
    cout<<"Number of final state particles = "<<E.size()<<"\n"<<endl;
    out_pythia<<E.size()<<endl;
     for(unsigned i=0;i<E.size();i++){
       out_pythia<<pID.at(i)<<" "<<px.at(i)<<" "<<py.at(i)<<" "<<pz.at(i)<<" "<<E.at(i)<<" "<<endl;
     }
    

  // End of event loop.
  }

  out_pythia.close();

  // Give statistics. Print histogram.
  pythia.statistics();


  } //for unsigned

  // Done.
  return 0;
}
