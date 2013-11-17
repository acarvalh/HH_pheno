// main11.cc is a part of the PYTHIA event generator.
// Copyright (C) 2011 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

#include <fstream>
#include "Pythia.h"

using namespace Pythia8;
int main() {
  // Generator. Shorthand for the event.
  Pythia pythia;

  //unsigned const nfiles=4;
  string path= "/afs/cern.ch/work/a/acarvalh/phenoHH/pythia8153/examples/";

/*
for(unsigned ifile=0; ifile<11; ifile++){
//unsigned ifile=0;
string namefile_in=path + "MR_";
if(ifile==0) namefile_in += "test-MR610.lhe";
if(ifile==1) namefile_in += "500_0.lhe";
if(ifile==2) namefile_in += "600_0.lhe";
if(ifile==3) namefile_in += "700_0.lhe";
if(ifile==4) namefile_in += "800_0.lhe";
if(ifile==5) namefile_in += "900_0.lhe";
if(ifile==6) namefile_in += "1000_0.lhe";
if(ifile==7) namefile_in += "1500_0.lhe";
if(ifile==8) namefile_in += "2000_0.lhe";
if(ifile==9) namefile_in += "2500_0.lhe";
if(ifile==10) namefile_in += "3000_0.lhe";
string namefile_out=namefile_in + ".pythia";
*/
    string namefile_in=path + "atEightTeV_events_patched.lhe";
    string namefile_out=namefile_in + ".pythia";
    //namefile_in += "test-MR610.lhe";
    cout<<"\n namefile_in = "<<namefile_in<<endl;
    cout<<"\n namefile_out = "<<namefile_out<<endl;

   
    // output file
    // we want to store the list of all final state particles
    ofstream out_pythia;
    // Highest precision required for jet clustering
    out_pythia.precision(15);

    // Generator. We here stick with default values, but changes
    // could be inserted with readString or readFile.

    
    // Initialize Les Houches Event File run. List initialization information.
    pythia.readString("Beams:frameType = 4");
     // the analysis program
    string sfile = "Beams:LHEF ="+namefile_in;
    pythia.readString(sfile.c_str());
    out_pythia.open(namefile_out.c_str());
 
  // Allow for possibility of a few faulty events.
  int nAbort = 10;
  int iAbort = 0;

  // pythia.settings.listAll();
  
  // Settings


  // turn of hadronization settings - for testing  

pythia.readString("PartonLevel:MI = off"); // Off multiple interactions
pythia.readString("PartonLevel:ISR = off"); // Shower on
pythia.readString("PartonLevel:FSR = off"); // Shower on
pythia.readString("HadronLevel:all = off"); // Of hadronization

//pythia.readString("SUSY:idB  = 0 ");
//pythia.readString("SUSY:idVecB  = 0 ");
/*
//pythia.readString("35:mayDecay = no");
pythia.readString("25:mayDecay = no");

pythia.readString("35:onMode = off");
// pythia.readString("25:onIfMatch = 5 -5"); // b
// pythia.readString("25:onIfMatch = 4 -4"); // c
// pythia.readString("25:onIfMatch = 21 21"); // gluon
pythia.readString("35:onIfMatch = 24 24"); // W
pythia.readString("24:mayDecay = no");


pythia.readString("21:mayDecay = no");
pythia.readString("4:mayDecay = no");
pythia.readString("-4:mayDecay = no");
pythia.readString("5:mayDecay = no");
pythia.readString("-5:mayDecay = no");
pythia.readString("5:mayDecay = no");
pythia.readString("-5:mayDecay = no");
pythia.readString("11:mayDecay = no");
pythia.readString("-11:mayDecay = no");
pythia.readString("13:mayDecay = no");
pythia.readString("-13:mayDecay = no");

*/
pythia.readString("SLHA:readFrom = 2");
pythia.readString("SLHA:file = Susy.txt ");

pythia.readString("SLHA:allowUserOverride = on ");
pythia.readString("24:onMode = off");
pythia.readString("-24:onMode = off");
pythia.readString("24:onIfMatch = 12 11"); // e ve
pythia.readString("24:onIfMatch = 14 13"); // mu numu
pythia.readString("-24:onIfMatch = 12 -11"); // e ve
pythia.readString("-24:onIfMatch = 14 -13"); // mu numu


  // turn decays
//  pythia.readString("HadronLevel:all = on"); // generic hadronization settings
  // hadronic Higgs decay
/*  pythia.readString("25:onMode = off");
  8pythia.readString("25:onIfMatch = 5 5"); // b
  pythia.readString("25:onIfMatch = 4 -4"); // c
  pythia.readString("25:onIfMatch = 21 21"); // gluon
  // WW "higgs" decay
  pythia.readString("35:onMode = off");
  pythia.readString("35:onIfMatch = 24 24"); // W
  // W decay
//  pythia.readString("24:onMode = off");
//  pythia.readString("24:onIfMatch = 12 11"); // e ve
//  pythia.readString("24:onIfMatch = 14 13"); // mu numu
*/
  // Initialization
  pythia.init();

  // Begin event loop; generate until none left in input file.
  for (int iEvent = 0; iEvent<10 ; ++iEvent) {

    cout<<"\n ievent = "<<iEvent<<"\n"<<endl;

    // Generate events, and check whether generation failed.
    if (!pythia.next()) {

      // If failure because reached end of file then exit event loop.
      if (pythia.info.atEndOfFile()) break;

      // First few failures write off as "acceptable" errors, then quit.
      if (++iAbort < nAbort) continue;
      break;
    }

    // Acess event record
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
        cout<<i<<" "<<particle_id<<" "<<particle_mother<<endl;
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

  //delete pythia;

 // } //for unsigned

  // Done.
  return 0;
}
