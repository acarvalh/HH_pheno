/// main11.cc is a part of the PYTHIA event generator.
// Copyright (C) 2011 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

#include "Pythia.h"
//#include <fstream>
using namespace Pythia8;
int main() {


  // Interface for conversion from Pythia8::Event to HepMC one. 
//  HepMC::I_Pythia8 ToHepMC;
  //  ToHepMC.set_crash_on_problem();

  // Specify file where HepMC events will be stored.
//  HepMC::IO_GenEvent ascii_io(argv, std::ios::out);
  // Following line is a deprecated alternative, removed in recent versions
  // HepMC::IO_Ascii ascii_io("hepmcout32.dat", std::ios::out);
  // Line below is an eye-readable one-way output, uncomment the include above
  // HepMC::IO_AsciiParticles ascii_io("hepmcout32.dat", std::ios::out);
  // Generator. Shorthand for the event.
  Pythia pythia;

  //unsigned const nfiles=4;
  string path= "/afs/cern.ch/work/a/acarvalh/phenoHH/pythia8153/examples/bulk_graviton/";


for(unsigned ifile=0; ifile<1; ifile++){
//unsigned ifile=0;
string namefile_in=path + "RSG_WBF_hh-Mhh";
if(ifile==0) namefile_in += "260.lhe";
if(ifile==1) namefile_in += "300.lhe";
if(ifile==2) namefile_in += "350.lhe";
if(ifile==3) namefile_in += "400.lhe";
if(ifile==4) namefile_in += "450.lhe";
if(ifile==5) namefile_in += "500.lhe";
if(ifile==6) namefile_in += "550.lhe";
if(ifile==7) namefile_in += "600.lhe";
if(ifile==8) namefile_in += "650.lhe";
if(ifile==9) namefile_in += "700.lhe";
if(ifile==10) namefile_in += "750.lhe";
string namefile_out=namefile_in + ".pythia";

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
    // could be inserted with readString or readFile.

    
    // Initialize Les Houches Event File run. List initialization information.
    pythia.readString("Beams:frameType = 4");
     // the analysis program
    string sfile = "Beams:LHEF ="+namefile_in;
    pythia.readString(sfile.c_str());
    out_pythia.open(namefile_out.c_str());
 
 // Interface for conversion from Pythia8::Event to HepMC event.
 // HepMC::I_Pythia8 ToHepMC;
 
 // Specify file where HepMC events will be stored.
 //HepMC::IO_GenEvent ascii_io(namefile_out.c_str(), std::ios::out);

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


//pythia.readString("HadronLevel:all = on"); // Of hadronization

// read decay table
pythia.readString("SLHA:readFrom = 2");
pythia.readString("SLHA:file = Susy.txt "); // input the decay table
// allow overwrite: only works for products of SM - like decays
pythia.readString("SLHA:allowUserOverride = on ");
//pythia.readString("25:mayDecay = no"); 
pythia.readString("24:onMode = off");
pythia.readString("-24:onMode = off");
pythia.readString("24:onIfMatch = 12 11"); // e ve
pythia.readString("24:onIfMatch = 14 13"); // mu numu
pythia.readString("-24:onIfMatch = 12 -11"); // e ve
pythia.readString("-24:onIfMatch = 14 -13"); // mu numu

  // Create an LHAup object that can access relevant information in pythia.
  LHAupFromPYTHIA8 myLHA(&pythia.event, &pythia.info);

  // Open a file on which LHEF events should be stored, and write header.
  myLHA.openLHEF("test4b.lhe"); 

  pythia.init();

  // Store initialization info in the LHAup object. 
  myLHA.setInit();

  // Write out this initialization info on the file.
  myLHA.initLHEF();

  // Begin event loop; generate until none left in input file.
  for (int iEvent = 0; iEvent<1; ++iEvent) {

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
    // text output -- for test
/*
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
  */  
    // Construct new empty HepMC event. Form with arguments is only
    // meaningful for HepMC 2.04 onwards, and even then unnecessary  
    // if HepMC was built with GeV and mm as units from the onset. 
    //HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
    //HepMC::GenEvent* hepmcevt = new HepMC::GenEvent(HepMC::Units::GEV, HepMC::Units::MM); 

    // Fill HepMC event, including PDF info.
    //ToHepMC.fill_next_event( pythia, hepmcevt );
    // This alternative older method fills event, without PDF info.
    // ToHepMC.fill_next_event( pythia.event, hepmcevt );

    // Write the HepMC event to file. Done with it.
    //ascii_io << hepmcevt;
    //delete hepmcevt;
    ///////////////////////////////////////////// 
    // Store event info in the LHAup object. 
    myLHA.setEvent();
    // Write out this event info on the file. 
    myLHA.eventLHEF();

  // End of event loop.
  }

  out_pythia.close();

  // Give statistics. Print histogram.
  pythia.statistics();

  // Update the cross section info based on Monte Carlo integration during run.
  myLHA.updateSigma();

  // Write endtag. Overwrite initialization info with new cross sections.
  myLHA.closeLHEF(true);
  //delete pythia;

  } //for unsigned

  // Done.
  return 0;
}
