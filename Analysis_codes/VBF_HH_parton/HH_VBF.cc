//////////////////
// to run:
// make HH_VBF.cc
/////////////////
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include <fstream>
using namespace fastjet;
using namespace std;
#include "Functions.h"
#include "choices.h"

int main() {
  cout<<"\n HI!!!!!!!!!!! \n "<< nevent <<endl;
  // Initialize random seed
  srand( time(NULL) );

}
