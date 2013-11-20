float HiggsMass = 125.;
// lepton - jet isolation
double DRmax = 0.2;
// Jets
float MINPTJET = 15.;  


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// Global parameters of Jet Finding
// Basic kinematical cuts
// To be applied only to hadron level events
double const rapmax=2.5;

bool const doplots = true;

// Parameters of the jet clustering
// Jet radius for the basic jet clustering
double const jetR = 0.7;

// Set parameters of mass drop
// mu = 0.67 and y = 0.09 are the default choice in FastJet
double const mu = 0.67;
double const ycut = 0.09;
// This parameters limits the efficiency as (1-ycut)^2
//double const ycut = 0.0;

// Separation in rapidity of the two Higgs candidates
double const Delta_y_max = 6;//  1.3; //optimized from dijet searches
double const fm=0.15;

// Parameters of jet filtering
double const Rfilt = 0.3;
int const n_subjet =3;

// Parameters of the jet substructure
double const Rsb = 1.1;

// Parameters of the b tagging
// Check pt >10 for at least one the bs

// photon isolation variables
