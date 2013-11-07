#include<math.h>
#include<stdio.h>
#include<unistd.h>
#include<sys/stat.h>
#include<sys/types.h>
#include <dlfcn.h>
#include <sys/wait.h> 
#include"num_in.h"
#include"num_out.h"
#include"VandP.h"
#include"dynamic_cs.h"
#include"rootDir.h" 
#include <time.h>

/* To run

./bin/make_main scan_example.c
./a.out

*/
int main(void)
{ int err,i;

	/* INTPUT PARAMETERS (to scan over) */
    double Mhh,	MRmin=100, 	MRmax=3000,   MRstep=10;
 
	/* OUTPUT PARAMETERS */
    // Higgs decay branching ratios
    double  wR,braa,brgg,brww,brzz,brhh,brtt,brbb,brqq,brtata,brll,brnee,brnmm,brntata,brnn;
    double  bru,brd,brc,brs,bree,brmm;
    txtList branchings;

//set model dir here
char mdldir[] = "models";

 // Set model number and number of points to collect, mdlnr is your model number
int mdlnr=10, npoints=500;
// 3 bulk RS
// 9 brane RS
// 10 RS1
//
//a model to switch between to reset values when reloading
 setModel(mdldir , mdlnr ); 

/*****************************************************************************/
 srand (time(NULL)); //this is used to seed the random number by the system time

 if (remove("scan.dat") == -1)
	perror("Error in deleting a file");

 FILE *file;
 file = fopen("scan.dat","a+"); /* apend file (add text to
					a file or create a file if it does not exist.*/

 // Writing parameter names at first line to keep track of columns:
 //input parameters (1)
 //output parameters (3)
 fprintf(file,"MGr\t\twR\t\tbraa\t\tbrgg\t\tbrww\t\tbrzz\t\tbrhh\t\tbrtt\t\tbrbb\t\tbrqq\t\tbrtata\t\tbrll\t\tbrnn\n");		
 fclose(file); /*done with header of file*/

 /*** Starting randomizing loop ***/
 for (i = 0; i <= npoints; i++){

 /********** generate random values for variables **********/
 Mhh     = MRmin+ i*(MRmax-MRmin)/npoints;

 /* Have to reset model every time, otherwise widths are not recalculated */
 setModel(mdldir , mdlnr ); 

 /******* assign variable values ********/
 /* the string is the calchep var name */
 //tan beta
 err=assignValW("Mhh", Mhh);

 // Calculation of public constraints  
 err=calcMainFunc();

 if(err!=0) { 
	  printf("Can not calculate constrained parameter %s\n",varNames[err]);i--;
 }
 else {
		// if the point survives the constraints collect more output values:
		// width and branchings of a particle
		wR     = pWidth("hh",&branchings);
		braa	= findBr(branchings,"A,A");
		brgg	= findBr(branchings,"G,G");	
		brww	= findBr(branchings,"W+,W-");	
		brzz	= findBr(branchings,"Z,Z");	
		brhh	= findBr(branchings,"h,h");
		brtt	= findBr(branchings,"t,T");	
		brbb	= findBr(branchings,"b,B");
                //
		bru	= findBr(branchings,"u,U");
		brd	= findBr(branchings,"d,D");
		brc	= findBr(branchings,"c,C");
		brs	= findBr(branchings,"s,S");
		brqq=bru+brd+brc+brs;
		//
		brtata	= findBr(branchings,"l,L");			
		bree	= findBr(branchings,"e,E");
		brmm	= findBr(branchings,"m,M");
		brll = brmm+bree;
		//
		brntata	= findBr(branchings,"nl,Nl");			
		brnee	= findBr(branchings,"ne,Ne");
		brnmm	= findBr(branchings,"nm,Nm");
		brnn = brnmm+brnee+brntata;
		// write values to file
  		file  = fopen("scan.dat","a+");
		//input parameters
  		fprintf(file,"%f\t",Mhh);
		//output parameters
  		fprintf(file,"%f\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",wR,braa,brgg,brww,brzz,brhh,brtt,brbb,brqq,brtata,brll,brnn);
  		fclose(file); 
 }
  
 }// *** end of rand loop ***

  return 0;
}
