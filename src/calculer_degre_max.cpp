

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "Graph.h"
#include <cstring>



/* namespace usage */
using namespace std;
	







int main(int argc, char** argv){



	char * name,*nameext; //, *nameextsol;
  int i;
	
	  if(argc!=2){
	    cerr<<"usage: "<<argv[0]<<" <DIMACS file name> "<<endl; 
	    return 1;
	  }


  srand(1);
  
  name= new char[40];
  nameext= new char[40];
  //nameextsol= new char[40];



  name=strcat(name,argv[1]);
  nameext=strcat(nameext,argv[1]);


  ifstream fic(nameext);
  
  C_Graph G;

  G.read_undirected_DIMACS(fic);
	int max = 0;
  for(int i = 0 ; i < G.nb_nodes ; i++){
	if(max < G.V_nodes[i].L_adjLinks.size()){
		max = G.V_nodes[i].L_adjLinks.size();

	}


   }
   cout<<max<<endl;



















	return 0;
}
