#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

/* scip includes */
#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"

/* user defined includes */
#include "Graph.h"
#include "SCIP_master.h"
#include "SCIP_pricer.h"
//#include "OddHoleSeparator.h"
#include "OddHoleSep.h"
#include "ChvatalGomorySep.h"
#include "ChvatalGomoryMultipleSep.h"
#include "SCIP_pricer_alternate.h"
#include <chrono>
#include "technique.h"


//#include "constraints.h"


/* namespace usage */
using namespace std;
using namespace scip;


#define NBMAXCONSITER 2000


//#define SCIP_DEBUG   // Si non commente, affiche le log de SCIP


int main(int argc, char** argv)
{
  char * name,*nameext; //, *nameextsol;
  int i;
	

  if(argc!=3){
    cerr<<"usage: "<<argv[0]<<" <DIMACS file name> <mode>"<<endl; 
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

/*
  StableSet s1;
  StableSet s2;
  StableSet s3;
	s1.initialize(7);
	s2.initialize(7);
  	s3.initialize(7);
	s1.add(0);
	s1.add(5);
	s1.add(3);
	s2.add(1);
	s2.add(6);
	s3.add(2);
	s3.add(4);
	s3.add(6);
  vector<StableSet*> ss;
  ss.push_back(&s1);
  ss.push_back(&s2);
  ss.push_back(&s3);

  list<list<int>*> rr;
   vector<float> cost;
	cost.push_back(0.5);
	cost.push_back(0.7);
	cost.push_back(0.5);
  G.findViolatedCycleIq( &ss, cost, rr, 3);
  cout<<"number of cycles : "<<rr.size()<<endl;
	  list<list<int>*>::iterator gege;
  list<int>::iterator hehe;
	for(gege = rr.begin(); gege != rr.end(); gege++){
		for(hehe = (*gege)->begin(); hehe != (*gege)->end(); hehe++){

			cout<<*hehe<<"";

		}
		cout<<endl;


	}

cout<<"\n\n\n\n\n\n\n\n\n\n"<<endl;*/
  fic.close();
  

  //////////////
  //////  SCIP INITIALIZATION
  //////////////

  SCIP *scip=NULL;
  
  C_master_coloring Master;
  Master._data.filename = nameext;

  Master._data.nbChvatal=0;
  Master._data.separationTime = 0;
  Master._data.pricingTime = 0;
  Master.separatingPhase = 1;
	
  
  // initialize SCIP environment 
  SCIPcreate(&scip);

   SCIPprintVersion(scip, NULL);
   //SCIPinfoMessage(scip, NULL, "\n");

   // include default plugins 
   SCIPincludeDefaultPlugins(scip);

   // set verbosity parameter 
   #ifdef SCIP_DEBUG
   SCIPsetIntParam(scip, "display/verblevel", 5);
   // SCIPsetBoolParam(scip, "display/lpinfo", TRUE); 
   #endif


  // create empty problem 
   SCIPcreateProb(scip, "COLORING", 0, 0, 0, 0, 0, 0, 0);

   SCIP_CALL( SCIPsetLongintParam(scip, "limits/totalnodes", 1) );
   SCIP_CALL( SCIPsetRealParam(scip, "limits/time", 3600) );

  // SCIP_CALL( SCIPsetRealParam(scip, "separating/minortho", 0.0001) );
   //SCIP_CALL( SCIPsetRealParam(scip, "separating/minorthoroot", 0.0001) );
   //SCIP_CALL( SCIPsetCharParam(scip, "separating/efficacynorm", 's') );
    //SCIP_CALL( SCIPsetIntParam(scip, "separating/maxstallrounds",10));
   //SCIP_CALL( SCIPsetIntParam(scip, "separating/maxruns", 1) );
   
  ////////////////////////
  //////  INEQUALITIES
  ////////////////////////


  
   // Set covering constraints
   

   Master.V_node_ineq.resize(G.nb_nodes, (SCIP_CONS*) NULL);
   
   char con_name[255];
   for (i = 0; i < G.nb_nodes; i++)
   {
      SCIP_CONS* con = NULL;
      (void) SCIPsnprintf(con_name, 255, "C%d", i);
      SCIPcreateConsLinear( scip, &con, con_name, 0, NULL, NULL,
			    1.0,   // lhs 
			    SCIPinfinity(scip),   // rhs  SCIPinfinity(scip) if >=1
			    true,  // initial
			    false, // separate 
			    true,  // enforce 
			    true,  // check 
			    true,  // propagate 
			    false, // local 
			    true,  // modifiable 
			    false, // dynamic 
			    false, // removable 
			    false  /* stickingatnode */);
      SCIPaddCons(scip, con);
      Master.V_node_ineq[i] = con; 
   }

  /////////////////////
  // Add variables corresponding to the glout coloring
  ////////////////////


   Master.L_var.clear();

   list<StableSet*> firstColor;
   G.getGloutColoring(firstColor);
   list<StableSet*>::iterator it;
   list<int>::iterator itt;
   char var_name[255];
   i = 0;
   for (it = firstColor.begin() ;it != firstColor.end();it++){
     SCIPsnprintf(var_name, 255, "V_%d",i);
     i++;
     SCIPdebugMsg(scip, "new variable <%s>\n", var_name);

     // create the new variable: Use upper bound of infinity such that we do not have to care about
     // the reduced costs of the variable in the pricing. The upper bound of 1 is implicitly satisfied
     // due to the set partitioning constraints.
      
     C_master_var* var=new C_master_var;
     
     SCIPcreateVar(scip, &(var->ptr), var_name,
		   0,                     // lower bound
		   SCIPinfinity(scip),      // upper bound
		   1,                     // objective
		   SCIP_VARTYPE_INTEGER, // variable type
		   true, false, NULL, NULL, NULL, NULL, NULL);

          
     // add new variable to the list of variables to price into LP (score: leave 1 here) 
     SCIPaddVar(scip, var->ptr);

    var->stable=(*it);

     // add coefficients into the set covering constraints 
    for(itt = (*it)->L_nodes.begin(); itt != (*it)->L_nodes.end(); itt++){
	
     	SCIPaddCoefLinear(scip, Master.V_node_ineq[*itt], var->ptr, 1.0);
    }


     
     
     
     Master.L_var.push_back(var);

   }
      Master._data.nbColonnes=firstColor.size();






    
   ////////////////////
   // Define pricer
   //////////////////
   
   static const char* PRICER_NAME = "Pricer_Coloring";

  
 bool* b;
 if( argv[2][0] ==  'c'  && argv[2][1] != 'a' ){
	    cout<<"yo"<<endl;
	    ChvatalGomoryMultipleSep* sep2 = new ChvatalGomoryMultipleSep(scip,&Master,G.nb_nodes); 
	    //ChvatalGomorySep* sep2 = new ChvatalGomorySep(scip,&Master,G.nb_nodes); 
	    SCIP_CALL(SCIPincludeObjSepa(scip,sep2,FALSE));

	    ObjPricerColoring* pricer_ptr;
	 pricer_ptr = new ObjPricerColoring(scip, PRICER_NAME, &G, &Master,true);
	 SCIPincludeObjPricer(scip, pricer_ptr, true);
	(*b) = true;

         


     } 
   else if(argv[2][1] == 'a' && argv[2][0] ==  'c' ){  // THIS IS THE SEPARATOR/PRICER COMBO THAT DOES NOT WORK
	   PricerAlternate * pricer_ptr;
	   pricer_ptr = new PricerAlternate(scip, PRICER_NAME, &G, &Master);
	   SCIPincludeObjPricer(scip, pricer_ptr, true);
	   SCIP_CALL(SCIPincludeObjSepa(scip,pricer_ptr->separator,FALSE));
	  b = &pricer_ptr->costIsNegative;
	   
   }
   else{
	 ObjPricerColoring* pricer_ptr;
	 pricer_ptr = new ObjPricerColoring(scip, PRICER_NAME, &G, &Master,true);
	 SCIPincludeObjPricer(scip, pricer_ptr, true);
	(*b)=true;
   }


   




		
    // activate pricer 
    SCIPactivatePricer(scip, SCIPfindPricer(scip, PRICER_NAME));




   
	////////////////SOLVE/////////////
	/////////////////////////////////
   auto start = chrono::steady_clock::now();
   
   SCIPsolve(scip);

   auto end = chrono::steady_clock::now();
   ofstream res;
   res.open("../result.txt",std::ios_base::app);
	
   
/*
   res << nameext << ";" <<  SCIPgetLPObjval(scip)   <<";"  <<Master.L_var.size() << ';' << G.nb_nodes << ';' << G.nb_links <<";"<<Master.additional_ineq.size() << ";"<< ((double)chrono::duration_cast<chrono::milliseconds>(end - start).count()) /1000 << ";"<< Master._data.separationTime/1000 <<";" << Master._data.pricingTime/1000 <<";" << isSolEntier(scip,&Master)<<";"<<(*b) << "\n";


   Master._data.writeData();
*/
   #ifdef SCIP_DEBUG
   cout<<"Write init pl"<<endl;
  SCIPwriteTransProblem(scip, "../init.lp", "lp", TRUE);master->_data.dualValue(0);
   #endif



   res.close();
   //writting columns
	/*
   ofstream col;
   col.open("../wholeLP/columns.txt",std::ios_base::out);
   list<C_master_var*>::iterator ite;
   list<int>::iterator stabit; 

   for(ite = Master.L_var.begin(); ite != Master.L_var.end(); ite++){  
		


		for(stabit = (*ite)->stable->L_nodes.begin();stabit != (*ite)->stable->L_nodes.end(); stabit++ ){
			col<< *stabit << " ";
		   
		}
		col<<endl;


	}

   



   col.close();
	*/

   //////////////STATISTICS//////
   //////////////////////////////
    //SCIPprintStatistics(scip, NULL);
   cout<< "Nb de contraintes générées "<< Master.additional_ineq.size()<< endl;

   cout<<"Nb de colonnes initialisation+pricing: "<<Master.L_var.size()<<endl;

   list<C_master_var*>::iterator itv;
   list<int>::iterator iti;   

   
   #ifdef SCIP_DEBUG
  
   SCIP_SOL *sol=SCIPgetBestSol(scip);
   Master.printSolution( scip,sol);
   #endif


   //  nameextsol=strcat(nameextsol,argv[1]); 
   // nameextsol=strcat(nameextsol,".color");
   //  ofstream ficsol(nameextsol);
   //for(i = 0; i < G.nb_nodes; i++) 
   // ficsol<<V_sol[i]<<" ";
   //ficsol.close();
   

   //DEINITIALISATION/////
   //////////////////////
  

   for (itv=Master.L_var.begin();itv!=Master.L_var.end();itv++){
     SCIPreleaseVar(scip, &((*itv)->ptr));
   }
   
   for (i = 0; i < G.nb_nodes; i++){     
     SCIPreleaseCons(scip, &Master.V_node_ineq[i]);
   }
	
   
   //SCIPfree(&scip);  // Ne devrait pas produire de bug... mais il y en a


   BMScheckEmptyMemory();

   return 0;
}


void C_master_coloring::printSolution(  SCIP * scip, SCIP_SOL *sol){
   list<C_master_var*>::iterator itv;
   list<int>::iterator iti;  
   //int nbcol = 0;

   cout<<"Affichage de la solution par nous"<<endl;
   for(itv = this->L_var.begin(); itv!=this->L_var.end(); itv++){
       (*itv)->printSolution( scip, sol);

   }

   cout<<"Affichage de la solution par SCIP"<<endl;
   SCIPprintBestSol(scip, NULL, FALSE);


}

void C_master_var::printSolution( SCIP * scip , SCIP_SOL *sol){


	if(SCIPgetSolVal (scip, sol, ptr) > SCIPepsilon(scip)){
		stable->print();

	}

	

}

