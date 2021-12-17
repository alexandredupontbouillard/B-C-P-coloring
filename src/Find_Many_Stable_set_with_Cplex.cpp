#include "Find_Many_Stable_set_with_Cplex.h"



using namespace std;

#define eps 1e-5
//#define CPLEX_LOG  //if non commented, prints logs

//#define WRITEPL






void Cplex_Multiple_pricing_algo::create_MIP(C_Graph *GG){

  
  #ifdef CPLEX_LOG
  cout<<"**************************** INIT CPLEX *********************"<<endl;
  nbIt =0;
  #endif
  
  G=GG;
  nb_additional_variable = 0;
  model= IloModel(env);


  ////////////////////////
  //////  VAR
  ////////////////////////

  int i,k;
  
  x = IloNumVarArray(env, G->nb_nodes, 0.0, 1.0, ILOINT);
  additional = IloNumVarArray(env, 0, 0.0, 1.0, ILOINT);
  for(i = 0; i < G->nb_nodes; i++) {
    ostringstream varname;
    IloNumVar v(env, 0.0, 1.0, ILOINT);
    varname.str("");
    varname<<"x_"<<i;
    x[i].setName(varname.str().c_str());

    
  }

  //////////////
  //////  EDGE INEQUALITIES
  //////////////

  IloRangeArray CC(env);
  int nbcst=0;

  // Cst x_i + x_j \le 1 for every edges ij in E
 
  for (k=0;k<G->nb_links;k++){
      IloExpr cst(env);
      cst+=x[G->V_links[k]->v1]+x[G->V_links[k]->v2];
      CC.add(cst<=1);
      ostringstream cstname;
      cstname.str("");
      cstname<<"CEdge_"<<G->V_links[k]->v1<<"_"<<G->V_links[k]->v2;
      CC[nbcst].setName(cstname.str().c_str());
      nbcst++;
  }
  

  model.add(CC);


  //////////////
  ////// OBJ
  //////////////


  // Initialization without any value
  obj = IloAdd(model, IloMaximize(env, 0.0));
  

   //////////////
  ////// CPLEX 
  //////////////

    


  cplex = IloCplex(model);

  #ifndef CPLEX_LOG
     cplex.setOut(env.getNullStream());
     cplex.setWarning(env.getNullStream());
   #endif

  #ifdef CPLEX_LOG
  cout<<"**************************** FIN INIT CPLEX *********************"<<endl;
  #endif



}



void Cplex_Multiple_pricing_algo::set_objective_coefficient(const vector<SCIP_Real>& obj_coeff, const vector<SCIP_Real>& coeff_additional){
 int i;
 
  for (i=0;i<G->nb_nodes;i++)
    obj.setLinearCoef(x[i],obj_coeff[i]);



	for(i = 0 ; i < coeff_additional.size(); i ++){
		obj.setLinearCoef(additional[i],coeff_additional[i]);
	
	}
}	

void Cplex_Multiple_pricing_algo::modify_pricer(list<Constraint*> ct){
		list<Constraint*>::iterator it;
		IloNumVar* ll;
		for(it = ct.begin() ; it !=ct.end(); it ++){
				
				ll = 	(*it)->modifyPricer(this);
				nb_additional_variable ++;
				

				additional.add(*ll);
				
		
		}


}


	
// return true if a stable set have been found by Cplex
bool  Cplex_Multiple_pricing_algo::find_stableset(StableSet* stable, SCIP_Real& objvalue, vector<int>* coefficients) {
  list<int>::const_iterator it;
   int i;

  #ifdef CPLEX_LOG
  cout<<endl<<" ************************* LAUCH PRICER with CPLEX"<<endl<<endl;
  #endif
 
	  
  #ifdef WRITEPL
  ostringstream filename;
  filename.str("");
  filename<<"../debug/pricing"<<nbIt<<".lp"<<"";
  nbIt++;
   

  cplex.exportModel(filename.str().c_str());
  #endif
     

  if ( !cplex.solve() ) {
    env.error() << "Failed to optimize Pricer with Cplex" << endl;
    exit(1);
  }
   

  if (cplex.getStatus()==CPX_STAT_INFEASIBLE){
    #ifdef CPLEX_LOG
    cout<<"NO SOLUTION BECAUSE INFEASIBILITY ---- OUCH Cela semble ne devoir jamais se produire----"<<endl;
    cout<<endl<<" ************************* END PRICER with CPLEX"<<endl<<endl;
    #endif
    
    return false;
  }
  else{
    
    

    objvalue=cplex.getObjValue();
    
    for(i=0;i<G->nb_nodes;i++){
      if (cplex.getValue(x[i])>1-eps)
	stable->add(i);
    }

    G->completeIntoMaxStable(stable);

    

   for(i=0; i < nb_additional_variable; i ++){	
   		
   		
   	(*coefficients)[i] = round(cplex.getValue(additional[i])); // 


   	
   }


    #ifdef CPLEX_LOG
    cout<<"Produced column : ";
    for (it=stable->L_nodes.begin();it!=stable->L_nodes.end();it++)
      cout<<*it<<" ";
    cout<<endl;
    cout<<"Objective value : "<<objvalue<<endl;
    cout<<endl<<" ************************* END PRICER with CPLEX"<<endl<<endl;     
    #endif
      
    return true;

  }
  

}





