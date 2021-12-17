#include "SCIP_pricer.h"
#include "scip/cons_linear.h"
#include <map>
#include <vector>
#include <iostream>

//#define OUTPUT_PRICERR  // Si non commente, affiche les logs
//#define OUTPUT_PRICER
#define epsilon 0.000001


#define separation 1

using namespace std;
using namespace scip;



SCIP_Real rounding(SCIP* scip, SCIP_Real r){
	//return r;
	if(r <= 0.0000001)
		return SCIPepsilon(scip);

	if(SCIPisEQ(scip,r,SCIPepsilon(scip)))	
		return SCIPepsilon(scip);


	
	else
		return r;
		return floor(r*1000000)/1000000;

}


/** Constructs the pricer object with the data needed
 *
 *  An alternative is to have a problem data class which allows to access the data.
 */
ObjPricerColoring::ObjPricerColoring(
   SCIP*                                scip,          /**< SCIP pointer */
   const char*                         pp_name,      /**< name of pricer */  
   C_Graph*                             GG,
   C_master_coloring*                      MM,
   bool inhibit_alternatee
   ):
  ObjPricer(scip, pp_name, "Finds stable set with negative reduced cost.", 0, TRUE)
{
  
  G=GG;  // save a pointer on data graph
  M=MM;
  nbNewConstaints = 0;
  inhibit_alternate = inhibit_alternatee;
  A_cplex.create_MIP(G);


 
}


/** Destructs the pricer object. */
ObjPricerColoring::~ObjPricerColoring()
{
  cout<<"Destructeur du pricer"<<endl;
}


/** initialization method of variable pricer (called after problem was transformed)
 *
 *  Because SCIP transformes the original problem in preprocessing, we need to get the references to
 *  the variables and constraints in the transformed problem from the references in the original
 *  problem.
 */
SCIP_DECL_PRICERINIT(ObjPricerColoring::scip_init)
{
  int i;
 
  for (i = 0; i < G->nb_nodes; i++){
    SCIPgetTransformedCons(scip, M->V_node_ineq[i], &(M->V_node_ineq)[i]);
  }
   
   return SCIP_OKAY;
}

/** Pricing of additional variables if LP is feasible.
 *
 *  - get the values of the dual variables you need
 *  - construct the reduced-cost arc lengths from these values
 *  - find the shortest admissible tour with respect to these lengths
 *  - if this tour has negative reduced cost, add it to the LP
 *
 *  possible return values for *result:
 *  - SCIP_SUCCESS    : at least one improving variable was found, or it is ensured that no such variable exists
 *  - SCIP_DIDNOTRUN  : the pricing process was aborted by the pricer, there is no guarantee that the current LP solution is optimal
 */
SCIP_DECL_PRICERREDCOST(ObjPricerColoring::scip_redcost)
{
   SCIPdebugMsg(scip, "call scip_redcost ...\n");

   /* set result pointer, see above */
   *result = SCIP_SUCCESS;
   
   /* call pricing routine */
   
   if(M->separatingPhase ||  isSolInteger(scip) || inhibit_alternate){
	
	auto start = chrono::steady_clock::now();
   	coloring_pricing(scip);
	auto end = chrono::steady_clock::now();
	M->_data.pricingTime+= chrono::duration_cast<chrono::milliseconds>(end - start).count();
	
   }
	else{
	//cout<<"yoyo"<<endl;
	M->separatingPhase = separation;
  	//*result = SCIP_DIDNOTRUN;
	}
   
   return SCIP_OKAY;
} 

bool ObjPricerColoring::isSolInteger(SCIP* scip){

	list<C_master_var*>::iterator it;
	//cout<<"Le point : "<<endl;
	for(it = M->L_var.begin(); it != M->L_var.end(); it++){

	//	cout<<SCIPgetSolVal(scip, NULL, (*it)->ptr)<<endl;
		if(! SCIPisIntegral(scip,SCIPgetSolVal(scip, NULL, (*it)->ptr))){
			return false;
			
		}
			
	
	}
	
	return true;

}



// Perform pricing
void ObjPricerColoring::coloring_pricing(
   SCIP*                 scip       /**< SCIP data structure */
   )
{
  #ifdef OUTPUT_PRICER
  cout<<"**************PRICER************"<<endl;
  #endif
  
  nbNewConstaints =nbNewConstaints + M->ctn.size();
 
  int i;
  vector<SCIP_Real> lambda;
  lambda.resize(G->nb_nodes);
  vector<SCIP_Real> rho;
  rho.resize(nbNewConstaints);
  
  for (i = 0; i < G->nb_nodes; i++)   // on récupère les coûts duaux des rhs1
    lambda[i]= rounding(scip, SCIPgetDualsolLinear(scip, M->V_node_ineq[i]));

  

   
 list<SCIP_ROW*>::iterator yy;
  i=0;
  for(yy = M->additional_ineq.begin();yy != M->additional_ineq.end();yy ++){ // on récupère les coûts duaux des contraintes ajoutées
	
	rho[i] = rounding(scip,SCIProwGetDualsol(*yy));
	i++;	
  }
  

    
  #ifdef OUTPUT_PRICER
  cout<<"dual solution"<<endl;
  for (i = 0; i < G->nb_nodes; i++){
    cout<<"Cte_"<<i<<"="<<lambda[i]<<" ";
  }
  cout<<endl;
  #endif

	  

    
      
  StableSet* stable = new StableSet();
  stable->initialize(G->nb_nodes);
  SCIP_Real cplex_objvalue;
 

  A_cplex.modify_pricer( M->ctn); // separated_constraints doit contenir la liste des contraintes générée depuis la dernière génération de colonnes

  M->ctn.clear();
	
  A_cplex.set_objective_coefficient(lambda,rho);
	
   vector<int>* coefficients = new vector<int>();  // contiendra la liste des coefficients du stable ajouté, dans les contraintes additionelles
  coefficients->resize(nbNewConstaints);
  

  


  bool stable_found = A_cplex.find_stableset(stable,cplex_objvalue,coefficients);

  


  //if((cplex_objvalue <=1+epsilon) && (cplex_objvalue >=1-epsilon)){
  //	cplex_objvalue = 1;
  //}
  //double reduced_cost = 1-cplex_objvalue;
  SCIP_Real reduced_cost = 1-cplex_objvalue;
  #ifdef OUTPUT_PRICERR
  if ( stable_found ){
    	    list<int>::const_iterator it;
	    cout<<"Stable found with reduced cost = "<<reduced_cost<<" : ";
	    for (it = stable->L_nodes.begin(); it != stable->L_nodes.end(); it++)
	      cout<<*it<<" ";
	    cout<<endl;
  }
  #endif
	
  /* add stable set variable */
 if ( stable_found && SCIPisNegative(scip, reduced_cost) && SCIPisSumLE(scip,1,cplex_objvalue) && (reduced_cost< -1e-3) ){
   //if(stable_found && reduced_cost < -1e-5 ) {
	    list<int>::const_iterator it;
	    M->_data.nbColonnes++;
	    M->separatingPhase --;
	    M->_data.objvalue.push_back( SCIPgetDualbound(scip));
	    M->_data.sep.push_back(M->separatingPhase);
	    M->_data.dualValue.push_back(0);
	    char var_name[255];
	    SCIPsnprintf(var_name, 255, "V_%d",M->L_var.size());
	    SCIPdebugMsg(scip, "new variable <%s>\n", var_name);

	    /* create the new variable: Use upper bound of infinity such that we do not have to care about
	     * the reduced costs of the variable in the pricing. The upper bound of 1 is implicitly satisfied
	     * due to the set partitioning constraints.
	     */

	    //ofstream res;
   	    //res.open("../stableSets/stable.txt",std::ios_base::app);
          


	    C_master_var *var = new C_master_var;
	    SCIPcreateVar(scip, &(var->ptr), var_name,
			  0.0,                     // lower bound
			  1.0,      // upper bound
			  1.0,                     // objective
			  SCIP_VARTYPE_INTEGER, // variable type
			  false, false, NULL, NULL, NULL, NULL, NULL);

	    /* add new variable to the list of variables to price into LP (score: leave 1 here) */
	    SCIPaddPricedVar(scip, var->ptr, 1.0);

	    /* add coefficient into the set covering constraints */
	    var->stable=stable;
		i = 0;

	
	    for (it = stable->L_nodes.begin(); it != stable->L_nodes.end(); it++){
	      SCIPaddCoefLinear(scip, M->V_node_ineq[*it], var->ptr, 1.0);
	
		i++;
	      
	    }
	
	    i = 0;
	    //res<<"nb constraints : "<<nbNewConstaints<<endl;
	    for(yy = M->additional_ineq.begin() ; yy != M->additional_ineq.end(); yy++){ // add coefficients of additional constraints
		 SCIPaddVarToRow(scip, (*yy), var->ptr, (*coefficients)[i]);
		//res<<(*coefficients)[i]<<" : " <<rho[i] << " ";
		
		i++;	

	    }
		//res<<endl;
		//res<<reduced_cost<<endl;
		
		
	    M->L_var.push_back(var);
     	    //res<<"\n"<<endl;
	    //res.close();


	  }
	  

#ifdef OUTPUT_PRICER
  SCIPwriteTransProblem(scip, "coloring.lp", "lp", FALSE);
  cout<<"************END PRICER******************"<<endl;
#endif


}


