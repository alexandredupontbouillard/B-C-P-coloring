#include "OddHoleSeparator.h"
#include <unistd.h>


using namespace std;

void OddHoleSeparator::init_handler(){ //  init_handler(list<SCIP_VAR *> &variables_param)

	return;
}


/** creates and captures a constraint */
SCIP_RETCODE SCIPcreate(
				   SCIP*                 scip,               /**< SCIP data structure */
				   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
				   const char*           name,               /**< name of constraint */
				   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
				   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing? */
				   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing? */
				   SCIP_Bool             check,              /**< should the constraint be checked for feasibility? */
				   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing? */
				   SCIP_Bool             local,              /**< is constraint only valid locally? */
				   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
				   SCIP_Bool             dynamic,            /**< is constraint dynamic? */
				   SCIP_Bool             removable           /**< should the constraint be removed from the LP due to aging or cleanup? */
				   )
{ 
	/* find the constraint handler */
	SCIP_CONSHDLR* conshdlr=NULL;
	conshdlr = SCIPfindConshdlr(scip, "odd Hole");
	if( conshdlr == NULL ){
		SCIPerrorMessage("J'ai pas trouve le handler\n");
		return SCIP_PLUGINNOTFOUND;
	} 

	/* create constraint */
	//SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, NULL, initial, separate, enforce, check, propagate,
 	//		    local, modifiable, dynamic, removable, FALSE) );
	cout<<"CREATE"<<endl;
	return SCIP_OKAY;
}



///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
/// Check if an integer solution is a solution of the problem
SCIP_RETCODE OddHoleSeparator::scip_check(SCIP* scip , SCIP_CONSHDLR* cshdlr, SCIP_CONS** cons, int x, SCIP_SOL* sol, unsigned int y, unsigned int w, unsigned int z, unsigned int h, SCIP_RESULT* result){

	
//SCIP_RETCODE OddHoleSeparator::scip_check(
//				    SCIP*              scip,               /**< SCIP data structure */
//				    SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
//				    SCIP_CONS**        conss,              /**< array of constraints to process */
//				    int                nconss,             /**< number of constraints to process */
//				    SCIP_SOL*          sol,                /**< the solution to check feasibility for */
//				    SCIP_Bool          checkintegrality,   /**< has integrality to be checked? */
//				    SCIP_Bool          checklprows,        /**< have current LP rows to be checked? */
//				    SCIP_Bool          printreason,        /**< should the reason for the violation be printed? */
//				    SCIP_RESULT*       result              /**< pointer to store the result of the feasibility checking call */
//				    ){
	/*
	cout<<"DEBUT CHECK\n";
	list<SCIP_VAR*> var_W;
	list<int>::iterator it_rhs;
	int objectif;

	//Pour vérifier la faisabilité, on regarde le coefficient de chacune des variables de chaque contrainte W et on le compare au rhs
	*result=SCIP_FEASIBLE;
	it_rhs = constraint_W_gamma->begin();
	for(list< list<int> >::iterator it_cons_sommets = constraint_W_sommets->begin() ; it_cons_sommets != constraint_W_sommets->end() ; it_cons_sommets++){
		var_W.clear();
		recup_variables_W(*it_cons_sommets,var_W);
		objectif = 0;
		list< list<int> >::iterator it_som = variables_sommets->begin();
		for(list<SCIP_VAR*>::iterator it_var = var_W.begin() ; it_var != var_W.end() ; it_var++){
			if(SCIPgetSolVal(scip, sol, *it_var) > 0.999){
				objectif = objectif + (ceil((*it_cons_sommets).size()/(*it_rhs)) - ceil((*it_cons_sommets).size()/(max(intersection(*it_cons_sommets,*it_som),*it_rhs))))+1;
			}
			it_som++;
		}
		if(objectif < *it_rhs){//Si une contrainte est violée par la solution, la solution est infaisable
			*result=SCIP_INFEASIBLE;
			cout<<"FIN CHECK\n";
			return SCIP_OKAY;
		}
		it_rhs++;
	}
	cout<<"FIN CHECK\n";*/
        
	// SCIP_CALL( separ_Contrainte(scip, NULL, result,conshdlr ) );
	*result = SCIP_INFEASIBLE;
	cout<<"FIN CHECK\n";
	return SCIP_OKAY;

}


SCIP_RETCODE OddHoleSeparator::scip_lock(SCIP*, SCIP_CONSHDLR*, SCIP_CONS*, SCIP_LOCKTYPE, int, int){



	return SCIP_OKAY;
}


SCIP_RETCODE OddHoleSeparator::separ_Contrainte(
		      SCIP*              scip,               /**< SCIP data structure */
		      SCIP_SOL*          sol,                /**< primal solution that should be separated */
		      SCIP_RESULT*       result,             /**< pointer to store the result of the separation call */
		      SCIP_CONSHDLR*     conshdlr          /**< the constraint handler itself */
		      )
{	

	vector<StableSet*> stables = vector<StableSet*>(master->L_var.size());
	vector<float> cost =  vector<float>(master->L_var.size());;
	
	list<C_master_var*>::iterator it;  // on récupère les valeurs des variables et les stables en base 
	int i = 0;
	for(it = master->L_var.begin(); it != master->L_var.end(); it++){
		stables[i] = (*it)->stable;
		cost[i]=(*it)->ptr->relaxsol;

		i++;
	}

	list<list<int>*> violatedCycles  ; // 
	
	G->findViolatedCycleIq(&stables,cost,violatedCycles,3);
	
	cout<< " Number of violated cycles : "<<violatedCycles.size()<<endl;
	list<list<int>*>::iterator itt;
	list<int>::iterator ittt;
	list<C_master_var*>::iterator ite; 
	for(itt = violatedCycles.begin(); itt != violatedCycles.end(); itt++){
		
		ostringstream namebuff; 
		namebuff.str(""); 
		namebuff << "cycleCst";
		namebuff << oddHole->size();


		SCIP_ROW* row;
            	SCIP_CALL( SCIPcreateEmptyRowConshdlr(scip, &row, conshdlr, "oddHole", 3.0, SCIPinfinity(scip), FALSE, TRUE, FALSE) );

            	SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
		

		SCIP_CONS *tmp_cons;
		SCIP_CALL_EXC( SCIPcreateConsLinear(scip, &tmp_cons, namebuff.str().c_str(), 0, NULL, NULL,
											3,//lhs
											SCIPinfinity(scip),//rhs
											TRUE, //initial
											TRUE, //separate
											TRUE, //enforce
											FALSE, //check
											TRUE, //propagate
											FALSE, //local
											TRUE, //modifable
											FALSE, //dynamic
											FALSE, //removable
											FALSE) );//stickinganode

		
		//rhsOddHole->push_back(3);
		//oddHole->push_back(tmp_cons);
		master->additional_ineq.push_back(tmp_cons);
	
		for(ite = master->L_var.begin(); ite != master->L_var.end(); ite++){  // tous les stables intersectant le cycle auront un coeffcient 1 dans la contrainte
			for(ittt = (*itt)->begin() ; ittt != (*itt)->end() ; ittt ++){
				
				if((*ite)->stable->contains(*ittt) ){

					SCIP_CALL_EXC( SCIPaddCoefLinear(scip, tmp_cons, (*ite)->ptr, 1) );
					SCIP_CALL( SCIPaddVarToRow(scip, row,  (*ite)->ptr, 1.0) );
				}

			}


		}
		
		if( SCIPisCutEfficacious(scip, sol, row) )
			    {
			       SCIP_Bool infeasible;
			       SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );
			       if ( infeasible )
				  *result = SCIP_CUTOFF;
			       else
				  *result = SCIP_SEPARATED;
			    }
			    //SCIP_CALL( SCIPreleaseRow(scip, &row) );

		SCIP_CALL( SCIPcreateOddHole(scip, &tmp_cons, "oddHOlecst", NULL, TRUE, TRUE, TRUE, FALSE, TRUE,
 			    FALSE, TRUE, FALSE) );
		
		//SCIPaddPoolCut(scip,row);
		OddHoleCst newcst = OddHoleCst(*itt);
		master->ctn.push_back(&newcst);
		*result = SCIP_SEPARATED; //SCIP_CONSADDED
		

		

	

	}

	
	return SCIP_OKAY;
}





/** separation method of constraint handler for LP solution */
SCIP_RETCODE OddHoleSeparator::scip_sepalp(
				     SCIP*              scip,               /**< SCIP data structure */
				     SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
				     SCIP_CONS**        conss,              /**< array of constraints to process */
				     int                nconss,             /**< number of constraints to process */
				     int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
				     SCIP_RESULT*       result              /**< pointer to store the result of the separation call */
				     )
{
  cout << "SEPALP"<<endl;
   	
  SCIP_CALL( separ_Contrainte(scip, NULL, result,conshdlr ) );
	
  cout<<"END SEPALP"<<endl;
  return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solution 
    The method is called outside the LP solution loop (e.g., by a relaxator or a primal heuristic), 
    which means that there is no valid LP solution. Instead, the method should produce cuts that 
    separate the given solution. **/
SCIP_RETCODE OddHoleSeparator::scip_sepasol(
				      SCIP*              scip,               /**< SCIP data structure */
				      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
				      SCIP_CONS**        conss,              /**< array of constraints to process */
				      int                nconss,             /**< number of constraints to process */
				      int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
				      SCIP_SOL*          sol,                /**< primal solution that should be separated */
				      SCIP_RESULT*       result              /**< pointer to store the result of the separation call */
				      )
{
  cout << "SEPASOL";
	
  //SCIP_CALL( separ_Contrainte(scip, NULL, result,conshdlr ) );
cout << "SEPASOL END";
  return SCIP_OKAY;
}


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/// ENFOLP is launched when there is no constraint left to had
/// The method is called at the end of the node processing loop for a node where the LP was solved. 
/// The LP solution has to be checked for feasibility. If possible, an infeasibility should be resolved 
/// by branching, reducing a variable's domain to exclude the solution or separating the solution with 
/// a valid cutting plane.

SCIP_RETCODE OddHoleSeparator::scip_enfolp(
				     SCIP*              scip,               /**< SCIP data structure */
				     SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
				     SCIP_CONS**        conss,              /**< array of constraints to process */
				     int                nconss,             /**< number of constraints to process */
				     int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
				     SCIP_Bool          solinfeasible,      /**< was the solution already declared infeasible by a constraint handler? */
				     SCIP_RESULT*       result              /**< pointer to store the result of the enforcing call */
				     ){
 


//  bool entier;
//  int k;
//  double tmp;
//  double val;

//  entier=true;
//  k=0;
//  while ((entier)&&(k<G->nbnode)){
//    tmp=SCIPgetSolVal(scip, NULL, (*v_varnode)[k]);
//    if ((tmp<=0.9999)&&(tmp>=0.00001))
//      entier=false;
//    k++;
//  }

//  if (entier){

//    SCIP_SOL* newsol;
//    SCIP_Bool success;

//    SCIP_CALL( SCIPcreateSol (scip, &newsol, NULL) );

//    for (k=0;k<G->nbnode;k++){
//      val= SCIPgetSolVal(scip, NULL, (*v_varnode)[k]);
//      SCIP_CALL( SCIPsetSolVal(scip, newsol, (*v_varnode)[k], val) );
//    }

//    SCIP_CALL( SCIPtrySol(scip, newsol, FALSE, FALSE, FALSE, FALSE, &success) );

//    *result=SCIP_CUTOFF;
//  }
//  else{
//    *result = SCIP_BRANCHED;
//  }
cout << "ENFOLP\n";
	

cout << "ENFOLP END\n";
	return SCIP_OKAY;
	
}

/*** The method is called at the end of the node processing loop for a node where the LP was not solved.
     The pseudo solution has to be checked for feasibility. If possible, an infeasibility should be 
     resolved by branching, reducing a variable's domain to exclude the solution or adding an additional 
     constraint. Separation is not possible, since the LP is not processed at the current node.  ***/

SCIP_RETCODE OddHoleSeparator::scip_enfops(
				     SCIP*              scip,               /**< SCIP data structure */
				     SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
				     SCIP_CONS**        conss,              /**< array of constraints to process */
				     int                nconss,             /**< number of constraints to process */
				     int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
				     SCIP_Bool          solinfeasible,      /**< was the solution already declared infeasible by a constraint handler? */
				     SCIP_Bool          objinfeasible,      /**< is the solution infeasible anyway due to violating lower objective bound? */
				     SCIP_RESULT*       result              /**< pointer to store the result of the enforcing call */
				     ){
	
  *result = SCIP_FEASIBLE; //indicates that all the constraints of the handler are feasible

cout << "ENFOPs\n";
	

cout << "ENFOPS END\n";

  return SCIP_OKAY;
}


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
// This method is called, after a constraint is added or removed from the transformed problem.
// If the constraint may get violated by changing the variable in any direction, it should call SCIPaddVarLocks(scip, var, nlockspos + nlocksneg, nlockspos + nlocksneg)

SCIP_RETCODE OddHoleSeparator::scip_lock(
				   SCIP*              scip,               /**< SCIP data structure */
				   SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
				   SCIP_CONS*         cons,               /**< the constraint that should lock rounding of its variables, or NULL if the
									   *   constraint handler does not need constraints */
				   int                nlockspos,          /**< no. of times, the roundings should be locked for the constraint */
				   int                nlocksneg           /**< no. of times, the roundings should be locked for the constraint's negation */
				   ){
cout << "LOCK\n";

// int i;


// for (i=0;i<G->nbnode;i++){
// 
//   SCIP_CALL( SCIPaddVarLocks(scip, (*v_varnode)[i] , nlocksneg, nlockspos)); //If the constraint may become violated by increasing the value of a variable,

//   //nlockspos, nlocksneg));  If the constraint may become violated by decreasing the value of a variable
//   //nlockspos + nlocksneg, nlockspos + nlocksneg)); If the constraint may become violated by changing the variable in any direction, 

// }
	
  return SCIP_OKAY;
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/// ??

/** transforms constraint data into data belonging to the transformed problem */
SCIP_RETCODE OddHoleSeparator::scip_trans(
				    SCIP*              scip,               //**< SCIP data structure *
				    SCIP_CONSHDLR*     conshdlr,           //**< the constraint handler itself *
				    SCIP_CONS*         sourcecons,         //**< source constraint to transform *
				    SCIP_CONS**        targetcons          //**< pointer to store created target constraint *
				    )
{
	
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, NULL,
 		    SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
			    SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
			    SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons), 
			    SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

 cout << "TRANS??\n";
	
  return SCIP_OKAY;
}



SCIP_RETCODE SCIPcreateOddHole(
				   SCIP*                 scip,               /**< SCIP data structure */
				   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
				   const char*           name,               /**< name of constraint */
				   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
				   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing? */
				   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing? */
				   SCIP_Bool             check,              /**< should the constraint be checked for feasibility? */
				   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing? */
				   SCIP_Bool             local,              /**< is constraint only valid locally? */
				   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
				   SCIP_Bool             dynamic,            /**< is constraint dynamic? */
				   SCIP_Bool             removable           /**< should the constraint be removed from the LP due to aging or cleanup? */
				   )
{ 
	/* find the constraint handler */
	SCIP_CONSHDLR* conshdlr=NULL;
	conshdlr = SCIPfindConshdlr(scip, "OddHole");
	if( conshdlr == NULL ){
		SCIPerrorMessage("J'ai pas trouve le handler de W\n");
		return SCIP_PLUGINNOTFOUND;
	} 
	

	
	SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, NULL, initial, separate, enforce, check, propagate,
 			    local, modifiable, dynamic, removable, FALSE) );
	cout<<"CREATE ODD HOLE"<<endl;
	return SCIP_OKAY;
}

