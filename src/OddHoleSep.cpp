
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/*  This file was written by Giacomo Nannicini,                              */
/*    Copyright (C) 2012 Singapore University of Technology and Design       */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sepa_OddHoleSep.c
 * @brief  Gomory Mixed-Integer Cuts
 * @author Giacomo Nannicini
 * @author Marc Pfetsch
 *
 * This file implements a Gomory Mixed-Integer (OddHoleSep) cuts generator that reads cuts from the simplex tableau, applying
 * the textbook formula:
 * \f[
 *    \sum_{j \in J_I : f_j \leq f_0} f_j x_j + \sum_{j \in J_I : f_j > f_0} f_0 \frac{1-f_j}{1 - f_0} x_j +
 *    \sum_{j \in J_C : a_j \geq 0}   a_j x_j - \sum_{j \in J_C : a_j < 0} f_0  \frac{a_j}{1-f_0} x_j \geq f_0.
 * \f]
 * Here, \f$J_I\f$ and \f$J_C \subseteq \{1, \ldots, n\}\f$ are the indices of integer and continuous non basic
 * variables, respectively. The tableaux row is given by \f$a_j\f$ and its right hand side is \f$a_0\f$. The values
 * \f$f_j\f$ for \f$j = 0, \ldots, n\f$ denote the fractional values of the tableaux row and rhs, i.e., \f$f_j = a_j -
 * \lfloor a_j \rfloor\f$.
 *
 * Here is a brief description of the simplex tableau that we can expect from the SCIP LP interfaces:
 *
 * - Nonbasic columns can be at the lower or upper bound, or they can be nonbasic at zero if they are free. Nonbasic columns
 *   at the upper bound must be flipped. Nonbasic free variables at zero are currently untested in the cut generator,
 *   but they should be handled properly anyway.
 *
 * - Nonbasic rows can be at lower or upper bound, depending on whether the lower or upper bound of the row is
 *   attained. SCIP always adds slack/surplus variables with a coefficient of +1: the slack variable is nonnegative in
 *   case of a <= constraint, it is nonpositive in case of a >= or ranged constraint. Therefore, slack variables
 *   corresponding to >= or ranged constraints must be flipped if the row is at its lower bound. (Ranged constraints at
 *   the upper bound do not have to be flipped, because the variable is nonpositive.)
 *
 * Generated cuts are modified and their numerical properties are checked before being added to the LP relaxation.
 * Default parameters for cut modification and checking procedures are taken from the paper
 *
 * G. Cornuejols, F. Margot, and G. Nannicini:@n
 * On the safety of Gomory cut generators.@n
 * Mathematical Programming Computation 5, No. 4 (2013), pp. 345-395.
 *
 * In addition to the routines described in the paper above, here we additionally check the support of the cutting
 * plane.
 *
 * @todo Check whether it is worth rescaling the cut to have integral coefficients on integer variables. This may lead
 * to an integral slack variable, that has stronger cut coefficients in subsequent rounds.
 */


/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/pub_misc.h"
#include "OddHoleSep.h"



/*SCIP_RETCODE SCIP_DECL_SEPAEXECLP(SCIP* scip, SCIP_SEPA* sepa, SCIP_RESULT* result, SCIP_Bool allowlocal){

	cout<<"separation begins"<<endl;

	SCIP_RETCODE rr =  separ_Contrainte(scip,result,sepa);

	cout<<"separation ends"<< endl;
	return rr;
}*/
/*
SCIP_RETCODE x (SCIP* scip, SCIP_SEPA* sepa, SCIP_RESULT* result, SCIP_Bool allowlocal){

cout<<"\n\n\nseparation begins\n\n\n"<<endl;

	//SCIP_RETCODE rr =  separ_Contrainte(scip,result,sepa);

	cout<<"separation ends"<< endl;
	//return rr;
	return SCIP_OKAY;
}
*/

SCIP_RETCODE OddHoleSep::scip_execlp(SCIP* scip, SCIP_SEPA* sepa, SCIP_RESULT* result, unsigned int i){
	
	cout<<"separation begins execlp"<<endl;
	
	SCIP_RETCODE rr =  separ_Contrainte(scip,result,sepa);
	
	cout<<"separation ends"<< endl;
	return rr;
	//	return SCIP_OKAY;

}

SCIP_RETCODE OddHoleSep::scip_execsol(SCIP* scip, SCIP_SEPA* sepa, SCIP_RESULT* result, unsigned int i){
	
	cout<<"separation begins execsol"<<endl;

	//SCIP_RETCODE rr =  separ_Contrainte(scip,result,sepa);

	cout<<"separation ends"<< endl;
	//return rr;
	return SCIP_OKAY;

}

SCIP_RETCODE OddHoleSep::separ_Contrainte(
		      SCIP*              scip,               /**< SCIP data structure */
		      SCIP_RESULT*       result,            /**< pointer to store the result of the separation call */
		      SCIP_SEPA* sepa
		      )
{	
	cout<<"yo"<<endl;
	vector<StableSet*> stables = vector<StableSet*>(master->L_var.size());
	vector<float> cost =  vector<float>(master->L_var.size());;
	
	list<C_master_var*>::iterator it;  // on récupère les valeurs des variables et les stables en base 
	int i = 0;

	for(it = master->L_var.begin(); it != master->L_var.end(); it++){
		stables[i] = (*it)->stable;
		cost[i]=SCIPgetSolVal(scip, NULL, (*it)->ptr);		

		i++;
	}	

		list<list<int>*> violatedCycles  ;  
	
	G->findViolatedCycleIq(&stables,cost,violatedCycles,3);
	
	cout<< " Number of violated cycles : "<<violatedCycles.size()<<endl;
	list<list<int>*>::iterator itt;
	list<int>::iterator ittt;
	list<C_master_var*>::iterator ite; 
	*result = SCIP_DIDNOTFIND;  
	for(itt = violatedCycles.begin(); itt != violatedCycles.end(); itt++){
		
		ostringstream namebuff; 
		namebuff.str(""); 
		namebuff << "cycleCst";
		namebuff << oddHole->size();


		SCIP_ROW* row;
            	SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &row, sepa, "oddHole",3.0, SCIPinfinity(scip) , FALSE, TRUE, TRUE) );

            	SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
		

		

		for(ite = master->L_var.begin(); ite != master->L_var.end(); ite++){  // tous les stables intersectant le cycle auront un coeffcient 1 dans la contrainte

			for(ittt = (*itt)->begin() ; ittt != (*itt)->end() ; ittt ++){
				
				if((*ite)->stable->contains(*ittt) ){

					SCIP_CALL( SCIPaddVarToRow(scip, row,  (*ite)->ptr, 1.0) );
					ittt = (*itt)->end();
					ittt--;
				}
				

			}

		}
		SCIP_CALL( SCIPflushRowExtensions(scip, row) );
						SCIPprintRow(scip,row,NULL);
		if( SCIPisCutEfficacious(scip, NULL, row) )
		{

			unsigned int b;
			SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE,&b) );

        		*result = SCIP_SEPARATED;
			rhsOddHole->push_back(3);
			oddHole->push_back(row);
			master->additional_ineq.push_back(row);
			
			//SCIPaddPoolCut(scip,row);
			OddHoleCst* newcst = new OddHoleCst(*itt);
			master->ctn.push_back(newcst);
    		}
		else{
			cout<<"niqué !! "<<endl;
			master->did_Separation_stop = true;
		}


		
	
		

	

	

	}

	
	return SCIP_OKAY;
}
	 




