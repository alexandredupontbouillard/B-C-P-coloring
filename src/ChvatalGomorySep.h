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

/**@file   OddHoleSep.h
 * @ingroup SEPARATORS
 * @brief  Gomory Mixed-Integer Cuts
 * @author Giacomo Nannicini
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SEPA_Chvatal_Gomory_Sep
#define __SCIP_SEPA_Chvatal_Gomory_Sep

#include "objscip/objscip.h"
#include "scip/scip.h"
#include "objscip/objsepa.h"
#include "scip/cons_linear.h"
#include "scip_exception.hpp"
#include <ilcplex/ilocplex.h>
#include <vector>
#include <chrono>
#include <list>
#include <algorithm>

#include "Graph.h"
#include "SCIP_master.h"
#include "Constraint.h"
#include "technique.h"

class C_master_coloring;

class ChvatalGomorySep : public scip::ObjSepa{
	public:
	/************************************************************/
	/*******************STRUCTURE DU SEPARATOR*******************/
	/************************************************************/

	
	int nb_nodes;
	int nbColumnAtLastIteration=0;
	int nbIt=0;
	list<SCIP_ROW *> *chvatalG;//Vecteur stockant les pointeurs sur les contraintes (les odd Holes ont rhs = 3 mais le rhs des anti hole dépend de la taille de l'anti hole)
	
	list<int>* rhs;//liste sotckant le rhs des contraintes
	C_master_coloring* master;

	IloEnv   env;
  	IloModel model;
  	IloObjective obj;
  	IloCplex cplex;
      
  	IloNumVarArray u;
	IloNumVarArray big_u;
	IloNumVarArray f;
	IloNumVarArray alpha;
	IloNumVar alpha_0;
	IloNumVar f_0;

	IloConstraint precision;
	

	ChvatalGomorySep(SCIP* scip,C_master_coloring* _master, int nb_nodes_	 )
		: scip::ObjSepa (scip, "Chvatal Gomory separator", "exact separation of first chvatal closure", 1500,1 , 1, FALSE, FALSE) // int sepapriority,  sepafrequency , maxdist 

	{
	master = _master;
	rhs = new list<int>();

	chvatalG = new list<SCIP_ROW *>() ;
	nbColumnAtLastIteration=0;
	nb_nodes= nb_nodes_;
	nbIt = 0;
	};
	//Destructeur
	virtual ~ChvatalGomorySep(){
	};


	
	void create_MIP();
	SCIP_RETCODE separ_Contrainte(
		      SCIP*              scip,               /**< SCIP data structure */
		      SCIP_RESULT*       result,            /**< pointer to store the result of the separation call */
		      SCIP_SEPA* sepa
		      );

	SCIP_RETCODE scip_execlp(SCIP*, SCIP_SEPA*, SCIP_RESULT*, unsigned int);
	SCIP_RETCODE scip_execsol(SCIP*, SCIP_SEPA*, SCIP_RESULT*, unsigned int);
	SCIP_RETCODE scip_init(SCIP*,SCIP_SEPA*);
	int makeConstraint(  SCIP* scip,SCIP_RESULT* result, SCIP_SEPA* sepa );

	 	
};	




#endif
