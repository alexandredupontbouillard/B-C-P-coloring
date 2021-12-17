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

#ifndef __SCIP_SEPA_OddHole
#define __SCIP_SEPA_OddHole

#include "objscip/objscip.h"
#include "scip/scip.h"
#include "objscip/objsepa.h"
#include "scip/cons_linear.h"
#include "scip_exception.hpp"

#include <list>
#include <algorithm>

#include "Graph.h"
#include "SCIP_master.h"
#include "Constraint.h"

class C_master_coloring;

class OddHoleSep : public scip::ObjSepa{
	public:
	/************************************************************/
	/*******************STRUCTURE DU SEPARATOR*******************/
	/************************************************************/
	C_Graph* G;//Graphe de base de l'instance
	C_Graph Cg; // graph complementaire
	int a = 0;
	//list<SCIP_VAR *> *variables;//Liste stockant les pointeurs sur les variables
	list<SCIP_ROW *> *oddHole;//Vecteur stockant les pointeurs sur les contraintes (les odd Holes ont rhs = 3 mais le rhs des anti hole dépend de la taille de l'anti hole)
	
	list<int>* rhsOddHole;//liste sotckant le rhs des contraintes
	C_master_coloring* master;
	OddHoleSep(SCIP* scip, C_Graph* G_,C_master_coloring* _master	 )
		: scip::ObjSepa (scip, "Odd Hole separator", "glout separation of cycle iq", 2000,1 , 1, FALSE, FALSE) // int sepapriority,  sepafrequency , maxdist 

	{G = G_;
	Cg = G->getComplementary();
	master = _master;
	rhsOddHole = new list<int>();

	oddHole = new list<SCIP_ROW *>() ;
	
		
	};
	//Destructeur
	virtual ~OddHoleSep(){
	};


	

	SCIP_RETCODE separ_Contrainte(
		      SCIP*              scip,               /**< SCIP data structure */
		      SCIP_RESULT*       result,            /**< pointer to store the result of the separation call */
		      SCIP_SEPA* sepa
		      );

	SCIP_RETCODE scip_execlp(SCIP*, SCIP_SEPA*, SCIP_RESULT*, unsigned int);
	SCIP_RETCODE scip_execsol(SCIP*, SCIP_SEPA*, SCIP_RESULT*, unsigned int);


	 	
};	




#endif
