#ifndef ODDHOLE_H
#define ODDHOLE_H

#include "objscip/objscip.h"
#include "scip/scip.h"
#include "scip/cons_linear.h"
#include "scip_exception.hpp"

#include <list>
#include <algorithm>

#include "Graph.h"
#include "SCIP_master.h"
#include "Constraint.h"

using namespace std;

class C_master_coloring;

class OddHoleSeparator : public scip::ObjConshdlr{
	public:
	/************************************************************/
	/*******************STRUCTURE DU HANDLER*********************/
	/************************************************************/
	C_Graph* G;//Graphe de base de l'instance
	C_Graph Cg; // graph complementaire

	list<SCIP_VAR *> *variables;//Liste stockant les pointeurs sur les variables
	list<SCIP_CONS *> *oddHole;//Vecteur stockant les pointeurs sur les contraintes (les odd Holes ont rhs = 3 mais le rhs des anti hole dépend de la taille de l'anti hole)
	
	list<int>* rhsOddHole;//liste sotckant le rhs des contraintes
	C_master_coloring* master;


	/************************************************************/
	/************************CONSTRUCTEUR************************/
	/************************************************************/
	//Constructeur
	OddHoleSeparator(SCIP* scip, C_Graph* G_,C_master_coloring* _master	 )
	: scip::ObjConshdlr(scip, "OddHole", "Handler For odd hole Constraints",
						1000000, -2000000, -2000000, //	int sepapriority, int enfopriority, int checkpriority, 
						1, -1, 1, 0,//	int sepafreq, int propfreq, int eagerfreq, int maxprerounds,
						FALSE, TRUE, FALSE, SCIP_PROPTIMING_BEFORELP, SCIP_PRESOLTIMING_FAST)//	SCIP_Bool delaysepa, SCIP_Bool delayprop, SCIP_Bool delaypresol, SCIP_Bool needscons
	{G = G_;
	Cg = G->getComplementary();
	master = _master;
	rhsOddHole = new list<int>();
	variables = new list<SCIP_VAR *>();
	oddHole = new list<SCIP_CONS *>() ;	
	};
	//Destructeur
	virtual ~OddHoleSeparator(){
	};

	/************************************************************/
	/*************************FONCTIONS**************************/
	/************************************************************/

	void init_handler(); //void init_handler(list<SCIP_VAR *> &variables_param)

SCIP_RETCODE scip_check(SCIP*, SCIP_CONSHDLR*, SCIP_CONS**, int, SCIP_SOL*, unsigned int, unsigned int, unsigned int, unsigned int, SCIP_RESULT*);
SCIP_RETCODE scip_lock(SCIP*, SCIP_CONSHDLR*, SCIP_CONS*, SCIP_LOCKTYPE, int, int);
	SCIP_RETCODE scip_check(
				  SCIP*              scip,               /**< SCIP data structure */
				  SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
				  SCIP_CONS**        conss,              /**< array of constraints to process */
				  int                nconss,             /**< number of constraints to process */
				  SCIP_SOL*          sol,                /**< the solution to check feasibility for */
				  SCIP_Bool          checkintegrality,   /**< has integrality to be checked? */
				  SCIP_Bool          checklprows,        /**< have current LP rows to be checked? */
				  SCIP_Bool          printreason,        /**< should the reason for the violation be printed? */
				  SCIP_RESULT*       result              /**< pointer to store the result of the feasibility checking call */
				  );

	SCIP_RETCODE scip_enfolp(
				   SCIP*              scip,               /**< SCIP data structure */
				   SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
				   SCIP_CONS**        conss,              /**< array of constraints to process */
				   int                nconss,             /**< number of constraints to process */
				   int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
				   SCIP_Bool          solinfeasible,      /**< was the solution already declared infeasible by a constraint handler? */
				   SCIP_RESULT*       result              /**< pointer to store the result of the enforcing call */
				   );

	SCIP_RETCODE scip_enfops(
				   SCIP*              scip,               /**< SCIP data structure */
				   SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
				   SCIP_CONS**        conss,              /**< array of constraints to process */
				   int                nconss,             /**< number of constraints to process */
				   int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
				   SCIP_Bool          solinfeasible,      /**< was the solution already declared infeasible by a constraint handler? */
				   SCIP_Bool          objinfeasible,      /**< is the solution infeasible anyway due to violating lower objective bound? */
				   SCIP_RESULT*       result              /**< pointer to store the result of the enforcing call */
				   );
	SCIP_RETCODE scip_lock(
				 SCIP*              scip,               /**< SCIP data structure */
				 SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */

				 SCIP_CONS*         cons,               /**< the constraint that should lock rounding of its variables, or NULL if the
									 *   constraint handler does not need constraints */
				 int                nlockspos,          /**< no. of times, the roundings should be locked for the constraint */
				 int                nlocksneg           /**< no. of times, the roundings should be locked for the constraint's negation */
				 );

	/** transforms constraint data into data belonging to the transformed problem */
	SCIP_RETCODE scip_trans(
				  SCIP*              scip,               //**< SCIP data structure *
				  SCIP_CONSHDLR*     conshdlr,           //**< the constraint handler itself *
				  SCIP_CONS*         sourcecons,         //**< source constraint to transform *
				  SCIP_CONS**        targetcons          //**< pointer to store created target constraint *
				  );

	/** separation method of constraint handler for LP solution
	*  possible return values for *result (if more than one applies, the first in the list should be used):
	*  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
	*  - SCIP_CONSADDED  : an additional constraint was generated
	*  - SCIP_REDUCEDDOM : a variable's domain was reduced
	*  - SCIP_SEPARATED  : a cutting plane was generated
	*  - SCIP_DIDNOTFIND : the separator searched, but did not find domain reductions, cutting planes, or cut constraints
	*  - SCIP_DIDNOTRUN  : the separator was skipped
	*  - SCIP_DELAYED    : the separator was skipped, but should be called again
	*/

	SCIP_RETCODE scip_sepalp(
				   SCIP*              scip,               /**< SCIP data structure */
				   SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
				   SCIP_CONS**        conss,              /**< array of constraints to process */
				   int                nconss,             /**< number of constraints to process */
				   int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
				   SCIP_RESULT*       result              /**< pointer to store the result of the separation call */
				   );

	/** separation method of constraint handler for arbitrary primal solution
	*  possible return values for *result (if more than one applies, the first in the list should be used):
	*  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
	*  - SCIP_CONSADDED  : an additional constraint was generated
	*  - SCIP_REDUCEDDOM : a variable's domain was reduced
	*  - SCIP_SEPARATED  : a cutting plane was generated
	*  - SCIP_DIDNOTFIND : the separator searched, but did not find domain reductions, cutting planes, or cut constraints
	*  - SCIP_DIDNOTRUN  : the separator was skipped
	*  - SCIP_DELAYED    : the separator was skipped, but should be called again
	*/
	SCIP_RETCODE scip_sepasol(
					SCIP*              scip,               /**< SCIP data structure */
					SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
					SCIP_CONS**        conss,              /**< array of constraints to process */
					int                nconss,             /**< number of constraints to process */
					int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
					SCIP_SOL*          sol,                /**< primal solution that should be separated */
					SCIP_RESULT*       result              /**< pointer to store the result of the separation call */
					);

	SCIP_RETCODE separ_Contrainte(
		      SCIP*              scip,               /**< SCIP data structure */
		      SCIP_SOL*          sol,                /**< primal solution that should be separated */
		      SCIP_RESULT*       result,             /**< pointer to store the result of the separation call */
		      SCIP_CONSHDLR*     conshdlr 
		      );



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
				  );
};

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
				   );

	
#endif
