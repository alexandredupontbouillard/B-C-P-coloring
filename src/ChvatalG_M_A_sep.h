
#ifndef ChvatalG_M_A_sep
#define ChvatalG_M_A_sep

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
#include <algorithm>
#include "technique.h"
class C_master_coloring;

class ChvatalGMAsep : public scip::ObjSepa{
	public:
	/************************************************************/
	/*******************STRUCTURE DU SEPARATOR*******************/
	/************************************************************/

	int nb_nodes;
	int nbColumnAtLastIteration=0;
	int nbIt=0;
	list<SCIP_ROW *> *chvatalG;//  list storing pointer of each added cut
	
	list<int>* rhs;//list storing rhs of each added cut
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
        vector<pair<vector<SCIP_Real>*,int>*>*  poolConstraints;


	IloConstraint precision;


	//Destructor
	virtual ~ChvatalGMAsep(){
	};


	
	void create_MIP();
	int separ_Contrainte(
		      SCIP*              scip,               /**< SCIP data structure */
		      SCIP_RESULT*       result            /**< pointer to store the result of the separation call */
		      );

	SCIP_RETCODE scip_execlp(SCIP*, SCIP_SEPA*, SCIP_RESULT*, unsigned int);
	SCIP_RETCODE scip_execsol(SCIP*, SCIP_SEPA*, SCIP_RESULT*, unsigned int);
	SCIP_RETCODE scip_init(SCIP*,SCIP_SEPA*);
	int makeConstraints(  SCIP* scip,SCIP_RESULT* result, SCIP_SEPA* sepa );
	vector<pair<vector<SCIP_Real>*,int>*>* getBestConstraintsPopulate();
	//void addPoolConstraintsToModel(SCIP* scip, SCIP_SEPA* sepa, SCIP_RESULT* result, unsigned int i);
	SCIP_Real distanceInfinie(int sol1, int sol2);
	bool isViolatedWhenRounding(pair<vector<SCIP_Real>*,int>* ct, SCIP* scip);
	ChvatalGMAsep(SCIP* scip,C_master_coloring* _master, int nb_nodes_)
		: scip::ObjSepa (scip, "ChvatalGomory separator MA", "exact separation of first chvatal closure", 1500,1 , 1, FALSE, FALSE) // int sepapriority,  sepafrequency , maxdist 

	{
	master = _master;
	rhs = new list<int>();

	chvatalG = new list<SCIP_ROW *>() ;
	nbColumnAtLastIteration=0;
	nb_nodes= nb_nodes_;
	nbIt = 0;
	
	};
};	




#endif
