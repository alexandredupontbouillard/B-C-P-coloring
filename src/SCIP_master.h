#ifndef __SCIP_MASTER_H__
#define __SCIP_MASTER_H__


#include <vector>
#include <list>
#include "Graph.h"
#include "Constraint.h"

/* scip includes */
#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"

#include "Data.h"

using namespace std;
using namespace scip;


class C_master_var{
 public:
  
   // Keep a pointer on every variables of the Master program
  SCIP_VAR * ptr;

  // Keep the list of nodes of the stable set corresponding to the variable
  StableSet* stable;
  void printSolution( SCIP * scip , SCIP_SOL *sol);
    
};



class C_master_coloring{
 public:


  // Keep a pointer on every constraint of the Master program
  vector<SCIP_CONS*> V_node_ineq;


  // Keep informations on every variables of the Master program
  list<C_master_var*> L_var;


  //Keep a pointer on every additional constraint of the Master program
  list<SCIP_ROW*> additional_ineq;

  list<Constraint*> ctn;

  Data _data;

   int separatingPhase;

   
	
  void printSolution( SCIP * scip , SCIP_SOL *sol);

};


#endif
