#ifndef FIND

#define FIND

#include <ilcplex/ilocplex.h>
#include <vector>
#include <list>
#include "objscip/objscip.h"
//#include "constraints.h"
#include "Graph.h"
#include "Constraint.h"
using namespace std;



class Constraint;


class Cplex_pricing_algo{
 public:
   unsigned int nbIt = 0;
   int nb_additional_variable=0;
  C_Graph *G;
  
  IloEnv   env;
  IloModel model;
  IloObjective obj;
  IloCplex cplex;
      
  IloNumVarArray x;

  IloNumVarArray additional;   // la liste des variables ajoutées dans la fonction objectif


 


  list<list<int> > *L_notequal; // Pointeur on a list of stable sets that
                                // have not to be found as solutions



  void create_MIP(C_Graph *G);
  
  

  void set_objective_coefficient(const vector<SCIP_Real>& obj_coeff, const vector<SCIP_Real>& additional_coeff);

  void modify_pricer(list<Constraint*> constraint);
  
  // Launch Cplex solver and get back an optimal stable set 
  bool find_stableset(StableSet* stable, SCIP_Real& objvalue, vector<int>* coefficients);


};

#endif
