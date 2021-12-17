#ifndef ODDHOLE
#define ODDHOLE

#include <ilcplex/ilocplex.h>
#include <vector>
#include <list>

//#include "constraints.h"
#include "Graph.h"


#include "Find_Stable_set_with_Cplex.h"
#include "objscip/objscip.h"
using namespace std;


class Cplex_pricing_algo;

class Constraint{ //abstract class wich has a method to change the pricer every time we add inequalities in master problem

	public : 
	virtual IloNumVar* modifyPricer(Cplex_pricing_algo* pricer) = 0;
	virtual ~Constraint(){}
};


class OddHoleCst : public Constraint {

	public : 
	unsigned int rhs;

	list<int>* c;
	
	OddHoleCst(list<int>* c_){
		c = c_;
		
	

	}
	IloNumVar* modifyPricer(Cplex_pricing_algo* pricer);



};



class ChvatalGomoryCst : public Constraint { 
	public : 
	unsigned int nb_nodes;
	unsigned int rhs;
	vector<SCIP_Real>* u;
	
	ChvatalGomoryCst(vector<SCIP_Real>* u_,int nb_nodes_){   //give to this constructor, vectors with weigth with less or equal than two decimals
	
		u = u_;
		nb_nodes = nb_nodes_;
		SCIP_Real sum=0;
		for(int i = 0 ; i < u->size(); i ++){
			sum+= (*u)[i];

		}
		rhs = ceil(sum-0.000001);
	}
	IloNumVar* modifyPricer(Cplex_pricing_algo* pricer);  


};

#endif
