
#include "Constraint.h"




int nbCstChavtal = 0;
int nbCstOddCycle = 0;

//////////// Constraints definition ///////////////////

IloNumVar* OddHoleCst::modifyPricer(Cplex_pricing_algo* pricer){


	IloNumVar* vc =new IloNumVar(pricer->env, 0.0, 1.0, ILOINT);
	ostringstream varname;	   
	varname.str("");
    	varname<<"OddHole_"<<nbCstOddCycle;
	nbCstOddCycle++;
        vc->setName(varname.str().c_str());
	IloRangeArray CC(pricer->env);
		
	list<int>::iterator it;	


	
	IloExpr* cst = new IloExpr(pricer->env);
	(*cst)+= (*vc);
	for (it=c->begin(); it !=   c->end();it++){

      		
      		(*cst)+=  -pricer->x[(*it)];
      		
      	
  	}
	
  	CC.add((*cst)<=0);
		  
  	pricer->model.add(CC);
	

	return vc;


}

IloNumVar* ChvatalGomoryCst::modifyPricer(Cplex_pricing_algo* pricer){
	IloNumVar* vc =new IloNumVar(pricer->env, 0.0, rhs, ILOINT);
	ostringstream varname;	   
	varname.str("");
	

    	varname<<"ChvatalGomory_"<<nbCstChavtal;
	nbCstChavtal++;
	vc->setName(varname.str().c_str());
	IloRangeArray CC(pricer->env);


	 //ofstream res;
   	 //res.open("../wholeLP/constraints.txt",std::ios_base::app);

	IloExpr* cst = new IloExpr(pricer->env);
	(*cst)+= (*vc);
	for (int i = 0; i <   nb_nodes ;i++){

      		
      		(*cst)+=  -pricer->x[i] * (*u)[i]; 
      	//	res<<(*u)[i]<<" ";	
      	
  	}

	CC.add((*cst)<=1-0.0001);
	//res<<"\n"<<"";
  	pricer->model.add(CC);
	
	//res.close();
	return vc;

	


}
