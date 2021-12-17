#include "Constraint.h"

#include "ChvatalG_M_A_sep.h"

#define DELTA 0.01

#define eps 1e-5

#define MAXNBCONSTRAINT 5
#define SIZEPOOL 20








void ChvatalGMAsep::create_MIP(){

  cout<<"**************************** INIT CPLEX *******************"<<endl;

  model= IloModel(env);

  nbIt = 0;


  u = IloNumVarArray(env, nb_nodes, 0.0, 1-0.0001  , ILOFLOAT);
  //big_u = IloNumVarArray(env, nb_nodes, 0, 100000  , ILOINT);
  alpha_0 = IloNumVar(env, 0.0, nb_nodes , ILOINT);
  alpha_0.setName("alpha_0");
  f_0= IloNumVar(env, 0.0, 1.0-DELTA, ILOFLOAT);
  f_0.setName("f_0");
  f = IloNumVarArray(env,0);
  alpha = IloNumVarArray(env,0);



  IloRangeArray CC(env);



  IloExpr cst(env);
  for (int  i =0; i <nb_nodes;i++){
      cst+= u[i];
    ostringstream varname;
 
    varname.str("");
    varname<<"u_"<<i;
    u[i].setName(varname.str().c_str());

    //IloExpr cst2(env);
    //cst2+= 100000* u[i]- big_u[i];
    //CC.add(cst2==0);
    
  }
  cst+= f_0 - alpha_0;
  CC.add(cst==0);

  model.add(CC);


  obj = IloAdd(model, IloMaximize(env, 0.0));


 obj.setLinearCoef(alpha_0,1);


  cplex = IloCplex(model);
  cplex.setParam(IloCplex::Param::MIP::Limits::Nodes, 1000);
  cplex.setParam(IloCplex::Param::MIP::Pool::Capacity, SIZEPOOL);
  cplex.setParam(IloCplex::Param::MIP::Pool::Replace,1 );  // fill solutionpool with bests solutions
     cplex.setOut(env.getNullStream());
     cplex.setWarning(env.getNullStream());

	 cout<<"**************************** FIN INIT CPLEX *********************"<<endl;


}
vector<pair<vector<SCIP_Real>*,int>*>* ChvatalGMAsep::getBestConstraintsPopulate(){
	

	vector<pair<vector<SCIP_Real>*,int>*>* result = new vector<pair<vector<SCIP_Real>*,int>*>();

	SCIP_Real c;
	list<int> solTaken;

	solTaken.push_back(cplex.IncumbentId);
	int lastCt  = cplex.IncumbentId;

	int maxCt = cplex.IncumbentId;
	SCIP_Real maxDiff;
	SCIP_Real r;


	for(int i = 0 ; i < MAXNBCONSTRAINT ; i ++){
		maxDiff = 0;
		


		for(int j= 0; j < cplex.getSolnPoolNsolns()  ; j++){
			
			if(! contains(solTaken, j)){
				r = distanceInfinie(lastCt,j);
				if( r > maxDiff){
					maxDiff = r;
					maxCt = j;
					
				}
				

			}


		}
		solTaken.push_back(maxCt);
		lastCt = maxCt;

	}
		
	list<int>::iterator it;
	for(it = solTaken.begin(); it != solTaken.end(); it++){
		SCIP_Real rhs = 0;
		vector<SCIP_Real>* v = new vector<SCIP_Real>(nb_nodes);
		for(int k = 0 ; k < nb_nodes ; k++){
			(*v)[k] = arrondis(cplex.getValue(u[k], (*it)));
			rhs+= arrondis(cplex.getValue(u[k], (*it)));

		}	
		pair<vector<SCIP_Real>*,int>* vv= new pair<vector<SCIP_Real>*,int>(v,arrondiSup(rhs));	
		result->push_back(vv);

	}	
	
	return result;

}	



SCIP_Real ChvatalGMAsep::distanceInfinie(int sol1, int sol2){
	SCIP_Real result= 0;
	for(int i = 0 ; i < nb_nodes; i ++){
		
		result+= abs(cplex.getValue(u[i], sol1), cplex.getValue(u[i], sol2));
		
	
	}

	return result;

}


int ChvatalGMAsep::makeConstraints(  SCIP* scip,SCIP_RESULT* result,SCIP_SEPA* sepa){
	
	SCIP_Real sum_u = 0;
	vector<SCIP_Real>* newct = new vector<SCIP_Real>();
	newct->resize(nb_nodes);
	SCIP_Real c;

	

	//cout<<cts->size()<<endl;
	for(int i = 0 ; i < poolConstraints->size(); i ++){

		SCIP_ROW* row;
    		SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &row, sepa, "ChvatalGomory",(*poolConstraints)[i]->second, SCIPinfinity(scip) , FALSE, TRUE, TRUE) );

    		SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
		list<C_master_var*>::iterator ite;
		list<int>::iterator stabit; 
		for(ite = master->L_var.begin(); ite != master->L_var.end(); ite++){  
			sum_u = 0; 


			for(stabit = (*ite)->stable->L_nodes.begin();stabit != (*ite)->stable->L_nodes.end(); stabit++ ){
				sum_u =sum_u + (*(	*poolConstraints)[i]->first)[(*stabit)];
			  
			}
		
		
		SCIP_CALL( SCIPaddVarToRow(scip, row,  (*ite)->ptr, arrondiSup(sum_u)) );

		}
		
		
		SCIP_CALL( SCIPflushRowExtensions(scip, row) );
				
		if( SCIPisCutEfficacious(scip, NULL, row) ) //|| objvalue >= 0.5)
		{
			cout<<"contrainte performante"<<endl;
			unsigned int b;
			if( SCIPaddCut(scip, NULL, row, TRUE,&b) ){
				
				//cout<<"cut added"<<endl;
				*result = SCIP_SEPARATED;
				rhs->push_back((*poolConstraints)[i]->second);
				chvatalG->push_back(row);
				master->additional_ineq.push_back(row);
				
				//SCIPaddPoolCut(scip,row);
				ChvatalGomoryCst* newcst = new ChvatalGomoryCst((*poolConstraints)[i]->first,nb_nodes);
				master->ctn.push_back(newcst);
				
				master->_data.nbChvatal++;
			

			}else{

				//cout<<"cut not added "<<endl;

			}
			
		
			
			
		}else{
			cout<<"il est dit que la coupe est nulle"<<endl;

		}
		
		
	

		

	


	}


	
	return 1;


}

bool ChvatalGMAsep::isViolatedWhenRounding(pair<vector<SCIP_Real>*,int>* ct, SCIP* scip){

	SCIP_Real sum = 0;

	SCIP_Real coeff_sum = 0;
		
	list<C_master_var*>::iterator ite;

	for(ite = master->L_var.begin(); ite != master->L_var.end(); ite++){  
		coeff_sum = 0; 


		for(list<int>::iterator stabit = (*ite)->stable->L_nodes.begin();stabit != (*ite)->stable->L_nodes.end(); stabit++ ){
			coeff_sum =coeff_sum + (*(ct->first))[(*stabit)];
		  
		}
		
		sum = sum + SCIPgetSolVal(scip, NULL, (*ite)->ptr)*arrondiSup(coeff_sum);	
	

	}
	if(sum <= ct->second-0.1){

		//cout<<"violated"<<endl;

		return true;
	}else{
		//cout<<"non violated "<<endl;

		 return false;
	}

}



int ChvatalGMAsep::separ_Contrainte(
		      SCIP*              scip,               /**< SCIP data structure */
		      SCIP_RESULT*       result            /**< pointer to store the result of the separation call */
		      ){
		
	master->_data.objvalue.push_back( SCIPgetDualbound(scip));
	master->_data.sep.push_back(master->separatingPhase);
	master->_data.dualValue.push_back(0);
	*result =SCIP_DIDNOTFIND ;
	//cout<<"separation \n \n \n "<<endl;
	vector<SCIP_Real> coeff; // store coefficients of the restricted master problem solution
	int nb_stables =master->L_var.size(); 
	coeff.resize(nb_stables);

	list<C_master_var*>::iterator it;  
	int i = 0;
	
	IloRangeArray CC(env);
	
	for(it = master->L_var.begin(); it != master->L_var.end(); it++){
		
		coeff[i]=SCIPgetSolVal(scip, NULL, (*it)->ptr);		
		
		if(i >= nbColumnAtLastIteration){   // we have to create new variables associated to new stables and their associated constraints			
			f.add(IloNumVar(env, 0.0, 1.0- DELTA, ILOFLOAT));
	 		alpha.add(IloNumVar(env, 0.0, nb_nodes, ILOINT));
			ostringstream varname;
    			varname.str("");
    			varname<<"f_"<<i+1;
		        f[i].setName(varname.str().c_str());
			
    			varname.str("");
    			varname<<"alpha_"<<i+1;
		        alpha[i].setName(varname.str().c_str());
			
			IloExpr cst(env);
			
			cst+= f[i]-alpha[i];
			//cst+= -alpha[i];
			for(int j = 0 ; j < nb_nodes ; j++){
				if((*it)->stable->contains(j)){

					cst+=  u[j];
				}

			}
			CC.add(cst==0);
	

		}



		obj.setLinearCoef(alpha[i],- arrondiS(coeff[i])); // update objective function of the separator


		i++;
	}
	nbColumnAtLastIteration =master->L_var.size();
	
	model.add(CC);
	
	if ( !cplex.populate() ) {
    		
		cout<<"seems like model is infeasible"<<endl;
		return SCIP_OKAY;
  	}
	
	//ostringstream filename;
    	//filename.str("");
    	//filename<<"../debug/sep"<<nbIt<<".lp"<<"";
	//nbIt++;

	
	//double yo;

	double objvalue=cplex.getObjValue();
       	

	list<C_master_var*>::iterator ite;
	list<int>::iterator stabit; 
	SCIP_Real c;
	

	if(objvalue >=0.1  ){


		poolConstraints =  getBestConstraintsPopulate();
		if(isViolatedWhenRounding((*poolConstraints)[0],scip)){
			return 1;
		}
	}
	else{

		poolConstraints = new vector<pair<vector<SCIP_Real>*,int>*>();
	}
	

	return 0;

}

/*
void ChvatalGMAsep::addPoolConstraintsToModel(SCIP* scip, SCIP_SEPA* sepa, SCIP_RESULT* result, unsigned int i){

	if(poolConstraints->size() ==0){
		*result = SCIP_DIDNOTFIND;


	}
	else{
		list<pair<SCIP_ROW*,pair<vector<SCIP_Real>*,int>*>*>::iterator it;
		
		for(int i = 0 ; i< poolConstraints->size(); i++){
			
			master->_data.nbChvatal++;
			unsigned int b;
			if( SCIPaddCut(scip, NULL, (*poolConstraints)[i]->first, TRUE,&b) ){
				
				
				*result = SCIP_SEPARATED;
				rhs->push_back((*it)->second->second);
				chvatalG->push_back((*it)->first);
				master->additional_ineq.push_back((*it)->first);
				
				//SCIPaddPoolCut(scip,row);
				ChvatalGomoryCst* newcst = new ChvatalGomoryCst((*it)->second->first,nb_nodes);
				master->ctn.push_back(newcst);
				
			}

		}
		
		poolConstraints.clear();
		
		*result = SCIP_SEPARATED;
		
	}


}
*/


SCIP_RETCODE ChvatalGMAsep::scip_execlp(SCIP* scip, SCIP_SEPA* sepa, SCIP_RESULT* result, unsigned int i){
	 *result = SCIP_DIDNOTFIND;
	
	
	cout<<"\n***************début séparation de chvatal*************\n "<<endl;
	

		int j  = makeConstraints(  scip, result, sepa);
	
	cout<<"\n***************fin séparation de chvatal*************** \n" << endl;
	

	


	return SCIP_OKAY;
}
SCIP_RETCODE ChvatalGMAsep::scip_execsol(SCIP* scip, SCIP_SEPA* sepa, SCIP_RESULT* result, unsigned int i){

	 *result = SCIP_DIDNOTFIND;
	cout<<"\n***************début séparation de chvatal*************\n "<<endl;

		int j  = makeConstraints(  scip, result, sepa);
		
	cout<<"\n***************fin séparation de chvatal*************** \n" << endl;

	return SCIP_OKAY;
}
SCIP_RETCODE ChvatalGMAsep::scip_init(SCIP* scip,SCIP_SEPA* sepa){



	return SCIP_OKAY;
}
