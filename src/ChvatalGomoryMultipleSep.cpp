#include "Constraint.h"

#include "ChvatalGomoryMultipleSep.h"

#define DELTA 0.01

#define eps 1e-5

#define MAXNBCONSTRAINT 5
#define SIZEPOOL 20








void ChvatalGomoryMultipleSep::create_MIP(){

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
  cplex.setParam(IloCplex::Param::MIP::Pool::Replace,1 );  // permet de demander des solutions différentes
     cplex.setOut(env.getNullStream());
     cplex.setWarning(env.getNullStream());

	 cout<<"**************************** FIN INIT CPLEX *********************"<<endl;


}
vector<pair<vector<SCIP_Real>*,int>*>* ChvatalGomoryMultipleSep::getBestConstraintsPopulate(){
	

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



SCIP_Real ChvatalGomoryMultipleSep::distanceInfinie(int sol1, int sol2){
	SCIP_Real result= 0;
	for(int i = 0 ; i < nb_nodes; i ++){
		
		result+= abs(cplex.getValue(u[i], sol1), cplex.getValue(u[i], sol2));
		
	
	}

	return result;

}


int ChvatalGomoryMultipleSep::makeConstraints(  SCIP* scip,SCIP_RESULT* result,SCIP_SEPA* sepa){
	SCIP_Real sum_u = 0;
	vector<SCIP_Real>* newct = new vector<SCIP_Real>();
	newct->resize(nb_nodes);
	SCIP_Real c;

	vector<pair<vector<SCIP_Real>*,int>*>* cts =  getBestConstraintsPopulate();
	
	//cout<<cts->size()<<endl;
	for(int i = 0 ; i < cts->size(); i ++){

		SCIP_ROW* row;
    		SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &row, sepa, "ChvatalGomory",(*cts)[i]->second, SCIPinfinity(scip) , FALSE, TRUE, TRUE) );

    		SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
		list<C_master_var*>::iterator ite;
		list<int>::iterator stabit; 
		
		for(ite = master->L_var.begin(); ite != master->L_var.end(); ite++){  
			sum_u = 0; 


			for(stabit = (*ite)->stable->L_nodes.begin();stabit != (*ite)->stable->L_nodes.end(); stabit++ ){
				sum_u =sum_u + (*(*cts)[i]->first)[(*stabit)];
			  
			}
		
		
		
		SCIP_CALL( SCIPaddVarToRow(scip, row,  (*ite)->ptr, arrondiSup(sum_u)) ); //seg fault

		}
		
		
		SCIP_CALL( SCIPflushRowExtensions(scip, row) );
				
		if( SCIPisCutEfficacious(scip, NULL, row) ) //|| objvalue >= 0.5)
		{
			
			
			master->_data.nbChvatal++;
			unsigned int b;
			if( SCIPaddCut(scip, NULL, row, TRUE,&b) ){

				
				*result = SCIP_SEPARATED;
				rhs->push_back((*cts)[i]->second);
				chvatalG->push_back(row);
				master->additional_ineq.push_back(row);
				
				//SCIPaddPoolCut(scip,row);
				ChvatalGomoryCst* newcst = new ChvatalGomoryCst((*cts)[i]->first,nb_nodes);
				master->ctn.push_back(newcst);
				
			}
		}
		
		
	

		

	}



	

	
	return 1;


}


SCIP_RETCODE ChvatalGomoryMultipleSep::separ_Contrainte(
		      SCIP*              scip,               /**< SCIP data structure */
		      SCIP_RESULT*       result,            /**< pointer to store the result of the separation call */
		      SCIP_SEPA* sepa
		      ){
	master->_data.objvalue.push_back( SCIPgetDualbound(scip));
	master->_data.sep.push_back(master->separatingPhase);
	master->_data.dualValue.push_back(0);
	*result =SCIP_DIDNOTFIND ;
	//cout<<"separation \n \n \n "<<endl;
	vector<SCIP_Real> coeff; // récupérer les coefficients des stables de la solution
	int nb_stables =master->L_var.size(); 
	coeff.resize(nb_stables);

	list<C_master_var*>::iterator it;  // on récupère les valeurs des variables et les stables en base 
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



		obj.setLinearCoef(alpha[i],- arrondiS(coeff[i])); // on met à jour la nouvelle fonction objectif


		i++;
	}
	nbColumnAtLastIteration =master->L_var.size();
	
	model.add(CC);
	// il faudrait ajouter la contrainte de borne
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

		makeConstraints(scip, result,sepa);
	}
	

	master->separatingPhase = 1;
	return SCIP_OKAY;

}



SCIP_RETCODE ChvatalGomoryMultipleSep::scip_execlp(SCIP* scip, SCIP_SEPA* sepa, SCIP_RESULT* result, unsigned int i){
	 *result = SCIP_DIDNOTFIND;
	
	auto start = chrono::steady_clock::now();
	//cout<<"\n***************début séparation de chvatal*************\n "<<endl;
	SCIP_CALL( separ_Contrainte(scip,result,sepa));
	//cout<<"\n***************fin séparation de chvatal*************** \n" << endl;
	auto end = chrono::steady_clock::now();

	

	master->_data.separationTime += chrono::duration_cast<chrono::milliseconds>(end - start).count();
	return SCIP_OKAY;
}
SCIP_RETCODE ChvatalGomoryMultipleSep::scip_execsol(SCIP* scip, SCIP_SEPA* sepa, SCIP_RESULT* result, unsigned int i){

	 *result = SCIP_DIDNOTFIND;
	//cout<<"\n***************début séparation de chvatal*************\n "<<endl;
	SCIP_CALL( separ_Contrainte(scip,result,sepa));
	//cout<<"\n***************fin séparation de chvatal*************** \n" << endl;
	return SCIP_OKAY;
}
SCIP_RETCODE ChvatalGomoryMultipleSep::scip_init(SCIP* scip,SCIP_SEPA* sepa){
	
	create_MIP();

	return SCIP_OKAY;
}
