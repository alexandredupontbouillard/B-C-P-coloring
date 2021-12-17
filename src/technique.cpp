


#include "technique.h"


SCIP_Real arrondiS(SCIP_Real b){  
	if(b< 1e-6)
		return 0;
	//return b;	
	return round(b*1000000.0)/1000000.0;


}


double arrondis(double b){ //arondi pour les coefficients de chvatal à 1e-6 près
	if(b < 1e-5)
		return 0;
	//return b;
	return floor(b*10000.0)/10000.0;

}

int arrondiSup(SCIP_Real b){
	SCIP_Real a = b - 0.000001;
	if(a < 1e-4)
		return 0;
	return ceil(a);


}

double abs(double r1, double r2){

	if(r1>r2){
		return r1-r2;
	}
	else{
		return r2-r1;
	}
}


bool contains(list<int> &l, int x){
	return find(l.begin(), l.end(), x) != l.end();
}


bool isSolEntier(SCIP* scip, C_master_coloring* Master){

	list<C_master_var*>::iterator it;
	//cout<<"Le point : "<<endl;
	for(it = Master->L_var.begin(); it != Master->L_var.end(); it++){

	//	cout<<SCIPgetSolVal(scip, NULL, (*it)->ptr)<<endl;
		if(! SCIPisIntegral(scip,SCIPgetSolVal(scip, NULL, (*it)->ptr))){
			return false;
			
		}
			
	
	}
	
	return true;

}
