

#include "Data.h"



void Data::writeData(){
	
	ofstream res;
	char nn[40]; // = new char[40];

	strcpy(nn, "../results_data/");
	strcat(nn,filename); 
	res.open(nn ,std::ios_base::out);
	
	list<float>::iterator v;
	list<float>::iterator dv = dualValue.begin();
	list<int>::iterator cORp = sep.begin();

	res<<"filename, nb colonnes, nb constraintes de chvatal , separation Time, pricingTime, totalTime  "<<endl;

	res<<filename << "," << nbColonnes << ","<< nbChvatal << "," <<  separationTime << ","<< pricingTime <<","<< totalTime <<endl; 

	res<<"objValue, dualValue, sep/pricing "<<endl;
	
	for(v = objvalue.begin() ; v != objvalue.end(); v++){
		
		res<< (*v) << ","<< (*dv) << "," << *cORp << endl;


		dv++;
		cORp++;
		

	}

	res.close();

}
