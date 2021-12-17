

#include <list>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>

using namespace std;



class Data{


	public : 
	// tracking values for sweet curves
	list<float> objvalue;
	list<float> dualValue;
	list<int> sep; // tracking if the value in position i corresponds to a cutting or a pricing phase, if value is negative or null it's pricing step, otherwise it's cutting

	// stats oer the size of LP
	
	int nbColonnes;
	int nbChvatal;


	// stats over the time resolution

	double separationTime; //seconds
	double pricingTime;
	double totalTime;	
	const char* filename;

	
	void writeData();



};




