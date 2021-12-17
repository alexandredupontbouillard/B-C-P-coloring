#ifndef Graph_H
#define Graph_H

#include <vector>
#include <list>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>

#include <math.h>
#include <unistd.h>
#include <utility>
using namespace std;


/****************************  Stable Set ****************************/

class StableSet{
public :
	list<int> L_nodes;
	vector<int> nodes;


	bool contains(int n);
	void add(int n);
	
	void initialize(int nb_nodes);

	int size();

	void print();



};


/****************************  C_edge  *******************************/
class C_link{
public:
  int num;      // Number of the edge
  int v1, v2;   // The two extremities of an edge v1v2 or of an arc (v1,v2)
  float weight;
  
  // return the extremity disctinc from v in O(1).
  int return_other_extrem(int v);
  
};


/***************************  C_node  *****************************/
class C_node{
public :
   int num;     // Number of the node
   float weight;
   
   list <C_link*> L_adjLinks;

   //Test if j is a neighbour of i in O(degre(i))
   bool test_neighbour(int j);

   //Test if j is a successor of i in O(degre(i))
   bool test_successor(int j);

};


/**************************  C_Graph  ******************************/
class C_Graph{
public:

  bool directed;  // True if directed / False if undirected
  int nb_nodes;   // Number of nodes
  int nb_links;   // Number of links
  vector<vector<int>> adjacency;  //

  // Encoding of the graph by adjacence list, i.e. a vector of list of edges 
  vector <C_node> V_nodes;

  // Additional encoding: a vector on the edges (on pointers over edges)
  vector <C_link*> V_links;

  /*********************************************/
  /*********** ALGORITHMS **********************/

  bool detect_cycle(vector<int>&sol);
  bool verifyChord(list<int> cycle);
  pair<float,int> ajouteSommetAuCycle(list<int>* cycle,vector<StableSet*>* stables , vector<int>* ss,vector<float>& cost ,  int v1, int v2, float placedispo,vector<int>* sommet_interdit);
  void findViolatedCycleIq(vector<StableSet*>* stables, vector<float>& cost,list<list<int>*> &result, int rhs);
  // to get max stable set from stable set (greedy)
  void completeIntoMaxStable(StableSet* stable);
  void complete(StableSet* stable, vector<int>* allV);
  void deleteNeighbourhood(int sommet, vector<int>* allV);
  bool verifyChord(list<int>* cycle,int v1,int v2,int pred1, int pred2);
  void getGloutColoring(list<StableSet*>& coloration);


  C_Graph  getComplementary();

  //list<list<int>>* findViolatedCycleIq(vector<StableSet*> stables, vector<float> cost);   

  /*********************************************/
  /*********** INPUT-OUTPUT FILES **************/
  
  // Read a DIMACS file and store the corresponding graph in C_Graph
  void read_undirected_DIMACS(istream & fic);

  // Read a TSP file and store the corresponding graph in C_Graph
  void read_undirected_TSP(istream & fic);
  
  // Read a directed "gra" format file and store the corresponding graph in C_Graph
  void read_directed_GRA(istream & in);

  // Write a Graphviz File with the DOT format
  void write_dot_G(string InstanceName);
  
  // Write a Graphviz File with the DOT format using an incidence vector of a stable set
  void write_dot_G_stableset(string InstanceName, vector<int> &stable);

  // Write a Graphviz File with the DOT format using a coloration vector
  void write_dot_G_color(string InstanceName, vector<int> &coloring);
  
};

void computeCostOfVertice(vector<StableSet*>* stables , vector<int>*  ss,vector<float>& cost, int v1,pair<float,list<int>>* m);
void majSS(list<int>* indices, vector<int>* ss);
bool cycleIn(list<list<int>*>& l1, list<int>* l2);
bool compareCycle(list<int>* l1, list<int>* l2);
int in(list<int> cycle, int vertice);
#endif
