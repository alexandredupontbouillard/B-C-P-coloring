#include "Graph.h"



#define epsilon 0.00001


/****************************  StableSet *****************************/

bool StableSet::contains(int n){

	return nodes[n]==1;

}

void StableSet::add(int n){
	
	if(nodes[n]==0){

		nodes[n]=1;
		L_nodes.push_back(n);
	}
	


}

void StableSet::print(){

	list<int>::iterator it;

	for(it = L_nodes.begin(); it != L_nodes.end(); it ++){
		cout<<(*it)<<" ";


	}
	cout<<endl;

}
	
void StableSet::initialize(int nb_nodes){
	L_nodes.clear();
	nodes = vector<int>(nb_nodes);

}

int StableSet::size(){

	return L_nodes.size();
}


/****************************  C_link  *******************************/

int C_link::return_other_extrem(int v){
	return (v==v1?v2:v1);
}






/***************************  C_node  *****************************/

bool C_node::test_neighbour(int j){
  list<C_link*>::iterator it;
  for(it=L_adjLinks.begin() ; it !=L_adjLinks.end() ; it++){
    if((*it)->return_other_extrem(num) == j)
      return true;
  }
  return false;
}

bool C_node::test_successor(int j){
  list<C_link*>::iterator it;
  for(it=L_adjLinks.begin() ; it !=L_adjLinks.end() ; it++){
    if((*it)->return_other_extrem(num) == j)
      return true;
  }
  return false;
}

/**************************  C_Graph  ******************************/




void C_Graph::completeIntoMaxStable(StableSet* stable){  // fonction à appeler pour compléter le stable

	vector<int> allV = vector<int>(nb_nodes);

	list<int>::iterator it;
	for (it = stable->L_nodes.begin(); it != stable->L_nodes.end(); ++it) {
	    allV[*it] = 1;
	    this->deleteNeighbourhood(*it, &allV);   // on supprime de notre liste de sommet allV, l'ensemble des voisins des sommets de stable
            
	    	

	}
	this->complete(stable,&allV);  // on complète le stable
		



}

void C_Graph::deleteNeighbourhood(int sommet, vector<int>* allV){  // fait passer à 1 la valeur de allV indicée par les sommets adjacents à sommet
	list<C_link*>::iterator it;
	
	for(it =V_nodes[sommet].L_adjLinks.begin() ; it != V_nodes[sommet].L_adjLinks.end() ; ++it){

		(*allV)[(*it)->return_other_extrem(sommet)]=1;

	}


}

void C_Graph::complete(StableSet* stable, vector<int>* allV){

	
	for(long unsigned int i = 0; i < allV->size(); i ++){ // dès lors qu'un sommet adjacent à aucun sommet déjà ajouté dans le stable est trouvé, on l'ajoute puis on supprime ses voisins des sommets restants
		if((*allV)[i]==0){
			stable->add(i);
			this->deleteNeighbourhood(i,allV);

		}


	}

}

void C_Graph::getGloutColoring(list<StableSet*>& coloration){
	
	StableSet* stablemax;

	vector<int> sommets = vector<int>(nb_nodes);


	list<int>::iterator it;

	for(int i = 0 ; i < nb_nodes; i ++){
		stablemax= new StableSet;
		stablemax->initialize(nb_nodes);
		





		for(int j = 0 ; j < nb_nodes ; j ++){   // on cherche le premier sommet non couvert de la liste de sommets
			if(sommets[j]==0){
				stablemax->add(j);
				j = nb_nodes;

			}


		}
		if(stablemax->size()>0){
			completeIntoMaxStable(stablemax);   // on complète le wingleton en stable max
			// on met à jour la liste de sommets 
			for(it = stablemax->L_nodes.begin() ; it != stablemax->L_nodes.end() ; it ++){
				sommets[*it] = 1;


			}
			coloration.push_back(stablemax);
		}
		else{
			i = nb_nodes;
		}
		
		
		


	}




}

C_Graph  C_Graph::getComplementary(){
	int i ,j , nb;
	C_Graph g;
	C_link* a;


	g.nb_nodes = nb_nodes;
	g.nb_links = nb_nodes*(nb_nodes-1) /2 - nb_links;
	g.adjacency.resize(nb_nodes);
	for( i = 0 ; i < nb_nodes; i++){

		g.adjacency[i]= vector<int>(nb_nodes,1);

    	}
	
	g.V_nodes.resize(nb_nodes);
        g.V_links.resize(g.nb_links);

        for (i=0;i<nb_nodes;i++){
        	g.V_nodes[i].num = i;
        	g.V_nodes[i].L_adjLinks.clear();
        	g.V_nodes[i].weight=0;
		g.adjacency[i][i] = 0;
        }
	

	for(i = 0 ; i < (int)V_links.size() ; i++){
		g.adjacency[V_links[i]->v1][V_links[i]->v2] = 0;
		g.adjacency[V_links[i]->v2][V_links[i]->v1] = 0;

	}
	
	nb = 0;
	
	for(i = 0; i < nb_nodes ; i ++ ) {

		for(j = 0 ; j <i ; j ++){
			if(g.adjacency[i][j]==1){
				a = new C_link;
				a->v1 = min(i,j);
				a->v2 = max(i,j);
				a->weight=0;
				a->num=nb;
				g.V_links[nb] = a;
				nb++;
			}

		}

	}
	
	return g;


    
	
}



// A FINIR



void C_Graph::findViolatedCycleIq(vector<StableSet*>* stables, vector<float>& cost, list<list<int>*> &result, int rhs=3){ // détecte de manière gloutonne si une inégalité de cycle a été détectée
	


	int pred1,pred2;
	list<int>* cycle;
	pair<float,list<int>> res;
	vector<int>  ss;
	vector<int> sommet_interdit;
	int v1,v2;
	float violation;
	pair<float,int> newSommet;
	bool oncontinue;
	
	for(int i = 0 ; i < nb_links ; i ++){

		
		cycle=new list<int>();
		sommet_interdit = vector<int>(nb_nodes,0);
		violation = 0;	
		v1 = V_links[i]->v1;
		v2 = V_links[i]->v2;
		
		cycle->push_back(v1);						//on ajoute au cycle le premier sommet de l'arête
		ss = vector<int> (stables->size(), 0);
		
		cycle->push_back(v2);						//on ajoute au cycle le second sommet de l'arête
		pred1 = v2;
		pred2 = v1;
		
		computeCostOfVertice( stables ,   &ss, cost, v1,&res);	//on calcul le coût de chaque sommet et on met à jour l'ensemble des stables
					
		majSS(&res.second, &ss);
		
		violation = violation + res.first;

	        computeCostOfVertice( stables ,   &ss, cost, v2, &res );
	       
		majSS(&res.second, &ss);

		
		violation = violation + res.first;
		sommet_interdit[v1]=1;
		sommet_interdit[v2]=1;
		newSommet= ajouteSommetAuCycle( cycle, stables ,  &ss, cost ,  v1, v2, rhs-violation,&sommet_interdit);
		oncontinue = false;
		violation = violation + newSommet.first; 

		if(newSommet.first !=-1){		
			sommet_interdit[newSommet.second]=1;
			oncontinue = adjacency[v1][newSommet.second] == 0 || adjacency[v2][newSommet.second] == 0 ;
			if(cycle->front()  == v1){
				pred2 = v2;
				v2 = cycle->back();
				

			}
			else{
				pred1 = v1;
				v1 = cycle->front();
			}
		}
		
	
		while(oncontinue){
			
			newSommet= ajouteSommetAuCycle( cycle, stables ,  &ss, cost ,  v1, v2, rhs-violation,&sommet_interdit);// on ajoute un sommet dans le cycle
			
			if(newSommet.second != -1 ){			
				sommet_interdit[newSommet.second]=1;									// on met à jour la liste de sommets interdits
				violation = violation + newSommet.first;
											// on met à jour la violation du cycle partiel
				if(cycle->front() ==newSommet.second ){									// on met à jour les extrémités du cycle
					pred1=v1;
					v1 = newSommet.second;
				}	
				else{
					pred2=v2;
					v2= newSommet.second;
				}
				
				
				
				newSommet= ajouteSommetAuCycle( cycle, stables ,  &ss, cost ,  v1, v2, rhs-violation,&sommet_interdit);// on ajoute un sommet dans le cycle
				if(newSommet.second!=-1){					
					sommet_interdit[newSommet.second]=1;									// on met à jour la liste de sommets interdits
					violation = violation + newSommet.first;								// on met à jour la violation du cycle partiel
					if(cycle->front() ==newSommet.second ){									// on met à jour les extrémités du cycle
						pred1=v1;
						v1 = newSommet.second;
					}	
					else{
						pred2=v2;
						v2= newSommet.second;
						
					}
			
					if(!  verifyChord(cycle,pred1,pred2,v1,v2)  ){
						
						if(adjacency[v1][v2] == 1){
							
							if(!cycleIn(result,cycle)){

								result.push_back(cycle);
							}
                                                        else
                                                                free(cycle);
							oncontinue = false;
						}
					}else{
						oncontinue=false;
						free(cycle);

					}

				}
				else{
				oncontinue = false;
				free(cycle);
				}
			
				
	 		}
			else{
				oncontinue  = false;
				free(cycle);

			}
			
			


		}
	}


	
	

	



}

bool C_Graph::verifyChord(list<int>* cycle,int pred1,int pred2,int v1, int v2){ //renvoit true si le cycle partiel contient une corde passant par ses extrémités
	//cout<<pred1<<" "<<pred2 << " "<<v1 << " "<<v2<<endl;

	list<int>::iterator itt;
	long unsigned int c= 1;
	for(itt = next(cycle->begin(),1); c<cycle->size(); itt++){
		if(adjacency[cycle->front()][(*itt)] == 1 && pred1!= (*itt) && v2 != (*itt) ){
		
			return true;

		}
		c++;

		
	}
	c =1 ;
	list<int>::reverse_iterator it;
	for(it = next(cycle->rbegin(),1); c < cycle->size(); it++){
		c++;
		if(adjacency[cycle->back()][(*it)] == 1 && pred2!= (*itt) && (*itt) != v1 ){
			
			return true;

		}

		
	}
	return false;


}









pair<float,int> C_Graph::ajouteSommetAuCycle(list<int>* cycle,vector<StableSet*>* stables , vector<int>* ss,vector<float>& cost ,  int v1, int v2, float placedispo,vector<int>* sommet_interdit){ //retourne le coût du sommet ajouté
	 list <C_link*>::iterator it;
	 
	list<int> mmm; 
	
	pair<float,list<int>> min = make_pair(100,mmm);

	
	pair<float,list<int>> m;
	int verticeMin2,k;
	int verticeMin1 = 0;
	 for( it= V_nodes[v1].L_adjLinks.begin() ; it != V_nodes[v1].L_adjLinks.end(); it++ ){
		
		k = (*it)->return_other_extrem(v1);

		if((*sommet_interdit)[k] != 1){
			
			 computeCostOfVertice(stables,ss,cost,k,&m);
			if(m.first < min.first){
				verticeMin1 = k;
				min.first = m.first;
				min.second = m.second;

			}
		}

	}

	verticeMin2 = -1;
	for( it= V_nodes[v2].L_adjLinks.begin() ; it != V_nodes[v2].L_adjLinks.end(); it++ ){
		
		k = (*it)->return_other_extrem(v2);

		if((*sommet_interdit)[k]!=1){		
			 computeCostOfVertice(stables,ss,cost,k,&m);
			
			if(m.first < min.first){
				verticeMin2 = k;
				min.first = m.first;
				min.second = m.second;

			}
		}

	}

	if(min.first < placedispo-epsilon){  // le sommet peut être ajouté

		if(verticeMin2 == -1){ // le sommet à ajouter est stocké dans verticeMin1, le sommet doit donc être ajouté à gauche dans la liste

			cycle->push_front(verticeMin1);
			
			verticeMin2= verticeMin1;
			
			
			

		}
		else{ // le sommet à ajouter est stocké dans verticeMin2, le sommet doit donc être ajouté à droite dans la liste.
			
			cycle->push_back(verticeMin2);
			
			

		}


	
		majSS(&min.second, ss);

		return make_pair(min.first,verticeMin2);


	}else{ // le sommet ne peut pas être ajouté, donc on renvoit stop
		return make_pair(-1,-1);

	}




	


}


// FINIT

bool cycleIn(list<list<int>*>& l1, list<int>* l2){  // return true si la liste l1 contient liste l2 

	list<list<int>*>::iterator it ;
	for(it = l1.begin() ; it != l1.end() ; it ++){
		if(compareCycle((*it),l2)) return true;

	}
	return false;
}




bool compareCycle(list<int>* l1, list<int>* l2){  // return true si les deux listes contiennent les mêmes sommets
	list<int>::iterator it1;
	list<int>::iterator it2;
	bool b;
	for(it1 = l1->begin() ; it1 != l1->end() ; it1++){
		b = false;
		for(it2 = l2->begin(); it2 != l2->end() ; it2++){
			if( (*it1) == (*it2) ) b = true;


		}

		if(!b) return false;
	

	}


	return true;
}

int in(list<int> cycle, int vertice){

	int position = 0;
	list<int>::iterator it;

	for(it = cycle.begin() ; it != cycle.end() ; it++){
		if((*it) == vertice) {

			return position;
		}
		position=position+1;
	}	
	return -1;


}

void majSS(list<int>* indices, vector<int>* ss){

	list<int>::iterator it;
	for(it= indices->begin(); it != indices->end(); it ++){
		(*ss)[(*it)]=1;

	}

}




void computeCostOfVertice(vector<StableSet*>* stables , vector<int>*  ss,vector<float>& cost, int v1, pair<float,list<int>>*  m){
	float sum = 0;
	list<int> ll;
	m->second.clear();
	
	for(long unsigned int i = 0 ; i < stables->size(); i ++){
		
		if((*ss)[i]==0 && (*stables)[i]->contains(v1)){
			sum = sum + cost[i];
			
			m->second.push_back(i);
		}



	}
		
	m->first = sum;
	


}




void C_Graph::read_undirected_DIMACS(istream & fic){
  if (!fic){
    cout<<"File Error"<<endl;
  }else{
    int k,i,j;
    string m1,m2;
    list<C_link>::iterator it;
    C_link *a;

    fic>>m1;
    fic>>m2;
		
    // Jump the file description and go to the data
    while (((m1!="p")&&(m2!="edge"))||((m1!="p")&&(m2!="col"))){
      m1=m2;
      fic>>m2;
    }

    directed=false;
    
    fic>>nb_nodes;
    fic>>nb_links;
		
    adjacency = vector<vector<int>>(nb_nodes);
    for( i = 0 ; i < nb_nodes; i++){

		adjacency[i]= vector<int>(nb_nodes);

    }
	

    V_nodes.resize(nb_nodes);
    V_links.resize(nb_links);

    for (i=0;i<nb_nodes;i++){
      V_nodes[i].num = i;
      V_nodes[i].L_adjLinks.clear();
      V_nodes[i].weight=1;
    }

    for (k=0;k<nb_links;k++){
      fic>>m1;
      fic>>i;
      fic>>j;
      adjacency[i-1][j-1]=1;
      adjacency[j-1][i-1]=1;
      a=new C_link;
      a->num=k;
      a->v1=min(i-1,j-1);
      a->v2=max(i-1,j-1);
      a->weight=0;
      V_nodes[i-1].L_adjLinks.push_back(a);
      V_nodes[j-1].L_adjLinks.push_back(a);
      V_links[k] = a;
    }
	
  }

}




void C_Graph::write_dot_G(string InstanceName){
  list<C_link>::iterator it;
  int i,k;

  ostringstream FileName; 
  FileName.str("");
  FileName <<InstanceName.c_str() << "_G.dot";

  ofstream fic(FileName.str().c_str());

  if (!directed) {
  
  fic<<"graph G {"<<endl;

  for(i=0 ; i<nb_nodes ; i++)
      fic<<"  "<<V_nodes[i].num<<"[shape = octagon]"<<endl;

  for(k=0 ; k<nb_links ; k++)
      fic<<"  \""<<V_links[k]->v1<<"\"--\""<<V_links[k]->v2<<"\";"<<endl;
  
  fic<<"}"<<endl;

  }
  else{


  }

  

  fic.close();
  ostringstream commande; 
  commande.str("");
  commande<<"dot -Tpdf -o "<<InstanceName.c_str() << "_G.pdf "<< FileName.str().c_str()<<endl;
  cout<<commande.str().c_str();
  if(system(commande.str().c_str())){cout<<"PDF generated successfully"<<endl;}
  return;
}



void C_Graph::write_dot_G_stableset(string InstanceName, vector<int>& stable){
  int i,k;
  ostringstream FileName; 
  FileName.str("");
  FileName <<InstanceName.c_str() << "_G_stable.dot";

  ofstream fic(FileName.str().c_str());
  
  fic<<"graph G {"<<endl;
  
  for(i=0 ; i<nb_nodes ; i++){
    if (stable[i]>1-epsilon)
      fic<<"  "<<V_nodes[i].num<<"[shape = octagon]"<<endl;
    else
      if (stable[i]<epsilon)
	fic<<"  "<<V_nodes[i].num<<"[shape = octagon, color=white]"<<endl;
  }
  
  for(k=0 ; k<nb_links ; k++){
    if ((stable[V_links[k]->v1]>1-epsilon)&&(stable[V_links[k]->v2]>1-epsilon))
      fic<<"  \""<<V_links[k]->v1<<"\"--\""<<V_links[k]->v2<<"\";"<<endl;
    else
      fic<<"  \""<<V_links[k]->v1<<"\"--\""<<V_links[k]->v2<<"\" [color=white];"<<endl;
  }
  
  fic<<"}"<<endl;
		
  fic.close();

  ostringstream commande; 
  commande.str("");
  commande<<"dot -Tpdf -o "<<InstanceName.c_str() << "_G_stable.pdf "<< FileName.str().c_str()<<endl;
  cout<<commande.str().c_str();
  if (system(commande.str().c_str())){cout<<"PDF generated successfully"<<endl;}


}





void C_Graph::write_dot_G_color(string InstanceName, vector<int>& coloring){
  int i,k;
  ostringstream FileName; 
  FileName.str("");
  FileName <<InstanceName.c_str() << "_G_color.dot";

  vector <string> colors;
  colors.push_back("green");
  colors.push_back("blue");
  colors.push_back("red");
  colors.push_back("cyan");
  colors.push_back("yellow");
  colors.push_back("magenta");
  colors.push_back("darkorchid");
  colors.push_back("darkorange");
  colors.push_back("deeppink");
  colors.push_back("forestgreen3");
  colors.push_back("indigo");
  colors.push_back("midnightblue");
  colors.push_back("violetred");

  int chi=0;
  for (i=0;i<nb_nodes;i++)
    if (chi<coloring[i]) chi=coloring[i];

  
  if(chi >= (int)colors.size()){
    cout<<"We only have 13 colors and this solutions needs "<<chi<<" colors... some nodes will have wrong colors!"<<endl;
  }
  
  ofstream fic(FileName.str().c_str());
  fic<<"graph G {"<<endl;
  
  for(i=0 ; i<nb_nodes ; i++){
    fic<<"  "<<V_nodes[i].num<<"[shape = octagon, style = filled , fillcolor = "<<colors[(coloring[V_nodes[i].num]) % colors.size()]<<" ]"<<endl;
  }
  
  for(k=0 ; k<nb_links ; k++){
      fic<<"  \""<<V_links[k]->v1<<"\"--\""<<V_links[k]->v2<<"\";"<<endl;
  }
  
  fic<<"}"<<endl;
		
  fic.close();

  ostringstream commande; 
  commande.str("");
  commande<<"dot -Tpdf -o "<<InstanceName.c_str() << "_G_color.pdf "<< FileName.str().c_str()<<endl;
  cout<<commande.str().c_str();
  if (system(commande.str().c_str())){cout<<"PDF generated successfully"<<endl;}
  

}
