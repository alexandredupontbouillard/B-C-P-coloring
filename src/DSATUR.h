#ifndef DSATUR_H
#define DSATUR_H
//////////////////////////////////////
// INCLUDE & DEFINE
//////////////////////////////////////
#include "Graphe.h"
#include <stdio.h>
#include <stdlib.h>
#include <string> 
#include <fstream>
#include <list>
#include <iostream>
#define ORDER_MAX 10000

//////////////////////////////////////
// STRUCTURE
//////////////////////////////////////
typedef struct t_graphe Type_G;
struct t_graphe { 
 int order;
 int matrice[ORDER_MAX][ORDER_MAX];
};

class DSATUR{
	public :
	
	vector<int> degre;
	vector<int> dsat;
	vector<int> coloration_issu_de_dsat;

	DSATUR(int n){
		degre.resize(n);
		dsat.resize(n);
		coloration_issu_de_dsat.resize(n);
	}

	void coloration(Type_G g,vector< list<int> > &coloration);
	Type_G init_bis(Type_G g,C_Graphe &g_);
	Type_G init(Type_G g);
	void init_dsat(Type_G g);
	void parcours_profondeur(Type_G g, unsigned int som_depart);
	void parcours_profondeur_rec(Type_G g, unsigned int s, int *marque);
	void parcours_largeur(Type_G g, unsigned int s);
	int * copyYdansZ(int * destination, int *origine, Type_G g);
	int premier_sommet_de_Z (int * vecteur,  Type_G g);
	int * retirer_s_et_ces_voisins (int s, int * vecteur,  Type_G g);
	int * glouton(Type_G g);
	int sommet_non_colorie_de_dsat_max(Type_G g);
	int coloration_avec_couleur_mininum(int t, Type_G g);
	int nb_de_couleur_differente_autour_de (int t, Type_G g);
	void maj_du_vecteur_dsat (int t, Type_G g);
	void color_by_dsatur(Type_G g);
	int verif_color(Type_G g);
	void display_graph(Type_G g);
	void display_vecteur(int *vecteur, Type_G g);
	int nb_chromatique(Type_G g);
	void lecture_dimacs_aretes();
};
#endif
