#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <chrono>
#include "fcts.h"
#include "methode.h"
#include "util.h"


// Paramètres de la simulation
const double tf = 1e1;                    // Temps final
const int nbSimMin = 1e4;
const int nbSimMax = 1e5;
const int nbFic = nbSimMax / 1e3;         // # de pts des fichiers (Pos,En), doir être supérieur à nbSim
const int nbNu = 2e0;                     // # de pts du fichier (Nu)
const double L = 11.;
const double rc = 2.5;
const double rc2 = 1. / (rc * rc);

int main(int argc, char *argv[]) {
	// Nom des fichiers
	const std::string nomSim = argv[1];                                  // Nom de la simulation
	const std::string nomCI = "./condition_initiale/"+nomSim+".csv";     // Nom du fichier à étudier
	const std::string nomPosX = "./resultats/"+nomSim+"_position_x.csv"; // Nom du fichier des positions suivant X
	const std::string nomPosY = "./resultats/"+nomSim+"_position_y.csv"; // Nom du fichier des positions suivant Y
	const std::string nomEn = "./resultats/"+nomSim+"_energie.csv";      // Nom du fichier des énergies
	const std::string nomNu = "./resultats/"+nomSim+"_nu.csv";           // Nom du fichier pour calculer nu

	// Calcul du # de particule(s)
	const int N = compteLigneFichier(nomCI);
	std::string methode;

	// Allocation mémoire des matrices position, vitesse et force (2D)
	// Colonne = particule, ligne = instant donné
	std::vector< std::vector<double> > X(nbSimMax, std::vector<double>(N,0));
	std::vector< std::vector<double> > Xp(nbSimMax, std::vector<double>(N,0));
	std::vector< std::vector<double> > Xpp(nbSimMax, std::vector<double>(N,0));
	std::vector< std::vector<double> > Y(nbSimMax, std::vector<double>(N,0));
	std::vector< std::vector<double> > Yp(nbSimMax, std::vector<double>(N,0));
	std::vector< std::vector<double> > Ypp(nbSimMax, std::vector<double>(N,0));

	// Implémentation des Conditions Initiales
	X[0] = extraireColonne(nomCI,0);
	Xp[0] = extraireColonne(nomCI,2);
	Y[0] = extraireColonne(nomCI,1);
	Yp[0] = extraireColonne(nomCI,3);

	// Actualisation des forces à l'instant 0
	for(int i = 0; i < N; i++) {
		Xpp[0][i] = FLJX(i,X[0],Y[0]);
		Ypp[0][i] = FLJY(i,X[0],Y[0]);
	}

	// Allocation mémoire des quantités qui permettent de calculer nu et l'énergie
	std::vector<double> deltaE(nbNu, -1. * EmLJ(X[0],Y[0],Xp[0],Yp[0])), dt(nbNu), nbSim(nbNu);
	std::vector<double> Em(nbFic), t(nbFic);

	// Boucle pour calculer nbSim
	if(nbNu <= 1) {
		nbSim[0] = nbSimMax;
	}
	else {
		double powMin = log10(nbSimMin);
		double powMax = log10(nbSimMax);
		for(int i = 0; i < nbNu; i++) {
			nbSim[i] = pow(10., powMin + (powMax - powMin) * i / (nbNu - 1));
		}
	}

	// Crétaion d'une variable pour calculer la durée de la boucle
	std::chrono::high_resolution_clock::time_point debut = std::chrono::high_resolution_clock::now();

	// Boucle principale de clacul
	#ifdef EULER
	methode = "Euler";
	std::cout << " > Méthode de Euler séléctionnée" << std::endl;
	// Boucle principale -> calcul des positions, vitesses et forces
	for(int k = 0; k < nbNu; k++) {
		dt[k] = tf / nbSim[k];

		for(int i = 0; i < nbSim[k]-1; i++) {
			for(int j = 0; j < N; j++) {
				// Actualisation des forces
				Xpp[i][j] = FLJX(j,X[i],Y[i]);
				Ypp[i][j] = FLJY(j,X[i],Y[i]);

				// Actualisation des positions
				X[i+1][j] = eulerQ(X[i][j],Xp[i][j],dt[k]);
				Y[i+1][j] = eulerQ(Y[i][j],Yp[i][j],dt[k]);

				// Replacement des positions dans la boite (Entre -L/2 et L/2)
				boite(X[i+1][j],L);
				boite(Y[i+1][j],L);

				// Actualisation des vitesses
				Xp[i+1][j] = eulerQp(Xp[i][j],Xpp[i][j],dt[k]);
				Yp[i+1][j] = eulerQp(Yp[i][j],Ypp[i][j],dt[k]);
			}
		}
		deltaE[k] += EmLJ(X[nbSim[k]-1],Y[nbSim[k]-1],Xp[nbSim[k]-1],Yp[nbSim[k]-1]);
	}

	#elif VV

	methode = "Velocity-Verlet";
	std::cout << " > Méthode de Velocity-Verlet séléctionnée" << std::endl;
	// Velocity-Verlet -> calcul des positions, vitesses et forces
	for(int k = 0; k < nbNu; k++) {
		dt[k] = tf / nbSim[k];

		for(int i = 0; i < nbSim[k]-1; i++) {
			for(int j = 0; j < N; j++) {
				// Actualisation des positions
				X[i+1][j] = velocityVerletQ(X[i][j],Xp[i][j],Xpp[i][j],dt[k]);
				Y[i+1][j] = velocityVerletQ(Y[i][j],Yp[i][j],Ypp[i][j],dt[k]);

				// Replacement des positions dans la boite (Entre -L/2 et L/2)
				boite(X[i+1][j],L);
				boite(Y[i+1][j],L);
			}
			for(int j = 0; j < N; j++) {
				// Actualisation des forces
				Xpp[i+1][j] = FLJX(j,X[i+1],Y[i+1]);
				Ypp[i+1][j] = FLJY(j,X[i+1],Y[i+1]);

				// Actualisation des vitesses
				Xp[i+1][j] = velocityVerletQp(Xp[i][j],Xpp[i][j],Xpp[i+1][j],dt[k]);
				Yp[i+1][j] = velocityVerletQp(Yp[i][j],Ypp[i][j],Ypp[i+1][j],dt[k]);
			}
		}
		deltaE[k] += EmLJ(X[nbSim[k]-1],Y[nbSim[k]-1],Xp[nbSim[k]-1],Yp[nbSim[k]-1]);
	}
	#endif

	// Calcul du temps mis par la boucle for
	std::chrono::high_resolution_clock::time_point fin = std::chrono::high_resolution_clock::now();
	double duree = std::chrono::duration_cast<std::chrono::duration<double> >(fin - debut).count();
	std::cout << " > durée du calcul : " << std::fixed << std::setprecision(2) << duree << " s" << std::endl;

	// Boucle pour créer les données sur l'énergie
	int valeurPasEnergie = nbSim[nbNu-1] / nbFic;
	for(int i = 0; i < nbFic; i++) {
		t[i] = i * valeurPasEnergie * dt[nbNu-1];
		Em[i] = EmLJ(X[i*valeurPasEnergie],Y[i*valeurPasEnergie],
		Xp[i*valeurPasEnergie],Yp[i*valeurPasEnergie]);
	}

	// Sauvegarder les resultats
	ecrirePos(nomPosX,X,nbFic,dt[nbNu-1],N,methode);
	ecrirePos(nomPosY,Y,nbFic,dt[nbNu-1],N,methode);
	ecrireEnNu(nomEn,t,Em,dt[nbNu-1],N,methode);
	ecrireEnNu(nomNu,dt,deltaE,dt[nbNu-1],N,methode);

	// Tests
	// std::vector<double> x(2), y(2), abs(100), ordo(100);
	// x[0] = 0.; y[0] = 0.;
	// x[1] = 2.5; y[1] = 0.;
	// std::cout << " > Force test : " << FLJX(0,x,y) << std::endl;

	return 0;
}