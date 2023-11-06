#include <cmath>
#include <math.h>
#include <random>
#pragma once
using namespace std;

class c_r //classe de variable al�atoire suivant une loi normale centr�e r�duite
{
    public:
    std::default_random_engine gen;
    std::uniform_real_distribution<double> dist;
    double simuler(); //fonction servant � faire un tirage suivant la loi N(0,1)
    c_r();
    ~c_r();
};
class W
{
public:
    double rho;
    double dt; //pas de temps sur lequel on simule le brownien
    double rac_dt; //racine carr�e de dt, on fait le choix ici de la prendre en constante au lieu de la prendre en argument de fonction pour �viter de calculer la racine carr� � chaque it�ration
    c_r loi_normale;
    double lambda1; //param�tre entrant en compte dans la factorisation de Choleski
    double lambda2; //param�tre entrant en compte dans la factorisation de Choleski
    double valeurs[2]; //vecteur W simul�
    void constantes(); //calcule les constantes utiles (rac_dt, lambda1, lambda2)
    void iterer(); //effectue la simulation de W
};

double normale(double moyenne, double variance, c_r& normale_c_r);
double N(double x);
double Put_price(double K,double S0,double sigma,double r, double T);
void ecriture(string nom_fichier, int n_points, int* X, double* Y, double* Y2=nullptr, double* Y3=nullptr, double* Y4=nullptr);
void ecriture(string nom_fichier, int n_points, double* X, double* Y, double* Y2=nullptr, double* Y3=nullptr, double* Y4=nullptr);
void ecriture_couple(string nom_fichier, int n_points, double* X, double** Y);
void Monte_Carlo_classique_echange(double r, double T, double sigma1, double sigma2, double rho, int n_tirages, double alphaS10, double betaS20, c_r& normale_c_r, double coefficient, double* variables);
void Monte_Carlo_condi_echange(double r, double T, double sigma1, double sigma2, double rho, int n_tirages, double alphaS10, double betaS20, c_r& normale_c_r, double coefficient, double* variables);
void Monte_Carlo_classique_spread(double r, double T, double sigma1, double sigma2, double rho, int n_tirages, double alphaS10, double betaS20, c_r& normale_c_r, double coefficient, double* variables, double K);
void Monte_Carlo_condi_spread(double r, double T, double sigma1, double sigma2, double rho, int n_tirages, double alphaS10, double betaS20, c_r& normale_c_r, double coefficient, double* variables, double K);
void Monte_Carlo_controle_spread(double r, double T, double sigma1, double sigma2, double rho, int n_tirages, double alphaS10, double betaS20, c_r& normale_c_r, double coefficient, double* variables, double K);
void Monte_Carlo_controle_spread2(double r, double T, double sigma1, double sigma2, double rho, int n_tirages, double alphaS10, double betaS20, c_r& normale_c_r, double coefficient, double* variables, double K);
void Monte_Carlo_Best_Of(double r, double T, double sigma1, double sigma2, double rho, int n_tirages, double alphaS10, double betaS20, c_r& normale_c_r, double coefficient, double* variables, double K);
void logspace(int debut, int fin, int n, int *valeurs);
void logspace(int debut, int fin, int n, int *valeurs);
void logspace(double debut, double fin, int n, double *valeurs);
void linspace(int debut, int fin, int n, int *valeurs);
void linspace(double debut, double fin, int n, double *valeurs);
