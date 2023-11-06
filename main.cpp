#include <iostream>
#include <ctime>
#include "Monte_Carlo.hpp"

using namespace std;
//ATTENTION : dans tous les tirages on ne prend pas en compte l'actualisation pour éviter une multiplication par tirage mais on la prend en compte dans les résultats (variances, moyennes...)
int main()
{
    //initialisations
    std::random_device rd;
    std::default_random_engine generator(time(0));
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    c_r normale_c_r;
    normale_c_r.gen = generator;
    normale_c_r.dist = distribution;
    double variables[4]; //contient prix, variance, intervalle bas, intervalle haut
    int n_tirages_max = 200000;
    int n_trajectoires[100];
    double liste_variances[100];
    double liste_prix[100];
    double liste_i_down[100];
    double liste_i_up[100];
    logspace(100, n_tirages_max, 100, n_trajectoires);

    //Option d'échange : Monte-Carlo classique
    for (int i=0; i<100; i++)
    {
        Monte_Carlo_classique_echange(0.01, 2.0, 0.25, 0.3, 0.5, n_trajectoires[i], 1.0, 1.0, normale_c_r, 1.645, variables);
        liste_prix[i] = variables[0];
        liste_variances[i] = variables[1];
        liste_i_down[i] = variables[2];
        liste_i_up[i] = variables[3];
    }
    ecriture("prix_variance_intervalles_90", 100, n_trajectoires, liste_prix, liste_variances, liste_i_down, liste_i_up);

    //Option d'échange : Monte-Carlo avec variable de conditionnement
    for (int i=0; i<100; i++)
    {
        Monte_Carlo_condi_echange(0.01, 2.0, 0.25, 0.3, 0.5, n_trajectoires[i], 1.0, 1.0, normale_c_r, 1.645, variables);
        liste_prix[i] = variables[0];
        liste_variances[i] = variables[1];
        liste_i_down[i] = variables[2];
        liste_i_up[i] = variables[3];

    }
    ecriture("prix_variance_intervalles_90_2", 100, n_trajectoires, liste_prix, liste_variances, liste_i_down, liste_i_up);

    //Option d'échange : Monte-Carlo avec réduction de variance par conditionnement
    double liste_rho[100];
    linspace(-0.99, 0.99, 100, liste_rho);
    for (int i=0; i<100; i++)
    {
        Monte_Carlo_condi_echange(0.01, 2.0, 0.25, 0.3, liste_rho[i], 27000, 1.0, 1.0, normale_c_r, 1.645, variables);
        liste_prix[i] = variables[0];
    }
    ecriture("prix_fonction_rho", 100, liste_rho, liste_prix);



    //Option sur Spread : Monte-Carlo classique
    for (int i=0; i<100; i++)
    {
        Monte_Carlo_classique_spread(0.01, 2.0, 0.25, 0.3, 0.5, n_trajectoires[i], 1.2, 1.0, normale_c_r, 1.645, variables, 0.2);
        liste_prix[i] = variables[0];
        liste_variances[i] = variables[1];
        liste_i_down[i] = variables[2];
        liste_i_up[i] = variables[3];
    }
    ecriture("prix_variance_intervalles_90_spread", 100, n_trajectoires, liste_prix, liste_variances, liste_i_down, liste_i_up);

    //Option sur Spread : Monte-Carlo avec réduction de variance par conditionnement
    for (int i=0; i<100; i++)
    {
        Monte_Carlo_condi_spread(0.01, 2.0, 0.25, 0.3, 0.5, n_trajectoires[i], 1.2, 1.0, normale_c_r, 1.645, variables, 0.2);
        liste_prix[i] = variables[0];
        liste_variances[i] = variables[1];
        liste_i_down[i] = variables[2];
        liste_i_up[i] = variables[3];

    }
    ecriture("prix_variance_intervalles_90_spread_2", 100, n_trajectoires, liste_prix, liste_variances, liste_i_down, liste_i_up);

    //Option d'échange : Monte-Carlo avec réduction de variance par conditionnement
    for (int i=0; i<100; i++)
    {
        Monte_Carlo_condi_spread(0.01, 2.0, 0.25, 0.3, liste_rho[i], 17000, 1.2, 1.0, normale_c_r, 1.645, variables, 0.2);
        liste_prix[i] = variables[0];
    }
    ecriture("prix_fonction_rho_spread", 100, liste_rho, liste_prix);

    //Option d'échange : Monte-Carlo avec réduction de variance par conditionnement
    double liste_K[100];
    linspace(-1.0, 1.0, 100, liste_K);
    for (int i=0; i<100; i++)
    {
        Monte_Carlo_condi_spread(0.01, 2.0, 0.25, 0.3, 0.5, 17000, 1.2, 1.0, normale_c_r, 1.645, variables, liste_K[i]);
        liste_prix[i] = variables[0];
        liste_i_down[i] = variables[2];
        liste_i_up[i] = variables[3];
    }
    ecriture("prix_fonction_K_spread", 100, liste_K, liste_prix, liste_i_down, liste_i_up);


    //Option de spread : Monte-Carlo variable de contrôle
    for (int i=0; i<100; i++)
    {
        Monte_Carlo_controle_spread(0.01, 2.0, 0.25, 0.3, 0.5, n_trajectoires[i], 1.2, 1.0, normale_c_r, 1.645, variables, 0.2);
        liste_prix[i] = variables[0];
        liste_variances[i] = variables[1];
        liste_i_down[i] = variables[2];
        liste_i_up[i] = variables[3];

    }
    ecriture("prix_variance_intervalles_90_spread_controle", 100, n_trajectoires, liste_prix, liste_variances, liste_i_down, liste_i_up);

    double liste_sigma[30];
    linspace(0.0, 0.8, 30, liste_sigma);
    double** liste_prix_couple = new double*[30];  //matrice des valeurs en fonction de (sigma1,sigma2)
    for(int i = 0; i < 30; i++)
    {
        liste_prix_couple[i] = new double[30];
    }
    //Option d'échange : Monte-Carlo avec réduction de variance par conditionnement
    for (int i=0; i<30; i++)
    {
        for (int j=0; j<30; j++)
        {
            Monte_Carlo_condi_spread(0.01, 2.0, liste_sigma[i], liste_sigma[j], 0.5, 17000, 1.2, 1.0, normale_c_r, 1.645, variables, 0.2);
            liste_prix_couple[i][j] = variables[0];
        }
    }
    ecriture_couple("prix_fonction_sigma_spread", 30, liste_sigma, liste_prix_couple);
    for(int i = 0; i < 30; i++)
    {
        delete[] liste_prix_couple[i];
    }
    delete[] liste_prix_couple;

    Monte_Carlo_Best_Of(0.01, 2.0, 0.25, 0.3, 0.5, 1, 1.2, 1.0, normale_c_r, 1.645, variables, 1.0);

    //Option de spread : Monte-Carlo variable de contrôle
    for (int i=0; i<100; i++)
    {
        Monte_Carlo_controle_spread2(0.01, 2.0, 0.25, 0.3, 0.5, n_trajectoires[i], 1.2, 1.0, normale_c_r, 1.645, variables, 0.2);
        liste_prix[i] = variables[0];
        liste_variances[i] = variables[1];
        liste_i_down[i] = variables[2];
        liste_i_up[i] = variables[3];

    }
    ecriture("prix_variance_intervalles_90_spread_controle_2", 100, n_trajectoires, liste_prix, liste_variances, liste_i_down, liste_i_up);
}


