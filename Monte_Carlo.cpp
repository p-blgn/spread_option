#include "Monte_Carlo.hpp"
#include <iostream>
#include <string>
#include <fstream>


c_r::c_r(){}
c_r::~c_r(){}
double c_r::simuler() //tirage d'une variable aléatoire suivant une N(0,1) à l'aide de l'algorithme de Box-Müller
{
    return sqrt(-2.0*log(1-dist(gen)))*cos(2.0*M_PI*dist(gen));
}
double normale(double moyenne, double variance, c_r& normale_c_r) //tirage d'une variable aléatoire suivant une N(moyenne,variance)
{
    return normale_c_r.simuler()*variance + moyenne;
}
void W::constantes() //calcule les constantes utiles (rac_dt, lambda1, lambda2)
{
    this->rac_dt = sqrt(this->dt);
    this->lambda1 = this->rac_dt*this->rho;
    this->lambda2 = this->rac_dt*sqrt(1-this->rho*this->rho);
}
void W::iterer() //effectue un tirage de W
{
    double U1 = normale(0.0,1.0,loi_normale);
    double U2 = normale(0.0,1.0,loi_normale);
    valeurs[0] = rac_dt*U1;
    valeurs[1] = lambda1*U1+lambda2*U2;
}
double N(double x) //calcule une approximation de Abramowitz & Stegun de N(x) où N est la fonction de répartition de la loi normale centrée réduite
{
    double b0 = 0.2316419;
    double b1 = 0.31938153;
    double b2 = -0.356563782;
    double b3 = 1.781477937;
    double b4 = -1.821255978;
    double b5 = 1.330274429;
    double t;
    if (x>=0)
    {
        t = 1.0/(1.0+b0*x);
        return 1.0-1/sqrt(2*M_PI)*exp(-x*x/2)*(b1*t+b2*t*t+b3*t*t*t+b4*t*t*t*t+b5*t*t*t*t*t);
    }
    else //on utilise le fait que N(-x)=1-N(x)
    {
        return 1.0-N(-x);
    }
}
double Put_price(double K, double S0, double sigma, double r, double T)
//calcul direct du prix du Put avec une formule fermée
{
    if (K<=0.0)
    {
        return 0.0;
    }
    else
    {
        double d1 = (log(K/S0)-(r-sigma*sigma/2)*T)/(sigma*sqrt(T));
        double d2 = d1-sigma*sqrt(T);
        double N1 = N(d1);
        double N2 = N(d2);
        return K*N1-S0*N2*exp(r*T);
    }
}
void Monte_Carlo_classique_echange(double r, double T, double sigma1, double sigma2, double rho, int n_tirages, double alphaS10, double betaS20, c_r& normale_c_r, double coefficient, double* variables)
{
    //Option d'échange : Monte-Carlo classique
    W mouvement; //brownien
    mouvement.loi_normale = normale_c_r; //Copie de la loi normale
    mouvement.rho = rho;
    mouvement.dt = T;
    mouvement.constantes();
    mouvement.valeurs[0] = 0.0;
    mouvement.valeurs[1] = 0.0;
    double S1;
    double S2;
    double somme = 0.0;
    double* valeurs = (double*) malloc(sizeof(double)*n_tirages); //contiendra les résultats des tirages
    if (valeurs==nullptr)
    {
        cout<<"Erreur d'allocation\n";
    }
    double Var_est = 0.0;
    double Moy_est = 0.0;
    double RMSE;
    for (int i=0; i<n_tirages; i++)
    {
        mouvement.iterer();
        S1 = alphaS10*exp((r-sigma1*sigma1/2)*T+sigma1*mouvement.valeurs[0]);
        S2 = betaS20*exp((r-sigma2*sigma2/2)*T+sigma2*mouvement.valeurs[1]);
        valeurs[i] = max(S1-S2,0.0);
        somme = somme + valeurs[i];
    }
    Moy_est = somme/n_tirages;
    for (int i=0; i<n_tirages; i++)
    {
        Var_est = Var_est + pow(Moy_est-valeurs[i],2);
    }
    Var_est = exp(-2*r*T)*Var_est/(n_tirages-1);
    RMSE = sqrt(Var_est)/sqrt(n_tirages);
    Moy_est = exp(-r*T)*Moy_est;
    variables[0] = Moy_est;
    variables[1] = Var_est;
    variables[2] = Moy_est-coefficient*RMSE;
    variables[3] = Moy_est+coefficient*RMSE;
    normale_c_r.gen = mouvement.loi_normale.gen; //comme on n'a fait une copie de la loi normale sa seed est toujours la même, on la remplace par la seed actuelle pour éviter de retirer les mêmes valeurs à chaque simulation
    free(valeurs);
}
void Monte_Carlo_condi_echange(double r, double T, double sigma1, double sigma2, double rho, int n_tirages, double alphaS10, double betaS20, c_r& normale_c_r, double coefficient, double* variables)
{
    //Option d'échange : Monte-Carlo avec réduction de variance par conditionnement
    W mouvement;
    mouvement.loi_normale = normale_c_r; //Copie de la loi normale
    mouvement.rho = rho;
    mouvement.dt = T;
    mouvement.constantes();
    mouvement.valeurs[0] = 0.0;
    mouvement.valeurs[1] = 0.0;
    double somme = 0.0;
    double Var_est = 0.0;
    double* valeurs = (double*) malloc(sizeof(double)*n_tirages);
    if (valeurs==nullptr)
    {
        cout<<"Erreur d'allocation\n";
    }
    double Moy_est = 0.0;
    double RMSE;
    double sigma = sigma2*sqrt(1.0-mouvement.rho*mouvement.rho);
    double K;
    double w;
    double S0;
    for (int i=0; i<n_tirages; i++)
    {
        mouvement.iterer();
        w = mouvement.valeurs[0];
        K = alphaS10*exp((r-sigma1*sigma1/2)*T+sigma1*w); //on se ramène à un put classique de strike K
        S0 = betaS20*exp(sigma2*mouvement.rho*w+(sigma*sigma-sigma2*sigma2)/2*T);
        valeurs[i] = Put_price(K,S0,sigma,r,T);
        somme = somme + valeurs[i];
    }
    Moy_est = somme/n_tirages;
    for (int i=0; i<n_tirages; i++)
    {
        Var_est = Var_est + pow(Moy_est-valeurs[i],2);
    }
    Var_est = exp(-2*r*T)*Var_est/(n_tirages-1);
    RMSE = sqrt(Var_est)/sqrt(n_tirages);
    Moy_est = exp(-r*T)*Moy_est;
    variables[0] = Moy_est;
    variables[1] = Var_est;
    variables[2] = Moy_est-coefficient*RMSE;
    variables[3] = Moy_est+coefficient*RMSE;
    normale_c_r.gen = mouvement.loi_normale.gen; //comme on n'a fait une copie de la loi normale sa seed est toujours la même, on la remplace par la seed actuelle pour éviter de retirer les mêmes valeurs à chaque simulation
    free(valeurs);
}
void Monte_Carlo_classique_spread(double r, double T, double sigma1, double sigma2, double rho, int n_tirages, double alphaS10, double betaS20, c_r& normale_c_r, double coefficient, double* variables, double K)
{
    //Option de spread : Monte-Carlo classique
    W mouvement;
    mouvement.loi_normale = normale_c_r; //Copie de la loi normale
    mouvement.rho = rho;
    mouvement.dt = T;
    mouvement.constantes();
    mouvement.valeurs[0] = 0.0;
    mouvement.valeurs[1] = 0.0;
    double somme = 0.0;
    double Var_est = 0.0;
    double* valeurs = (double*) malloc(sizeof(double)*n_tirages);
    if (valeurs==nullptr)
    {
        cout<<"Erreur d'allocation\n";
    }
    double S1;
    double S2;
    double Moy_est = 0.0;
    double RMSE;
    for (int i=0; i<n_tirages; i++)
    {
        mouvement.iterer();
        S1 = alphaS10*exp((r-sigma1*sigma1/2)*T+sigma1*mouvement.valeurs[0]);
        S2 = betaS20*exp((r-sigma2*sigma2/2)*T+sigma2*mouvement.valeurs[1]);
        valeurs[i] = max(S1-S2-K,0.0);
        somme = somme + valeurs[i];
    }
    Moy_est = somme/n_tirages;
    for (int i=0; i<n_tirages; i++)
    {
        Var_est = Var_est + pow(Moy_est-valeurs[i],2);
    }
    Var_est = exp(-2*r*T)*Var_est/(n_tirages-1);
    RMSE = sqrt(Var_est)/sqrt(n_tirages);
    Moy_est = exp(-r*T)*Moy_est;
    variables[0] = Moy_est;
    variables[1] = Var_est;
    variables[2] = Moy_est-coefficient*RMSE;
    variables[3] = Moy_est+coefficient*RMSE;
    normale_c_r.gen = mouvement.loi_normale.gen; //comme on n'a fait une copie de la loi normale sa seed est toujours la même, on la remplace par la seed actuelle pour éviter de retirer les mêmes valeurs à chaque simulation
    free(valeurs);
}
void Monte_Carlo_condi_spread(double r, double T, double sigma1, double sigma2, double rho, int n_tirages, double alphaS10, double betaS20, c_r& normale_c_r, double coefficient, double* variables, double K)
{
    //Option sur Spread : Monte-Carlo avec réduction de variance par conditionnement
    W mouvement;
    mouvement.loi_normale = normale_c_r; //Copie de la loi normale
    mouvement.rho = rho;
    mouvement.dt = T;
    mouvement.constantes();
    mouvement.valeurs[0] = 0.0;
    mouvement.valeurs[1] = 0.0;
    double somme = 0.0;
    double Var_est = 0.0;
    double* valeurs = (double*) malloc(sizeof(double)*n_tirages);
    if (valeurs==nullptr)
    {
        cout<<"Erreur d'allocation\n";
    }
    double Moy_est = 0.0;
    double RMSE;
    double sigma = sigma2*sqrt(1.0-mouvement.rho*mouvement.rho);
    double K2;
    double w;
    double S0;
    for (int i=0; i<n_tirages; i++)
    {
        mouvement.iterer();
        w = mouvement.valeurs[0];
        K2 = alphaS10*exp((r-sigma1*sigma1/2)*T+sigma1*w)-K;
        S0 = betaS20*exp(sigma2*mouvement.rho*w+(sigma*sigma-sigma2*sigma2)/2*T);
        valeurs[i] = Put_price(K2,S0,sigma,r,T);
        somme = somme + valeurs[i];
    }
    Moy_est = somme/n_tirages;
    for (int i=0; i<n_tirages; i++)
    {
        Var_est = Var_est + pow(Moy_est-valeurs[i],2);
    }
    Var_est = exp(-2*r*T)*Var_est/(n_tirages-1);
    RMSE = sqrt(Var_est)/sqrt(n_tirages);
    Moy_est = exp(-r*T)*Moy_est;
    variables[0] = Moy_est;
    variables[1] = Var_est;
    variables[2] = Moy_est-coefficient*RMSE;
    variables[3] = Moy_est+coefficient*RMSE;
    normale_c_r.gen = mouvement.loi_normale.gen; //comme on n'a fait une copie de la loi normale sa seed est toujours la même, on la remplace par la seed actuelle pour éviter de retirer les mêmes valeurs à chaque simulation
    free(valeurs);
}

void Monte_Carlo_controle_spread(double r, double T, double sigma1, double sigma2, double rho, int n_tirages, double alphaS10, double betaS20, c_r& normale_c_r, double coefficient, double* variables, double K)
{
    W mouvement;
    mouvement.loi_normale = normale_c_r; //Copie de la loi normale
    mouvement.rho = rho;
    mouvement.dt = T;
    mouvement.constantes();
    mouvement.valeurs[0] = 0.0;
    mouvement.valeurs[1] = 0.0;
    mouvement.rho = rho;
    mouvement.dt = T;
    mouvement.constantes();
    double somme = 0.0; //somme liée à la variable à estimer
    double somme2 = 0.0; //somme liée à la variable de contrôle
    double Var_est = 0.0;
    double covar = 0.0;
    double Var_est2 = 0.0;
    double c_opt;
    double somme3 = 0.0;
    double Var_est3 = 0.0;
    double Moy_est3;
    double correl;
    double S1;
    double S2;
    double RMSE;
    double Moy_est;
    double* valeurs = (double*) malloc(sizeof(double)*n_tirages);
    if (valeurs==nullptr)
    {
        cout<<"Erreur d'allocation\n";
    }
    for (int i=0; i<n_tirages; i++)
    {
        mouvement.iterer();
        S1 = alphaS10*exp((r-sigma1*sigma1/2)*T+sigma1*mouvement.valeurs[0]);
        S2 = betaS20*exp((r-sigma2*sigma2/2)*T+sigma2*mouvement.valeurs[1]);
        valeurs[i] = S1-S2-K;
        somme2 = somme2 + valeurs[i];
        if (valeurs[i]>=0)
        {
            somme = somme + valeurs[i];
        }
    }
    Moy_est = somme/n_tirages;
    for (int i=0; i<n_tirages; i++)
    {
        Var_est = Var_est + pow(max(valeurs[i],0.0)-Moy_est,2);
        Var_est2 = Var_est2 + pow(exp(r*T)*(alphaS10-betaS20)-K-valeurs[i],2);
        covar = covar + (Moy_est-max(valeurs[i],0.0))*(exp(r*T)*(alphaS10-betaS20)-K-valeurs[i]);
    }
    Var_est = exp(-2*r*T)*Var_est/(n_tirages-1);
    Var_est2 = exp(-2*r*T)*Var_est2/(n_tirages-1);
    covar = exp(-2*r*T)*covar/(n_tirages-1);
    correl = covar/sqrt(Var_est*Var_est2);
    c_opt = sqrt(Var_est/Var_est2)*correl;
    Moy_est = exp(-r*T)*Moy_est;
    for (int i=0; i<n_tirages; i++)
    {
        valeurs[i] = max(valeurs[i],0.0)-c_opt*(valeurs[i]-(exp(r*T)*(alphaS10-betaS20)-K));
        somme3 = somme3 + valeurs[i];
    }
    Moy_est3 = somme3/n_tirages;
    for (int i=0; i<n_tirages; i++)
    {
        Var_est3 = Var_est3 + pow(Moy_est3-valeurs[i],2);
    }
    Var_est3 = exp(-2*r*T)*Var_est3/(n_tirages-1);
    Moy_est3 = exp(-r*T)*Moy_est3;
    RMSE = sqrt(Var_est3)/sqrt(n_tirages);
    double variance_calcul = (1-correl*correl)*Var_est;
    variables[0] = Moy_est3;
    variables[1] = Var_est3;
    variables[2] = Moy_est3-coefficient*RMSE;
    variables[3] = Moy_est3+coefficient*RMSE;
    normale_c_r.gen = mouvement.loi_normale.gen;
    free(valeurs);
}

void Monte_Carlo_controle_spread2(double r, double T, double sigma1, double sigma2, double rho, int n_tirages, double alphaS10, double betaS20, c_r& normale_c_r, double coefficient, double* variables, double K)
{
    W mouvement;
    mouvement.loi_normale = normale_c_r; //Copie de la loi normale
    mouvement.rho = rho;
    mouvement.dt = T;
    mouvement.constantes();
    mouvement.valeurs[0] = 0.0;
    mouvement.valeurs[1] = 0.0;
    mouvement.rho = rho;
    mouvement.dt = T;
    mouvement.constantes();
    double somme = 0.0;
    double somme2 = 0.0;
    double Var_est = 0.0;
    double covar = 0.0;
    double Var_est2 = 0.0;
    double c_opt;
    double somme3 = 0.0;
    double Var_est3 = 0.0;
    double Moy_est3;
    double correl;
    double RMSE;
    double S1;
    double S2;
    double Moy_est;
    double Moy_est2;
    double* valeurs = (double*) malloc(sizeof(double)*n_tirages);
    double* valeurs2 = (double*) malloc(sizeof(double)*n_tirages);
    double sigma = sqrt(sigma1*sigma1 + sigma2*sigma2 - 2*sigma1*sigma2*rho);
    double d1 = 1.0/(sigma*sqrt(T))*(log(alphaS10/betaS20)+sigma*sigma*T/2);
    double d2 = d1 - sigma*sqrt(T);
    double esperance = alphaS10*N(d1) - betaS20*N(d2);
    if (valeurs==nullptr)
    {
        cout<<"Erreur d'allocation\n";
    }
    for (int i=0; i<n_tirages; i++)
    {
        mouvement.iterer();
        S1 = alphaS10*exp((r-sigma1*sigma1/2)*T+sigma1*mouvement.valeurs[0]);
        S2 = betaS20*exp((r-sigma2*sigma2/2)*T+sigma2*mouvement.valeurs[1]);
        valeurs[i] = max(S1-S2-K,0.0);
        valeurs2[i] = max(S1-S2,0.0);
        somme2 = somme2 + valeurs2[i];
        somme = somme + valeurs[i];
    }
    Moy_est = somme/n_tirages;
    Moy_est2 = somme2/n_tirages;
    for (int i=0; i<n_tirages; i++)
    {
        Var_est = Var_est + pow(valeurs[i]-Moy_est,2);
        Var_est2 = Var_est2 + pow(valeurs2[i]-Moy_est2,2);
        covar = covar + (Moy_est-valeurs[i])*(Moy_est2-valeurs2[i]);
    }
    Var_est = exp(-2*r*T)*Var_est/(n_tirages-1);
    Var_est2 = exp(-2*r*T)*Var_est2/(n_tirages-1);
    covar = exp(-2*r*T)*covar/(n_tirages-1);
    correl = covar/sqrt(Var_est*Var_est2);
    c_opt = sqrt(Var_est/Var_est2)*correl;
    Moy_est = exp(-r*T)*Moy_est;
    Moy_est2 = exp(-r*T)*Moy_est2;
    for (int i=0; i<n_tirages; i++)
    {
        valeurs[i] = valeurs[i]-c_opt*(valeurs2[i]-esperance*exp(r*T));
        somme3 = somme3 + valeurs[i];
    }
    Moy_est3 = somme3/n_tirages;
    for (int i=0; i<n_tirages; i++)
    {
        Var_est3 = Var_est3 + pow(Moy_est3-valeurs[i],2);
    }
    Var_est3 = exp(-2*r*T)*Var_est3/(n_tirages-1);
    Moy_est3 = exp(-r*T)*Moy_est3;
    RMSE = sqrt(Var_est3)/sqrt(n_tirages);
    double variance_calcul = (1-correl*correl)*Var_est;
    normale_c_r.gen = mouvement.loi_normale.gen;
    variables[0] = Moy_est3;
    variables[1] = Var_est3;
    variables[2] = Moy_est3-coefficient*RMSE;
    variables[3] = Moy_est3+coefficient*RMSE;
    free(valeurs);
    free(valeurs2);
}
void ecriture(string nom_fichier, int n_points, int* X, double* Y, double* Y2, double* Y3, double* Y4)
//permet d'écrire dans un fichier .txt des couples (int,double) (ou des triplets (int, double, double)...)
{
    string const chemin("C:/Users/pierr/Documents/Cours/PRB222/Spread option/" + nom_fichier + ".txt");
    ofstream monFlux(chemin.c_str());
    if(monFlux)
    {
    string chaine;
    for(int i = 0 ; i < n_points ; i++)
    {
        chaine += to_string(X[i]);
        chaine += " ";
        chaine += to_string(Y[i]);
        chaine += " ";
        if (Y2!=nullptr)
        {
            chaine += to_string(Y2[i]);
            chaine += " ";
        }
        if (Y3!=nullptr)
        {
            chaine += to_string(Y3[i]);
            chaine += " ";
        }
        if (Y4!=nullptr)
        {
            chaine += to_string(Y4[i]);
            chaine += " ";
        }
        chaine += "\n";
    }
    monFlux << chaine << endl ;
    }
    else
    cout << "Impossible d'ouvrir le fichier." << endl;
}
void ecriture(string nom_fichier, int n_points, double* X, double* Y, double* Y2, double* Y3, double* Y4)
//permet d'écrire dans un fichier .txt des couples (double,double) (ou des triplets (double, double, double)...)
{
    string const chemin("DossierSortie" + nom_fichier + ".txt");
    ofstream monFlux(chemin.c_str());
    if(monFlux)
    {
    string chaine;//la chaîne contenant le texte à inscrire dans le fichier
    for(int i = 0 ; i < n_points ; i++)
    {
        chaine += to_string(X[i]);
        chaine += " ";
        chaine += to_string(Y[i]);
        chaine += " ";
        if (Y2!=nullptr)
        {
            chaine += to_string(Y2[i]);
            chaine += " ";
        }
        if (Y3!=nullptr)
        {
            chaine += to_string(Y3[i]);
            chaine += " ";
        }
        if (Y4!=nullptr)
        {
            chaine += to_string(Y4[i]);
            chaine += " ";
        }
        chaine += "\n"; //on écrit un résultat par ligne
    }
    monFlux << chaine << endl;
    }
    else
    cout << "Impossible d'ouvrir le fichier." << endl;
}
void linspace(double debut, double fin, int n, double *valeurs)
//crée un ensemble de valeurs de type double espacées de manière linéaire (équivalent au linspace Python)
{
    double pas = (fin-debut)/(n-1);
    for (int i=0; i<n; i++)
    {
        valeurs[i] = i*pas + debut;
    }
}
void linspace(int debut, int fin, int n, int *valeurs)
//crée un ensemble de valeurs de type int espacées de manière linéaire (équivalent au linspace Python)
{
    double pas = (fin-debut)*1.0/(n-1);
    for (int i=0; i<n; i++)
    {
        valeurs[i] = i*pas + debut;
    }
}
void logspace(double debut, double fin, int n, double *valeurs)
//crée un ensemble de valeurs de type double espacées de manière exponentielle (selon une échelle log) (équivalent au numpy.logspace Python)
{
    double pas = log(fin/debut)/(n-1);
    for (int i=0; i<n; i++)
    {
        valeurs[i] = debut*exp(i*pas);
    }
}
void logspace(int debut, int fin, int n, int *valeurs)
//crée un ensemble de valeurs de type int espacées de manière exponentielle (selon une échelle log) (équivalent au numpy.logspace Python)
{
    double pas = log(fin*1.0/debut)/(n-1);
    for (int i=0; i<n; i++)
    {
        valeurs[i] = debut*exp(i*pas);
    }
}
void ecriture_couple(string nom_fichier, int n_points, double* X, double** Y)
//Permet d'écrire dans un fichier .txt des valeurs de fonctiondes de R^2 dans R (ici les valeurs en abscisse sont dans A^2 avec A un certain ensemble
{
    string const chemin("DossierSortie" + nom_fichier + ".txt");
    ofstream monFlux(chemin.c_str());
    if(monFlux)
    {
    string chaine;
    for(int i = 0 ; i < n_points; i++)
    {
        for (int j=0; j<n_points; j++)
        {
            chaine += to_string(X[i]);
            chaine += " ";
            chaine += to_string(X[j]);
            chaine += " ";
            chaine += to_string(Y[i][j]);
            chaine += " "; //on écrit toutes les valeurs ayant la même abscisse sur une même ligne
        }
        chaine += "\n";
    }
    monFlux << chaine << endl;
    }
    else
    cout << "Impossible d'ouvrir le fichier." << endl;
}
void Monte_Carlo_Best_Of(double r, double T, double sigma1, double sigma2, double rho, int n_tirages, double alphaS10, double betaS20, c_r& normale_c_r, double coefficient, double* variables, double K)
{
    //Option sur Best Of : Monte-Carlo classique
    cout<<"######################"<<"\n";
    cout<<"Monte-Carlo classique Best Of"<<"\n";
    cout<<"######################"<<"\n";
    W mouvement;
    mouvement.loi_normale = normale_c_r; //Copie de la loi normale
    mouvement.rho = rho;
    mouvement.dt = T;
    mouvement.constantes();
    mouvement.valeurs[0] = 0.0;
    mouvement.valeurs[1] = 0.0;
    double somme = 0.0;
    double Var_est = 0.0;
    double* valeurs = (double*) malloc(sizeof(double)*n_tirages);
    if (valeurs==nullptr)
    {
        cout<<"Erreur d'allocation\n";
    }
    double S1;
    double S2;
    double Moy_est = 0.0;
    double RMSE;
    for (int i=0; i<n_tirages; i++)
    {
        mouvement.iterer();
        S1 = alphaS10*exp((r-sigma1*sigma1/2)*T+sigma1*mouvement.valeurs[0]);
        S2 = betaS20*exp((r-sigma2*sigma2/2)*T+sigma2*mouvement.valeurs[1]);
        valeurs[i] = max(S1,S2) - K;
        somme = somme + valeurs[i];
    }
    Moy_est = somme/n_tirages;
    for (int i=0; i<n_tirages; i++)
    {
        Var_est = Var_est + pow(Moy_est-valeurs[i],2);
    }
    Var_est = exp(-2*r*T)*Var_est/(n_tirages-1);
    RMSE = sqrt(Var_est)/sqrt(n_tirages);
    Moy_est = exp(-r*T)*Moy_est;
    variables[0] = Moy_est;
    variables[1] = Var_est;
    variables[2] = Moy_est-coefficient*RMSE;
    variables[3] = Moy_est+coefficient*RMSE;
    cout<<"Moyenne = "<<Moy_est<<"\n";
    cout<<"Variance = "<<Var_est<<"\n";
    cout<<"Intervalle 95% = ["<<Moy_est-1.96*RMSE<<","<<Moy_est+1.96*RMSE<<"]"<<"\n";
    cout<<"Nombre de trajectoires = "<<n_tirages<<"\n";
    normale_c_r.gen = mouvement.loi_normale.gen;
    free(valeurs);
}
