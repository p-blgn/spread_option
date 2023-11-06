%Classique
fid=fopen("prix_variance_intervalles_90.txt",'r');
if fid <=0
   msg=['Le fichier : ' nomfile ' n a pas ete trouve'];
   error(msg);
end
n = 100;
chaine_pos = zeros(n,5);
status = 0;
for a=1:n
    ligne = str2num(fgetl(fid));
    chaine_pos(a,:) = ligne;
    status = feof(fid);
end
fclose(fid);
semilogx(chaine_pos(:,1),chaine_pos(:,2),'b');
hold on;
semilogx(chaine_pos(:,1),chaine_pos(:,3),'m');
hold on;
semilogx(chaine_pos(:,1),chaine_pos(:,4),'b-.');
hold on;
semilogx(chaine_pos(:,1),chaine_pos(:,5),'b-.');
hold on;

%Conditionnement
fid=fopen("prix_variance_intervalles_90_2.txt",'r');
if fid <=0
   msg=['Le fichier : ' nomfile ' n a pas ete trouve'];
   error(msg);
end
n = 100;
chaine_pos = zeros(n,5);
vrai_prix = zeros(n,1);
prix = Prix(0.25,0.3,0.5,1.0,1.0,2);
status = 0;
for a=1:n
    ligne = str2num(fgetl(fid));
    chaine_pos(a,:) = ligne;
    status = feof(fid);
    vrai_prix(a) = prix;
end
fclose(fid);
semilogx(chaine_pos(:,1),chaine_pos(:,2),'r');
hold on;
semilogx(chaine_pos(:,1),chaine_pos(:,3),'black');
hold on;
semilogx(chaine_pos(:,1),chaine_pos(:,4),'r-.');
hold on;
semilogx(chaine_pos(:,1),chaine_pos(:,5),'r-.');
hold on;
semilogx(chaine_pos(:,1),vrai_prix,'green');
title("Prix et variance option d'échange");
xlabel('Nombre de trajectoires') ;
ylabel('Variance/Prix');
legend({'P classique','Variance - Classique','IC 90% haut - Classique','IC 90% bas - Classique','P Condi','Variance - Condi','IC 90% haut - Condi','IC 90% bas - Condi','Prix théorique'},'Location','northeast');

figure;
%Rho
fid=fopen("prix_fonction_rho.txt",'r');
if fid <=0
   msg=['Le fichier : ' nomfile ' n a pas ete trouve'];
   error(msg);
end
n = 100;
chaine_pos = zeros(n,2);
status = 0;
valeurs2 = zeros(n,1);
for a=1:n
    ligne = str2num(fgetl(fid));
    chaine_pos(a,:) = ligne;
    status = feof(fid);
    valeurs2(a) = Prix(0.3,0.25,chaine_pos(a,1),1,1,2);
end
fclose(fid);
plot(chaine_pos(:,1),valeurs2,'r-.');
hold on;
plot(chaine_pos(:,1),chaine_pos(:,2),'b');
hold on;
title("Prix d'une option d'échange en fonction de \rho");
xlabel('\rho');
ylabel('Prix');
legend({'Prix calculé','Prix simulé'},'Location','northeast');


figure;
%Classique
fid=fopen("prix_variance_intervalles_90_spread.txt",'r');
if fid <=0
   msg=['Le fichier : ' nomfile ' n a pas ete trouve'];
   error(msg);
end
n = 100;
chaine_pos = zeros(n,5);
status = 0;
for a=1:n
    ligne = str2num(fgetl(fid));
    chaine_pos(a,:) = ligne;
    status = feof(fid);
end
fclose(fid);
chaine_classique = chaine_pos;
semilogx(chaine_pos(:,1),chaine_pos(:,2),'b');
hold on;
semilogx(chaine_pos(:,1),chaine_pos(:,3),'m');
hold on;
semilogx(chaine_pos(:,1),chaine_pos(:,4),'b-.');
hold on;
semilogx(chaine_pos(:,1),chaine_pos(:,5),'b-.');
hold on;

%Conditionnement
fid=fopen("prix_variance_intervalles_90_spread_2.txt",'r');
if fid <=0
   msg=['Le fichier : ' nomfile ' n a pas ete trouve'];
   error(msg);
end
n = 100;
chaine_pos = zeros(n,5);
status = 0;
for a=1:n
    ligne = str2num(fgetl(fid));
    chaine_pos(a,:) = ligne;
    status = feof(fid);
end
fclose(fid);
semilogx(chaine_pos(:,1),chaine_pos(:,2),'r');
hold on;
semilogx(chaine_pos(:,1),chaine_pos(:,3),'black');
hold on;
semilogx(chaine_pos(:,1),chaine_pos(:,4),'r-.');
hold on;
semilogx(chaine_pos(:,1),chaine_pos(:,5),'r-.');
title('Prix et variance option sur Spread');
xlabel('Nombre de trajectoires') ;
ylabel('Variance/Prix');
legend({'P classique','Variance - Classique','IC 90% haut - Classique','IC 90% bas - Classique','P Condi','Variance - Condi','IC 90% haut - Condi','IC 90% bas - Condi'},'Location','northeast');



figure;
%Rho
fid=fopen("prix_fonction_rho_spread.txt",'r');
if fid <=0
   msg=['Le fichier : ' nomfile ' n a pas ete trouve'];
   error(msg);
end
n = 100;
chaine_pos = zeros(n,2);
status = 0;
for a=1:n
    ligne = str2num(fgetl(fid));
    chaine_pos(a,:) = ligne;
    status = feof(fid);
end
fclose(fid);
plot(chaine_pos(:,1),chaine_pos(:,2),'b');
hold on;
title("Prix d'une option sur Spread en fonction de \rho");
xlabel('\rho');
ylabel('Prix');

figure;
%K
fid=fopen("prix_fonction_K_spread.txt",'r');
if fid <=0
   msg=['Le fichier : ' nomfile ' n a pas ete trouve'];
   error(msg);
end
n = 100;
chaine_pos = zeros(n,4);
status = 0;
for a=1:n
    ligne = str2num(fgetl(fid));
    chaine_pos(a,:) = ligne;
    status = feof(fid);
end
fclose(fid);
plot(chaine_pos(:,1),chaine_pos(:,2),'b');
hold on;
plot(chaine_pos(:,1),chaine_pos(:,3),'b-.');
hold on;
plot(chaine_pos(:,1),chaine_pos(:,4),'b-.');
title("Prix d'une option sur Spread en fonction de K");
xlabel('K') ;
ylabel('Prix');
legend({'Prix','IC 90% haut','IC 90% bas'},'Location','northeast');

figure;
%Controle
fid=fopen("prix_variance_intervalles_90_spread_controle.txt",'r');
if fid <=0
   msg=['Le fichier : ' nomfile ' n a pas ete trouve'];
   error(msg);
end
n = 100;
chaine_pos = zeros(n,5);
status = 0;
for a=1:n
    ligne = str2num(fgetl(fid));
    chaine_pos(a,:) = ligne;
    status = feof(fid);
end
fclose(fid);
semilogx(chaine_pos(:,1),chaine_pos(:,2),'r');
hold on;
semilogx(chaine_pos(:,1),chaine_pos(:,3),'black');
hold on;
semilogx(chaine_pos(:,1),chaine_pos(:,4),'r-.');
hold on;
semilogx(chaine_pos(:,1),chaine_pos(:,5),'r-.');
hold on;
semilogx(chaine_classique(:,1),chaine_classique(:,2),'b');
hold on;
semilogx(chaine_classique(:,1),chaine_classique(:,3),'m');
hold on;
semilogx(chaine_classique(:,1),chaine_classique(:,4),'b-.');
hold on;
semilogx(chaine_classique(:,1),chaine_classique(:,5),'b-.');
title("Prix et variance option sur Spread - Monte-Carlo contrôle avec Put");
xlabel('Nombre de trajectoires') ;
ylabel('Variance/Prix');
legend({'P Contrôle','Variance - Contrôle','IC 90% haut - Contrôle','IC 90% bas - Contrôle','P Classique','Variance - Classique','IC 90% haut - Classique','IC 90% bas - Classique'},'Location','northeast');



figure;
%sigma
fid=fopen("prix_fonction_sigma_spread.txt",'r');
if fid <=0
   msg=['Le fichier : ' nomfile ' n a pas ete trouve'];
   error(msg);
end
n = 30;
chaine_pos = zeros(n,90);
status = 0;
for a=1:n
    ligne = str2num(fgetl(fid));
    chaine_pos(a,:) = ligne;
    status = feof(fid);
end
X = zeros(1,n);
Y = zeros(n,n);
for a=1:n
    X(a) = chaine_pos(a,1);
    for b=1:n
        Y(b,a) = chaine_pos(a,3*b);
    end
end
fclose(fid);
surf(X,X,Y);
xlabel('\sigma_1');
ylabel('\sigma_2');
zlabel('Prix');
title("Prix en fonction des volatilités d'une option sur Spread");

figure;
%Controle
fid=fopen("prix_variance_intervalles_90_spread_controle_2.txt",'r');
if fid <=0
   msg=['Le fichier : ' nomfile ' n a pas ete trouve'];
   error(msg);
end
n = 100;
chaine_pos = zeros(n,5);
status = 0;
for a=1:n
    ligne = str2num(fgetl(fid));
    chaine_pos(a,:) = ligne;
    status = feof(fid);
end
fclose(fid);
semilogx(chaine_pos(:,1),chaine_pos(:,2),'r');
hold on;
semilogx(chaine_pos(:,1),chaine_pos(:,3),'black');
hold on;
semilogx(chaine_pos(:,1),chaine_pos(:,4),'r-.');
hold on;
semilogx(chaine_pos(:,1),chaine_pos(:,5),'r-.');
hold on;
semilogx(chaine_classique(:,1),chaine_classique(:,2),'b');
hold on;
semilogx(chaine_classique(:,1),chaine_classique(:,3),'m');
hold on;
semilogx(chaine_classique(:,1),chaine_classique(:,4),'b-.');
hold on;
semilogx(chaine_classique(:,1),chaine_classique(:,5),'b-.');
title("Prix et variance option sur Spread - Monte-Carlo contrôle avec option d'échange");
xlabel('Nombre de trajectoires'); 
ylabel('Variance/Prix');
legend({'P Contrôle','Variance - Contrôle','IC 90% haut - Contrôle','IC 90% bas - Contrôle','P Classique','Variance - Classique','IC 90% haut - Classique','IC 90% bas - Classique'},'Location','northeast');

Prix(0.25,0.3,0.5,1.2,1,2) - exp(-0.02) + 1; %Prix théorique du Best-Of