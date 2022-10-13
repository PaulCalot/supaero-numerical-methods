%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROGRAMME PRINCIPAL - ELASTICITE 2D STATIONNAIRE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lecture et chargement du maillage %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Chargement du maillage

load 'maillage_1.mat';
%load 'maillage_2.mat';
%load 'maillage_3.mat';
%load 'maillage_trou.mat';

% Coordonnees des noeuds du maillage

X=p(1,:);
Y=p(2,:);

% Nombre de noeuds

Noeuds=length(X);

% Nombre de triangles

Ntri=length(t(1,:));

% Initialisation de Z

Z=zeros(1,Noeuds);

% Visualisation du maillage en 2D (Z = Z(X,Y))

tri=t(1:3,:)';
figure;
trisurf(tri,X,Y,Z);
view(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcul de la matrice de rigidite %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parametres mecaniques de la plaque

E = 1000.;     % Module d'Young
v = 0.3;       % Coefficient de Poisson

% Matrice D de la loi de comportement


% Boucle sur les triangles


% Fin de boucle sur les triangles

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcul du second membre %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialisation du second membre

G=zeros(2*Noeuds,1);

% Valeur du chargement (force concentree) sur le bord

gx = 10.0;
gy = 0.0;

% Boucle sur les noeuds recevant des forces concentrees             

for i=1:Noeuds
    
    if X(i)==2      % Arete recevant des forces concentrees
        
        G(2*i-1)=gx;
        G(2*i)=gy;
        
        if Y(i)==0  % Cas du premier coin (bas)
            G(2*i-1)=gx/2;
            G(2*i)=gy/2;
        end
        
        if Y(i)==1  % Cas du second coin (haut)
            G(2*i-1)=gx/2;
            G(2*i)=gy/2;
        end
        
    end
    
end

% Fin de boucle sur les noeuds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Imposition des conditions aux limites de Dirichlet %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Boucle sur les noeuds bloques



% Fin de boucle sur les noeuds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resolution de systeme et visualisation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Resolution

% Sol=A\G;
Sol=zeros(1,2*Noeuds);

% Calcul de la position des noeuds apres deplacement

for i=1:Noeuds
    XX(i)=X(i)+Sol(2*i-1);
    YY(i)=Y(i)+Sol(2*i);
end

% Visualisation 2D du maillage deforme

figure;
Z=zeros(1,Noeuds);
trisurf(tri,XX,YY,Z);
view(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcul des contraintes de Von Mises %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Boucle sur les triangles

sigmaVM=zeros(1,Noeuds);

% Visualisation des contraintes

figure;
trisurf(tri,XX,YY,sigmaVM)
view(2);
