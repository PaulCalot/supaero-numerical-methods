%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROGRAMME PRINCIPAL - ELASTICITE 2D STATIONNAIRE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lecture et chargement du maillage %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Chargement du maillage

%load 'maillage_1.mat';
%load 'maillage_2.mat';
%load 'maillage_3.mat';
load 'maillage_trou.mat';

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
D = E/(1-v^2)*[1 v       0;
               v 1       0;
               0 0 (1-v)/2];

% Boucle sur les triangles
Elm_Mat_Cell = cell(1,Ntri);
A=zeros(2*Noeuds);

for K = 1:Ntri

    [sk, s1v, s2v, s3v] = Tri2Surf(K,tri,X,Y,Z);
    Bp = matBp(s1v,s2v,s3v);
    Ak = (1/(4*sk))*Bp'*D*Bp;
    Elm_Mat_Cell{K} = Ak;
    
    for i = 1:3
        for j = 1:3
        SKi = tri(K,i);
        SKj = tri(K,j);
        A(2*SKi-1, 2*SKj-1) = A(2*SKi-1, 2*SKj-1) + Ak(2*i-1, 2*j-1);
        A(2*SKi-1, 2*SKj) = A(2*SKi-1, 2*SKj) + Ak(2*i-1, 2*j);
        A(2*SKi, 2*SKj-1) = A(2*SKi, 2*SKj-1) + Ak(2*i, 2*j-1);
        A(2*SKi, 2*SKj) = A(2*SKi, 2*SKj) + Ak(2*i, 2*j);
       
        end % OF j FOR

    end  % OF i FOR
    
end % OF k FOR


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
for i = 1:Noeuds
    if isequal(X(i),0)
        A(2*i,:) = 0;
        A(2*i-1,:) = 0;
        A(:,2*i) = 0;
        A(:,2*i-1) = 0;
        A(2*i,2*i) = 1;
        A(2*i-1,2*i-1) = 1;       
    end
end


% Fin de boucle sur les noeuds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resolution de systeme et visualisation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Resolution

Sol=zeros(1,2*Noeuds);
Sol=A\G;

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


%% Fonction de creation du vecteur de sommet

function [Sk,s1v,s2v,s3v] = Tri2Surf(K,tri,X,Y,Z)
    S1 = tri(K,1);
    S2 = tri(K,2);
    S3 = tri(K,3);
    
    s1v = Sommet2Vect(S1,X,Y,Z);
    s2v = Sommet2Vect(S2,X,Y,Z);
    s3v = Sommet2Vect(S3,X,Y,Z);

    s12v = s1v - s2v;
    s13v = s1v - s3v;
    Sk = 0.5*norm(cross(s12v, s13v));
end


function [vect_S1] = Sommet2Vect(S,X,Y,Z)

vect_S1 = [X(S); Y(S); Z(S)];

end


function Bp = matBp(s1v,s2v,s3v)

x1 = s1v(1);
x2 = s2v(1);
x3 = s3v(1);

y1 = s1v(2);
y2 = s2v(2);
y3 = s3v(2);

z1 = s1v(3);
z2 = s2v(3);
z3 = s3v(3);

Bp = [y2 - y3,       0, y3 - y1,       0, y1 - y2,       0;
            0, x3 - x2,       0, x1 - x3,       0, x2 - x1;
      x3 - x2, y2 - y3, x1 - x3, y3 - y1, x2 - x1, y1 - y2];
 

end