clear all
close all 
format long
addpath('POINTS_DE_HESTHAVEN')
addpath('MAILLAGES_TP')
addpath('QUADRATURES')
addpath('OUTILS')
addpath('MATRICES_ELEMENTAIRES')

%%%%%%%%%%%%%%%%%%%%%
%PARAMETRES PHYSIQUES
%%%%%%%%%%%%%%%%%%%%%

%paramètres du vide
c0   = 3*1E8;
mu0  = 4*pi*1E-7;
eps0 = 1/mu0/c0^2;
Z0   = sqrt(mu0/eps0);
Y0   = 1/Z0;

%permittivité et perméabilité relatives pour l'exemple hétérogène
eps_diel = 1;
mu_diel  = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%
%PARAMETRES DE SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%
fich_mesh = 'cavity_pas03.mat';
%fich_mesh = 'cavite_ptsource_pas01.mat';
ordre   = 3;          %ordre de l'approximation polynomiale
Tf      = 1E-9;          %Temps de simulation
isource = 0;        %=0 si source modale et =1 si source ponctelle
penE    = 0*1/(2*Z0)*1;  %pénalisation pour l'éqution en E
penH    = 0*1/(2*Y0)*1;  %pénalisation pour l'éqution en H

CL_SM   = 0;          %Parametre pour prendre en compte (=1) ou pas (=0) une condition aux limites de Silver-Muller sur le bord extérieur
nref    = [1/sqrt(2) -1 0;1/sqrt(2) 0 -1];%Normales unitaires sortantes sur le triangle de référence

%Formules de quadratures
ordreq           = 2*ordre; %ordre de la quadrature sur les triangles
[Nq,xq,yq,wq]    = quadrature_triangle_reference(ordreq); 
[Ng1D,xg1D,wg1D] = quadrature_segment_reference(ordre);

%Parametres de post-traitement
p_visu    = 1;
mesh_ref  = creation_maillage_raffine_reference(p_visu);%redécoupage du triangle de référence
flag_film = 0; % =0 pas d animation et =1 animation
Xmin = 0;  Xmax = 1;
Ymin = 0;  Ymax = 1;
Zmin = -1; Zmax = 1;

%%%%%%%%%%%%%%%%%%%%%%%%
%%Definition des sources 
%%%%%%%%%%%%%%%%%%%%%%%%
if isource == 0 
%Source modale
mx = 1; my=1;
else 
%Source ponctuelle
xs      = 0.5;
ys      = 0.5;
rs      = 0.2;
t0      = 1E-10;
%t0 = 5E-9/40;
sigma_t = 2E-12;
end
 



 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Determination des (ordre+1)*(ordre+2)/2 points d'interpolation et 
%fonctions de base sur le triangle de reference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[N,x,y,A] = ddl_ET_fct_base_SUR_triangle_reference(ordre);


%%%%%%%%%%%%%%%%%%%%%
%Lecture du maillage
%%%%%%%%%%%%%%%%%%%%%
[mesh,Type_arete] = lecture_mesh(fich_mesh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ASSEMBLAGE MATRICES DE MASSE, RIGIDITE ET SAUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Preparation de la vectorisation des matrices élémentaires 
%Mat(1:P,1:Q) devient Mat_V(1:P*Q) + ii_Mat(1:P*Q) +jj_Mat(1:P*Q)
%%ii et jj pour les matrices NxN
%%iiv et jjv pour les matrices 2Nx2N
%%iir et jjr pour les matrices 2NxN
[ii,jj,iiv,jjv,iir,jjr]=indices_vectorisation_matrices_elementaires(N);

[kk,kkv,kkr,kkpe,kkph,kk_SM] = initialisation_indices_avancee_assemblage(N);
[NNZ_IME,NNZ_IMH,NNZ_R,NNZ_P_E,NNZ_P_H,NNZ_SM] = NB_NON_ZEROS_MATRICES(N,mesh,Type_arete,CL_SM);


Nddl_E = 2*N*mesh.Nbtri; %Nbre de degrés de liberté vectoriels
Nddl_H = N*mesh.Nbtri;   %Nbre de degrés de liberté scalaires


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation des matrices sparse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inverse matrice de masse H (scalaire)
I_IMH = zeros(NNZ_IMH,1);J_IMH = zeros(NNZ_IMH,1);IMH  = zeros(NNZ_IMH,1);
%Inverse matrice de masse E (vectorielle)
I_IME = zeros(NNZ_IME,1);J_IME = zeros(NNZ_IME,1);IME  = zeros(NNZ_IME,1);
%Inverse matrice de masse E (vectorielle)
I_SQRT_IME = zeros(NNZ_IME,1);J_SQRT_IME = zeros(NNZ_IME,1);SQRT_IME  = zeros(NNZ_IME,1);
%Inverse matrice de masse H (scalaire)
I_SQRT_IMH = zeros(NNZ_IMH,1);J_SQRT_IMH = zeros(NNZ_IMH,1);SQRT_IMH  = zeros(NNZ_IMH,1);
%Matrice de rigidité
I_R = zeros(NNZ_R,1);J_R = zeros(NNZ_R,1);R = zeros(NNZ_R,1);
%matrice de penalisation pour l'équation en E (vectorielle)
I_P_E = zeros(NNZ_P_E,1);J_P_E = zeros(NNZ_P_E,1);P_E = zeros(NNZ_P_E,1);
%matrice de penalisation pour l'équation en  H (scalaire) 
I_P_H = zeros(NNZ_P_H,1);J_P_H = zeros(NNZ_P_H,1);P_H = zeros(NNZ_P_H,1);
%Matrice pour tenir compte de la condition de Silver-Muller
I_SM = zeros(NNZ_SM,1);J_SM = zeros(NNZ_SM,1);SM = zeros(NNZ_SM,1);



%Matrice de masse H et son inverse sur le triangle de reference
MELEM_H  = MASSE_ELEM(A,ordre,Nq,xq,yq,wq);
IMELEM_H = MELEM_H\diag(ones(N,1));
%Matrice de masse E et son inverse sur le triangle de reference
MELEM_E  = MASSE_ELEM_V(A,ordre,Nq,xq,yq,wq);
IMELEM_E = MELEM_E\diag(ones(2*N,1));
%Racine carrée des matrices de masse élementaires pour le calcul de la CFL
SQRT_IMELEM_H = IMELEM_H^0.5; 
SQRT_IMELEM_E = IMELEM_E^0.5; 

%Initialisation des autres matrices élementaires
SELEME_aux = zeros(3,2*N,N);
PENE_E_aux = zeros(3,2*N,2*N);
PENE_H_aux = zeros(3,N,N);
SELEME = zeros(3,2*N^2);
PENE_E = zeros(3,4*N^2);
PENE_H = zeros(3,N^2);
ZELEM_SM = zeros(3,N^2);
ZELEM_SM_aux = zeros(3,N,N);

%%%%%%%%%%%%%%%%%%%%
%Boucle d'assemblage
%%%%%%%%%%%%%%%%%%%%
for T = 1:mesh.Nbtri
 xT = mesh.coor(1,mesh.Tri(1:3,T));
 yT = mesh.coor(2,mesh.Tri(1:3,T));
 DF = DFT(xT,yT);
 INV_DFT = transpose(inv(DF));
 JAC = abs(det(DF));
 n1 =mesh.Tri(1,T);
 n2 =mesh.Tri(2,T);
 n3 =mesh.Tri(3,T);
 %%Récupération des informations sur les voisins de T
 voisT = mesh.voisin(:,T);
 na = [n2 n3 n1; n3 n1 n2]; % na(1:2,i)  donne des noeuds de la i eme aerete de T
 xa = [xT(2) xT(3) xT(1); xT(3) xT(1) xT(2)];
 ya = [yT(2) yT(3) yT(1); yT(3) yT(1) yT(2)];
 la1 = norm([xT(3)-xT(2);yT(3)-yT(2)]);
 la2 = norm([xT(1)-xT(3);yT(1)-yT(3)]);
 la3 = norm([xT(2)-xT(1);yT(2)-yT(1)]);
 la = [la1 la2 la3];
 nT = normales_triangle(INV_DFT,nref);
 
 %%Parametres dielectriques
 if mesh.Tri(4,T) == 1
  epsT = eps0;
  muT  = mu0;
  ZT = Z0;
 else
  epsT = eps_diel*eps0;
  muT  = mu_diel*mu0;
  ZT = Z0*sqrt(mu_diel/eps_diel);
 end 

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Calcul des matrice élémentaires
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 RELEM  = RIGIDITE_ELEM(A,ordre,Nq,xq,yq,wq,INV_DFT);
 SELEMI = SAUTI_ELEM(A,ordre,Ng1D,xg1D,wg1D,xT,yT,voisT,xa,ya,la,nT); 
 PENI_E = PENI_E_ELEM(A,ordre,Ng1D,xg1D,wg1D,xT,yT,voisT,xa,ya,la,nT);
 PENI_H = PENI_H_ELEM(A,ordre,Ng1D,xg1D,wg1D,xT,yT,voisT,xa,ya,la,nT);
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%Vectorisation des matrices élémentaires
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 IMELEM_E_V = IMELEM_E(:)/(JAC*epsT);
 IMELEM_H_V = IMELEM_H(:)/(JAC*muT);
 RELEM_V  = RELEM(:)*JAC-SELEMI(:);
 SQRT_IMELEM_H_V = SQRT_IMELEM_H(:)/sqrt(JAC*muT);
 SQRT_IMELEM_E_V = SQRT_IMELEM_E(:)/sqrt(JAC*epsT);
 PENI_E_V = PENI_E(:);
 PENI_H_V = PENI_H(:);
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Interaction de T avec ses voisins
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for i=1:3
  if voisT(i)~=0  %il y a un voisin
   xvT = mesh.coor(1,mesh.Tri(1:3,voisT(i)));
   yvT = mesh.coor(2,mesh.Tri(1:3,voisT(i)));
   SELEME_aux(i,:,:) = SAUTE_ELEM(A,ordre,Ng1D,xg1D,wg1D,xT,yT,xvT,yvT,nT(:,i),xa(:,i),ya(:,i),la(i));
   SELEME(i,:) = SELEME_aux(i,:); %vectorisation
   PENE_E_aux(i,:,:) = PENE_E_ELEM(A,ordre,Ng1D,xg1D,wg1D,xT,yT,xvT,yvT,nT(:,i),xa(:,i),ya(:,i),la(i));
   PENE_E(i,:) = PENE_E_aux(i,:);%vectorisation
   PENE_H_aux(i,:,:) = PENE_H_ELEM(A,ordre,Ng1D,xg1D,wg1D,xT,yT,xvT,yvT,nT(:,i),xa(:,i),ya(:,i),la(i));
   PENE_H(i,:) = PENE_H_aux(i,:);%vectorisation
  else
   %Calcul matrice élementaire pour la condition de Sliver-Muller
   if CL_SM ==1 && Type_arete(na(1,i),na(2,i))==1
   ZELEM_SM_aux(i,:,:) = SM_ELEM(A,ordre,Ng1D,xg1D,wg1D,xT,yT,xa(:,i),ya(:,i),la(i));
   ZELEM_SM(i,:) = ZELEM_SM_aux(i,:);%vectorisation
   end  
  end
 end 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Assemblage des matrices globales
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%Inverse de la matrice de masse en H
I_IMH(kk) = N*(T-1)+ii; 
J_IMH(kk) = N*(T-1)+jj; 
IMH(kk) = IMELEM_H_V(:);

%Inverse de la matrice de masse en E
I_IME(kkv) = 2*N*(T-1)+iiv; 
J_IME(kkv) = 2*N*(T-1)+jjv; 
IME(kkv) = IMELEM_E_V(:);

%Inverse de la racine de la matrice de masse en E
I_SQRT_IME(kkv) = 2*N*(T-1)+iiv; 
J_SQRT_IME(kkv) = 2*N*(T-1)+jjv; 
SQRT_IME(kkv)  = SQRT_IMELEM_E_V(:);


%Inverse de la racine de la matrice de masse en E
I_SQRT_IMH(kk) = N*(T-1)+ii; 
J_SQRT_IMH(kk) =  N*(T-1)+jj; 
SQRT_IMH(kk)  = SQRT_IMELEM_H_V(:);


%Matrice de rigidité : interactions internes au triangle
I_R(kkr) = 2*N*(T-1)+iir; 
J_R(kkr) = N*(T-1)+jjr; 
R(kkr) = RELEM_V(:);
kkr = kkr + 2*N^2;

%Matrice de pénalisation en E : interactions internes au triangle
I_P_E(kkpe) = 2*N*(T-1)+iiv; 
J_P_E(kkpe) = 2*N*(T-1)+jjv;
P_E(kkpe) = PENI_E_V(:);
kkpe = kkpe + 4*N^2;

%Matrice de pénalisation en H : interactions internes au triangle
I_P_H(kkph) = N*(T-1)+ii; 
J_P_H(kkph) = N*(T-1)+jj;
P_H(kkph) = PENI_H_V(:);
kkph = kkph + N^2;

for i=1:3
  if voisT(i)~=0
   %Matrice de rigidité : interactions avec le triangle voisT(i)
   I_R(kkr) = 2*N*(T-1)+iir; 
   J_R(kkr) = N*(voisT(i)-1)+jjr; 
   R(kkr) = SELEME(i,:);
   kkr = kkr + 2*N^2;
   
   %Matrice de pénalisation en E : interactions avec le triangle voisT(i)
   I_P_E(kkpe) = 2*N*(T-1)+iiv; 
   J_P_E(kkpe) = 2*N*(voisT(i)-1)+jjv;
   P_E(kkpe) = -PENE_E(i,:);
   kkpe = kkpe + 4*N^2;

   %Matrice de pénalisation en H : interactions avec le triangle voisT(i)
   I_P_H(kkph) = N*(T-1)+ii; 
   J_P_H(kkph) = N*(voisT(i)-1)+jj;
   P_H(kkph) = -PENE_H(i,:);
   kkph = kkph + N^2;
 else 
 if CL_SM ==1 && Type_arete(na(1,i),na(2,i))==1
   %Matrice de la condition de Silver-Muller
  I_SM(kk_SM) = N*(T-1)+ii; 
  J_SM(kk_SM) = N*(T-1)+jj;
  SM(kk_SM) = -Z0*ZELEM_SM(i,:);
  kk_SM = kk_SM + N^2;
  end 
 end  
end
kk = kk + N^2;
kkv = kkv + 4*N^2;
 
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sparsification des matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IMH = sparse(I_IMH,J_IMH,IMH,Nddl_H,Nddl_H);
IME = sparse(I_IME,J_IME,IME,Nddl_E,Nddl_E);
R = sparse(I_R,J_R,R,Nddl_E,Nddl_H);
P_E = sparse(I_P_E,J_P_E,P_E,Nddl_E,Nddl_E);
P_H = sparse(I_P_H,J_P_H,P_H,Nddl_H,Nddl_H);
if CL_SM ==1
SM = sparse(I_SM,J_SM,SM,Nddl_H,Nddl_H);
end 
SQRT_IMH = sparse(I_SQRT_IMH,J_SQRT_IMH,SQRT_IMH,Nddl_H,Nddl_H);
SQRT_IME = sparse(I_SQRT_IME,J_SQRT_IME,SQRT_IME,Nddl_E,Nddl_E);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Determination de la condition CFL%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[hmin, hmax] = calcul_pas_maillage(mesh);
ACFL = SQRT_IMH*R'*SQRT_IME;
ACFL = ACFL*ACFL';
CFL =2/sqrt(eigs(ACFL,1))*(c0/hmin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Parametre de discrétisation temporelle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dt = CFL/4*hmin^3/c0
dt = CFL*hmin/c0
niter = ceil(Tf/dt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Initialisation de vecteurs%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E0 = zeros(2*N*mesh.Nbtri,1); %Contiendra E_n
E1 = zeros(2*N*mesh.Nbtri,1); %Contiendra E_(n+1)
H0 = zeros(N*mesh.Nbtri,1); %Contiendra H_(n+1/2)
H1 = zeros(N*mesh.Nbtri,1); %Contiendra H_(n+3/2)
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Prise en compte des conditions initiale%%%%%%%
%%%%%%%%%%%%%et de la source%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isource == 0
%Prise en compte des conditions initiales
E0 = interpolation_mode_E(0,mx,my,c0,mesh,ordre,N,x,y);
H0 = interpolation_mode_H(0,mx,my,c0,mesh,ordre,N,x,y);

%Calcul d une approximation d'ordre 2 de H(dt/2) 
dt_H0  = -IMH*(R'*E0+penH*P_H*H0); 
H_half = H0+(dt/2)*dt_H0;
H0     = H_half;
else
if isource == 1
%Prise en compte du point source
RHS =  zeros(2*N*mesh.Nbtri,1);
for T = 1:mesh.Nbtri
 xT = mesh.coor(1,mesh.Tri(1:3,T));
 yT = mesh.coor(2,mesh.Tri(1:3,T));
 DF = DFT(xT,yT);
 JAC = abs(det(DF));
 INV_DFT = transpose(inv(DF)); 
 JELEM = RHS_ELEM(xs,ys,rs,A,ordre,Nq,xq,yq,wq,xT,yT);
 JELEM = JAC*JELEM;
 %Assemblage
  for k = 1:N
   Ix = 2*(T-1)*N + 2*(k-1)+1;  %%numero global du kieme ddl Ex du triangle T
   Iy = 2*(T-1)*N + 2*(k-1)+2;  %%numero global du kieme ddl Ey du triangle T
   RHS(Ix,1) = RHS(Ix,1) + JELEM(2*(k-1)+1,1);
   RHS(Iy,1) = RHS(Iy,1) + JELEM(2*(k-1)+2,1);
  end     
end 
end
end


%%%%%%%%%%%%%%%%%%%%%%%
%%%%Itérations en temps
%%%%%%%%%%%%%%%%%%%%%%%

if flag_film == 1
figure;
 Film(1) = visu_solution_Film_E(E0,1,A,ordre,mesh,mesh_ref,Xmin,Xmax,Ymin,Ymax,Zmin,Zmax);
nf=1;
end 


for n=1:niter
 if(mod(n,100) ==0) 
  disp('nombre d iterations')
  n
 end     
 tn = (n-1)*dt;

 %Avancée en E
 Vaux1(:,1) = R*H0-penE*P_E*E0;
 if isource == 1
  fct = exp(-(tn+dt/2-t0)^2/sigma_t^2);
  Vaux1(:,1) = Vaux1(:,1) + fct*RHS;
 end
 E1(:,1) = E0(:,1) + dt*IME*Vaux1(:,1);

 %Avancée en H
 Vaux2(:,1) = R'*E1(:,1)+penH*P_H*H0;
 if CL_SM ==1 
  Vaux2(:,1) =Vaux2(:,1)-SM*H0;
 end 
 H1(:,1) = H0(:,1) - dt*IMH*Vaux2(:,1);
 
 %Mise à jour
 E0 = E1;
 H0 = H1;
 
 if flag_film == 1 &&  mod(n,20) == 0 
  nf =nf +1;
 Film(nf) = visu_solution_Film_E(E1,1,A,ordre,mesh,mesh_ref,Xmin,Xmax,Ymin,Ymax,Zmin,Zmax);
 end 
 
end 


[erreurL2_E,erreurL2_H]  = erreurL2_mode(E0,H0,mx,my,c0,niter,dt,A,ordre,mesh)
hmax
visu_solution(E0,H0,A,ordre,mesh,mesh_ref);
if flag_film == 1
 movie(Film);
end

var = sqrt(mesh.Nbtri*N)