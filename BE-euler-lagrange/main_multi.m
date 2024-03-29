%% Programme principal MULTIDOMAINES

Nxint = [20,10,20];
Nyint = [20,10,20];

Nx = sum(Nxint);
Ny = sum(Nyint);

Xmin = 0.0;
Ymin = 0.0;

LongX = 1.0;
LongY = 1.0;

Dx = LongX/Nx;
Dy = LongY/Ny;

%% Temps caracteristique pour la trainee
tau = 1;

%% Choix du schema utilise par domaine 0 = Eulerien - 1 = Lagrangien
%% 2 = multi-domain (pre-set values)
ModSchem_ = 5;
ModSchem = zeros(3);
if(ModSchem_ == 1)
    ModSchem(1, 1)= 1;
    ModSchem(1, 2)= 1;
    ModSchem(1, 3)= 1;
    ModSchem(2, 1)= 1;
    ModSchem(2, 2)= 1;
    ModSchem(2, 3)= 1;
    ModSchem(3, 1)= 1;
    ModSchem(3, 2)= 1;
    ModSchem(3, 3)= 1;
elseif(ModSchem_ == 2)
    ModSchem(1, 1)= 0;
    ModSchem(1, 2)= 0;
    ModSchem(1, 3)= 0;
    ModSchem(2, 1)= 0;
    ModSchem(2, 2)= 1;
    ModSchem(2, 3)= 0;
    ModSchem(3, 1)= 0;
    ModSchem(3, 2)= 0;
    ModSchem(3, 3)= 0;
elseif(ModSchem_ == 3)
    ModSchem(1, 1)= 0;
    ModSchem(1, 2)= 0;
    ModSchem(1, 3)= 0;
    ModSchem(2, 1)= 0;
    ModSchem(2, 2)= 1;
    ModSchem(2, 3)= 1;
    ModSchem(3, 1)= 0;
    ModSchem(3, 2)= 1;
    ModSchem(3, 3)= 1;
elseif(ModSchem_ == 4)
    ModSchem(1, 1)= 0;
    ModSchem(1, 2)= 0;
    ModSchem(1, 3)= 1;
    ModSchem(2, 1)= 0;
    ModSchem(2, 2)= 1;
    ModSchem(2, 3)= 1;
    ModSchem(3, 1)= 1;
    ModSchem(3, 2)= 1;
    ModSchem(3, 3)= 1;
elseif(ModSchem_ == 5)
    ModSchem(1, 1)= 0;
    ModSchem(1, 2)= 1;
    ModSchem(1, 3)= 0;
    ModSchem(2, 1)= 1;
    ModSchem(2, 2)= 1;
    ModSchem(2, 3)= 1;
    ModSchem(3, 1)= 0;
    ModSchem(3, 2)= 1;
    ModSchem(3, 3)= 1;
end

%% Initialisation des variables Euleriennes
Density = zeros(Nx,Ny);
Vitesse = zeros(Nx,Ny,2);
%% Initialisation des variables Lagrangiiennes
flagr = cell(3,3);
for i=1:3
  for j=1:3
    flagr{i,j}=[];
  end;
end;

Ugaz = ones(Nx,Ny,2);

%% Pas de temps pour la simulation
dt = min(tau/10,min(Dx,Dy)/5.);

%% CONDITIONS AUX LIMITES (Nord - Sud - Est - Ouest) POUR TOUT LE DOMAINE
CLN = zeros(Nx,3);
CLS = zeros(Nx,3);
CLE = zeros(Ny,3);
CLO = zeros(Ny,3);

%% A COMPLETER : Conditions aux limites
fd = 1;
fqdm = 1; 
CLS(floor(Nx/2)-2:floor(Nx/2)+2,1) = fd; %% Flux de densite en kg.m^(-1).s^(-1) -> flux linéique
CLS(floor(Nx/2)-2:floor(Nx/2)+2,2) = 0; %% Flux de quantite de mouvement (en x) en kg.s^(-2)
CLS(floor(Nx/2)-2:floor(Nx/2)+2,3) = fqdm; %% Flux de quantite de mouvement (en y) en kg.s^(-2)
CLO(floor(Ny/2)-2:floor(Ny/2)+2,1) = fd;
CLO(floor(Ny/2)-2:floor(Ny/2)+2,2) = fqdm;
CLO(floor(Ny/2)-2:floor(Ny/2)+2,3) = 0;

%% Utilisation du schema Euler Explicite en temps

%%% Initialisation des CL PAR DOMAINE
clin = cell(3);
clis = cell(3);
clie = cell(3);
clio = cell(3);

clin{1,1} = zeros(Nxint(1),3);
clin{1,2} = zeros(Nxint(1),3);
clin{1,3} = CLN(1:Nxint(1),:);
clin{2,1} = zeros(Nxint(2),3);
clin{2,2} = zeros(Nxint(2),3);
clin{2,3} = CLN(Nxint(1)+1:Nxint(1)+Nxint(2),:);
clin{3,1} = zeros(Nxint(3),3);
clin{3,2} = zeros(Nxint(3),3);
clin{3,3} = CLN(Nxint(1)+Nxint(2)+1:Nxint(1)+Nxint(2)+Nxint(3),:);

clis{1,1} = CLS(1:Nxint(1),:);
clis{1,2} = zeros(Nxint(1),3);
clis{1,3} = zeros(Nxint(1),3);
clis{2,1} = CLS(Nxint(1)+1:Nxint(1)+Nxint(2),:);
clis{2,2} = zeros(Nxint(2),3);
clis{2,3} = zeros(Nxint(2),3);
clis{3,1} = CLS(Nxint(1)+Nxint(2)+1:Nxint(1)+Nxint(2)+Nxint(3),:);
clis{3,2} = zeros(Nxint(3),3);
clis{3,3} = zeros(Nxint(3),3);

clie{1,1} = zeros(Nyint(1),3); 
clie{1,2} = zeros(Nyint(2),3);
clie{1,3} = zeros(Nyint(3),3);
clie{2,1} = zeros(Nyint(1),3); 
clie{2,2} = zeros(Nyint(2),3);
clie{2,3} = zeros(Nyint(3),3);
clie{3,1} = CLE(1:Nyint(1),:);
clie{3,2} = CLE(Nyint(1)+1:Nyint(1)+Nyint(2),:);
clie{3,3} = CLE(Nyint(1)+Nyint(2)+1:Nyint(1)+Nyint(2)+Nyint(3),:);

clio{1,1} = CLO(1:Nyint(1),:);
clio{1,2} = CLO(Nyint(1)+1:Nyint(1)+Nyint(2),:);
clio{1,3} = CLO(Nyint(1)+Nyint(2)+1:Nyint(1)+Nyint(2)+Nyint(3),:);
clio{2,1} = zeros(Nyint(1),3); 
clio{2,2} = zeros(Nyint(2),3);
clio{2,3} = zeros(Nyint(3),3);
clio{3,1} = zeros(Nyint(1),3); 
clio{3,2} = zeros(Nyint(2),3);
clio{3,3} = zeros(Nyint(3),3);


%% Affichage du champ de vitesse gazeux
figure(5); quiver(Ugaz(:,:,1)',Ugaz(:,:,2)');


%% Boucle en temps
for i=1:20
  figure(1); contourf(Density', 20); % [0. 0.1 0.2 0.3 0.4 0.5 0.6 0.7 1 1.2 1.5 2 2.5 3 3.5 4 5 6]
  xl1 = xlabel('npx');
  yl1 = ylabel('npy');
  colorbar;
  figure(2); contourf(Vitesse(:,:,1)', 20); %'
  xl2 = xlabel('npx');
  yl2 = ylabel('npy');
  colorbar;
  figure(3); contourf(Vitesse(:,:,2)', 20); %'
  xl3 = xlabel('npx');
  yl3 = ylabel('npy');
  colorbar;
  figure(4); quiver(Vitesse(:,:,1)',Vitesse(    :,:,2)');
  xl4 = xlabel('npx');
  yl4 = ylabel('npy');
  drawnow;
  Npart = 0;
  for k=1:3
    for l=1:3
      Npart = Npart + size(flagr{k,l},1);
    end;
  end;
  Npart
  %display('Appuyez sur Entree pour continuer ...');
  %pause; 

  for j=1:50

    for k=1:3
    for l=1:3

    NxLoc = Nxint(k);
    NyLoc = Nyint(l);
    LxLoc = LongX*NxLoc/Nx;
    LyLoc = LongY*NyLoc/Ny;
    DebNx = sum(Nxint(1:k-1));
    DebNy = sum(Nyint(1:l-1));
    DensLoc = Density(DebNx+1:DebNx+NxLoc, DebNy+1:DebNy+NyLoc);
    VitLoc = Vitesse(DebNx+1:DebNx+NxLoc, DebNy+1:DebNy+NyLoc,:);
    UgLoc = Ugaz(DebNx+1:DebNx+NxLoc, DebNy+1:DebNy+NyLoc,:);

    if (ModSchem(k,l) > 0)
      [d,v,floc,cl11,cl21,cl31,cl41] = Schema_Lagrange(NxLoc,NyLoc,LxLoc,LyLoc,flagr{k,l},clin{k,l},clis{k,l},clie{k,l},clio{k,l},dt, UgLoc, tau);
      flagr{k,l} = floc;
    else
      [d,v,cl11,cl21,cl31,cl41] = Schema_Euler(NxLoc,NyLoc,LxLoc,LyLoc,DensLoc,VitLoc,clin{k,l},clis{k,l},clie{k,l},clio{k,l},dt, UgLoc, tau);
    end;
    Density(DebNx+1:DebNx+NxLoc, DebNy+1:DebNy+NyLoc) = d;
    Vitesse(DebNx+1:DebNx+NxLoc, DebNy+1:DebNy+NyLoc,:) = v;

    %% A COMPLETER
    %% Redistribution des Conditions Limites (ce qui sort d'un domaine rentrera dans le domaine voisin)
    %% Est -> Ouest
    if (k<3)
     clio{k+1,l} = cl31; % clio{k, l};
    end
    %% Ouest -> Est
    if (k>1) 
      clie{k-1,l} = cl41; % clie{k, l};
    end;
    %% Nord -> Sud
    if (l<3)
      clis{k,l+1} = cl11; % clis{k, l};
    end;
    %% Sud -> Nord
    if (l>1) 
      clin{k,l-1} = cl21; %clin{k, l};
    end;
    end;
    end;

  end;
end;
axis([0 50 0 50]);
grid minor
saveas(figure(1), "results/multi_density_kind:" + num2str(ModSchem_) + "_tau:" + num2str(tau) + ".png");
saveas(figure(2), "results/multi_vx_kind:" + num2str(ModSchem_) + "_tau:" + num2str(tau) + ".png");
saveas(figure(3), "results/multi_vy_kind:" + num2str(ModSchem_) + "_tau:" + num2str(tau) + ".png");
saveas(figure(4), "results/multi_velocity_quiver_kind:" + num2str(ModSchem_) + "_tau:" + num2str(tau) + ".png");


