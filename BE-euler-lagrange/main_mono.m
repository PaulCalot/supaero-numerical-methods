%% Programme principal MONODOMAINE

Nx = 50;
Ny = 50;

Xmin = 0.0;
Ymin = 0.0;

LongX = 1.0;
LongY = 1.0;

Dx = LongX/Nx;
Dy = LongY/Ny;

%% Type de schema utilise : 0 = Euler, 1 = Lagrange
ModSchem = 0;

%% Valeur du temps caracteristique de trainee
tau = 1;

%% Initialisation
Density = zeros(Nx,Ny);
Vitesse = zeros(Nx,Ny,2);
flagr = [];

Ugaz = ones(Nx,Ny,2);

%% Pas de temps pour la simulation
dt = min(tau/10,min(Dx,Dy)/5.);


%% Conditions aux limites (Nord - Sud - Est - Ouest) - Initialisation des tableaux
CLN = zeros(Nx,3);
CLS = zeros(Nx,3);
CLE = zeros(Ny,3);
CLO = zeros(Ny,3);

%% A COMPLETER : Conditions aux limites
fd = 1;
fqdm = 1; 
CLS(floor(Nx/2)-2:floor(Nx/2)+2,1) = fd; %% Flux de densite en kg.m^(-1).s^(-1)
CLS(floor(Nx/2)-2:floor(Nx/2)+2,2) = 0; %% Flux de quantite de mouvement (en x) en kg.s^(-2)
CLS(floor(Nx/2)-2:floor(Nx/2)+2,3) = fqdm; %% Flux de quantite de mouvement (en y) en kg.s^(-2)
CLO(floor(Ny/2)-2:floor(Ny/2)+2,1) = fd;
CLO(floor(Ny/2)-2:floor(Ny/2)+2,2) = fqdm;
CLO(floor(Ny/2)-2:floor(Ny/2)+2,3) = 0;


%% Utilisation du schema Euler Explicite en temps
figure(5); quiver(Ugaz(:,:,1)',Ugaz(:,:,2)');

for i=1:20
 figure(1); contourf(Density', 10); %'
 xl1 = xlabel('npx');
 yl1 = ylabel('npy');
 colorbar;
 figure(2); contourf(Vitesse(:,:,1)', 10); %'
 xl2 = xlabel('npx');
 yl2 = ylabel('npy');
 colorbar;
 figure(3); contourf(Vitesse(:,:,2)', 10); %'
 xl3 = xlabel('npx');
 yl3 = ylabel('npy');
 colorbar;
 figure(4); quiver(Vitesse(:,:,1)',Vitesse(:,:,2)');
 xl4 = xlabel('npx');
 yl4 = ylabel('npy');
 drawnow;
 Npart = size(flagr,1);
 %display('Appuyez sur Entree pour continuer ...');
 %pause; 

 for j=1:50

   if (ModSchem > 0)
     [d,v,floc,cl11,cl21,cl31,cl41] = Schema_Lagrange(Nx,Ny,LongX,LongY,flagr,CLN,CLS,CLE,CLO,dt,Ugaz,tau);
     flagr = floc;
   else
     [d,v,cl11,cl21,cl31,cl41] = Schema_Euler(Nx,Ny,LongX,LongY,Density,Vitesse,CLN,CLS,CLE,CLO,dt,Ugaz,tau);
   end;
   Density = d;
   Vitesse = v;

 end;
end;

axis([0 50 0 50]);
grid minor
saveas(figure(1), "results/mono_density_kind:" + num2str(ModSchem) + "_tau:" + num2str(tau) + ".png");
saveas(figure(2), "results/mono_vx_kind:" + num2str(ModSchem) + "_tau:" + num2str(tau) + ".png");
saveas(figure(3), "results/mono_vy_kind:" + num2str(ModSchem) + "_tau:" + num2str(tau) + ".png");
saveas(figure(4), "results/mono_velocity_quiver_kind:" + num2str(ModSchem) + "_tau:" + num2str(tau) + ".png");