function [d,v,cl1,cl2,cl3,cl4] = Schema_Euler(Nx,Ny,Lx,Ly,din,vin,cln,cls,cle,clo,dt,Ug,tau)

d = din;
v = vin;

dx = Lx/Nx;
dy = Ly/Ny;

%% Nord-Sud-Est-Ouest <-> 1-2-3-4
cl1 = zeros(Nx,3);
cl2 = zeros(Nx,3);
cl3 = zeros(Ny,3);
cl4 = zeros(Ny,3);

FluxdV = zeros(Nx,Ny+1);
FluxdH = zeros(Nx+1,Ny);
FluxvV = zeros(Nx,Ny+1,2);
FluxvH = zeros(Nx+1,Ny,2);

%% Calcul des flux horizontaux
for i=1:Ny
for j=1:Nx
  FluxdH(j,i) = FluxdH(j,i) + dy*min(0.,vin(j,i,1))*din(j,i);
  FluxdH(j+1,i) = FluxdH(j+1,i) + dy*max(0.,vin(j,i,1))*din(j,i);
  FluxvH(j,i,1) = FluxvH(j,i,1) + dy*min(0.,vin(j,i,1))*din(j,i)*vin(j,i,1);
  FluxvH(j,i,2) = FluxvH(j,i,2) + dy*min(0.,vin(j,i,1))*din(j,i)*vin(j,i,2);
  FluxvH(j+1,i,1) = FluxvH(j+1,i,1) + dy*max(0.,vin(j,i,1))*din(j,i)*vin(j,i,1);
  FluxvH(j+1,i,2) = FluxvH(j+1,i,2) + dy*max(0.,vin(j,i,1))*din(j,i)*vin(j,i,2);
end;
end;

%% A COMPLETER
%% Integrations des CL Horizontales (Est et Ouest) dans les flux 
for i=1:Ny
  %% A COMPLETER : Flux sortants pour les CL Sortantes
  cl3(i,1) = FluxdH(Nx+1,i)/...;
  cl3(i,2:3) = FluxvH(Nx+1,i,:)/...;
  cl4(i,1) = FluxdH(1,i)/...;
  cl4(i,2:3) = FluxvH(1,i,:)/...;

  %% A COMPLETER : CL Entrantes pour les flux entrants
  FluxdH(1,i) =  FluxdH(1,i) + ...
  FluxdH(Nx+1,i) = FluxdH(Nx+1,i) + ...
  FluxvH(1,i,1) =   FluxvH(1,i,1) + ...
  FluxvH(1,i,2) =   FluxvH(1,i,2) + ...
  FluxvH(Nx+1,i,1) = FluxvH(Nx+1,i,1) + ...
  FluxvH(Nx+1,i,2) =   FluxvH(Nx+1,i,2) + ...
end;

%% Calcul des flux verticaux
for i=1:Nx
for j=1:Ny
  FluxdV(i,j) = FluxdV(i,j) + dx*min(0.,vin(i,j,2))*din(i,j);
  FluxdV(i,j+1) = FluxdV(i,j+1) + dx*max(0.,vin(i,j,2))*din(i,j);
  FluxvV(i,j,1) = FluxvV(i,j,1) + dx*min(0.,vin(i,j,2))*din(i,j)*vin(i,j,1);
  FluxvV(i,j,2) = FluxvV(i,j,2) + dx*min(0.,vin(i,j,2))*din(i,j)*vin(i,j,2);
  FluxvV(i,j+1,1) = FluxvV(i,j+1,1) + dx*max(0.,vin(i,j,2))*din(i,j)*vin(i,j,1);
  FluxvV(i,j+1,2) = FluxvV(i,j+1,2) + dx*max(0.,vin(i,j,2))*din(i,j)*vin(i,j,2);
end;
end;

%% Integrations des CL Verticales (Nord et Sud) dans les flux
for i=1:Nx
  %% A COMPLETER : Flux sortants pour les CL Sortantes
  cl1(i,1) = ... ;
  cl1(i,2:3) = ... ;
  cl2(i,1) = ... ;
  cl2(i,2:3) = ... ;

  %% A COMPLETER : CL Entrantes pour les flux entrants
  FluxdV(i,1) = ..
  FluxdV(i,Ny+1) = ..
  FluxvV(i,1,1) = ..
  FluxvV(i,1,2) = ..
  FluxvV(i,Ny+1,1) = ..
  FluxvV(i,Ny+1,2) = ..
end;

%% Mise a jour des densites et vitesses
for i=1:Nx
for j=1:Ny
  qmvt(i,j,1) = din(i,j)*vin(i,j,1);
  qmvt(i,j,2) = din(i,j)*vin(i,j,2);
  d(i,j) = din(i,j) - dt/dx/dy*(FluxdH(i+1,j)+FluxdV(i,j+1)-FluxdH(i,j)-FluxdV(i,j));
  qmvt(i,j,1) = qmvt(i,j,1) - dt/dx/dy*(FluxvH(i+1,j,1)+FluxvV(i,j+1,1)-FluxvH(i,j,1)-FluxvV(i,j,1));
  qmvt(i,j,2) = qmvt(i,j,2) - dt/dx/dy*(FluxvH(i+1,j,2)+FluxvV(i,j+1,2)-FluxvH(i,j,2)-FluxvV(i,j,2));

  %% A COMPLETER : Integration temporelle de la quantite de mouvement avec le terme de rappel vers le gaz
  qmvt(i,j,:) = ..

%% Reconstruction de la vitesse a partir de la quantite de mouvement(pour affichage)
  if (d(i,j) > 0.00000001) 
    v(i,j,1) = qmvt(i,j,1)/d(i,j);
    v(i,j,2) = qmvt(i,j,2)/d(i,j);
  else
    v(i,j,:) = 0.;
  end;
end;
end;
