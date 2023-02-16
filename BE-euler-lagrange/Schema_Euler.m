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
end
end

%% A COMPLETER
%% Integrations des CL Horizontales (Est et Ouest) dans les flux 
for i=1:Ny
  %% A COMPLETER : Flux sortants pour les CL Sortantes
  cl3(i,1) = FluxdH(Nx+1,i)/dy;
  cl3(i,2:3) = FluxvH(Nx+1,i,:)/dy;
  cl4(i,1) = FluxdH(1,i)/dy;
  cl4(i,2:3) = FluxvH(1,i,:)/dy;

  %% A COMPLETER : CL Entrantes pour les flux entrants
  FluxdH(1,i) =  FluxdH(1,i) + clo(i,1)*dy;
  FluxdH(Nx+1,i) = FluxdH(Nx+1,i) + cle(i,1)*dy;
  FluxvH(1,i,1) =   FluxvH(1,i,1) + clo(i,2)*dy;
  FluxvH(1,i,2) =   FluxvH(1,i,2) + clo(i,3)*dy;
  FluxvH(Nx+1,i,1) = FluxvH(Nx+1,i,1) + cle(i,2)*dy;
  FluxvH(Nx+1,i,2) =   FluxvH(Nx+1,i,2) + cle(i,3)*dy;
end

%% Calcul des flux verticaux
for i=1:Nx
for j=1:Ny
  FluxdV(i,j) = FluxdV(i,j) + dx*min(0.,vin(i,j,2))*din(i,j);
  FluxdV(i,j+1) = FluxdV(i,j+1) + dx*max(0.,vin(i,j,2))*din(i,j);
  FluxvV(i,j,1) = FluxvV(i,j,1) + dx*min(0.,vin(i,j,2))*din(i,j)*vin(i,j,1);
  FluxvV(i,j,2) = FluxvV(i,j,2) + dx*min(0.,vin(i,j,2))*din(i,j)*vin(i,j,2);
  FluxvV(i,j+1,1) = FluxvV(i,j+1,1) + dx*max(0.,vin(i,j,2))*din(i,j)*vin(i,j,1);
  FluxvV(i,j+1,2) = FluxvV(i,j+1,2) + dx*max(0.,vin(i,j,2))*din(i,j)*vin(i,j,2);
end
end

%% Integrations des CL Verticales (Nord et Sud) dans les flux
for i=1:Nx
  %% A COMPLETER : Flux sortants pour les CL Sortantes
  cl1(i,1) = FluxdV(i,Ny+1)/dx;
  cl1(i,2:3) = FluxvV(i,Ny+1,:)/dx ;
  cl2(i,1) = FluxdV(i,1)/dx ;
  cl2(i,2:3) = FluxvV(i,1,:)/dx ;

  %% A COMPLETER : CL Entrantes pour les flux entrants
  FluxdV(i,1) = FluxdV(i,1) + cls(i,1)*dx;
  FluxdV(i,Ny+1) = FluxdV(i,Ny+1) + cln(i,1)*dx;
  FluxvV(i,1,1) =  FluxvV(i,1,1) + cls(i,2)*dx;
  FluxvV(i,1,2) = FluxvV(i,1,2) + cls(i,3)*dx;
  FluxvV(i,Ny+1,1) = FluxvV(i,Ny+1,1) + cln(i,2)*dx;
  FluxvV(i,Ny+1,2) = FluxvV(i,Ny+1,2) + cln(i,3)*dx;
end

%% Mise a jour des densites et vitesses
for i=1:Nx
for j=1:Ny
  qmvt(i,j,1) = din(i,j)*vin(i,j,1);
  qmvt(i,j,2) = din(i,j)*vin(i,j,2);
  d(i,j) = din(i,j) - dt/dx/dy*(FluxdH(i+1,j)+FluxdV(i,j+1)-FluxdH(i,j)-FluxdV(i,j));
  qmvt(i,j,1) = qmvt(i,j,1) - dt/dx/dy*(FluxvH(i+1,j,1)+FluxvV(i,j+1,1)-FluxvH(i,j,1)-FluxvV(i,j,1));
  qmvt(i,j,2) = qmvt(i,j,2) - dt/dx/dy*(FluxvH(i+1,j,2)+FluxvV(i,j+1,2)-FluxvH(i,j,2)-FluxvV(i,j,2));

  %% A COMPLETER : Integration temporelle de la quantite de mouvement avec le terme de rappel vers le gaz
  qmvt(i,j,:) = qmvt(i,j,:)+dt/tau*(din(i,j)*Ug(i,j,:)-qmvt(i,j,:));

%% Reconstruction de la vitesse a partir de la quantite de mouvement(pour affichage)
  if (d(i,j) > 0.00000001) 
    v(i,j,1) = qmvt(i,j,1)/d(i,j);
    v(i,j,2) = qmvt(i,j,2)/d(i,j);
  else
    v(i,j,:) = 0.;
  end
end
end
