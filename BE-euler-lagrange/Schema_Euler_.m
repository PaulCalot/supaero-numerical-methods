function [d,v,clN1,clS2,clE3,clO4] = Schema_Euler(Nx,Ny,Lx,Ly,din,vin,cln,cls,cle,clo,dt,Ug,tau)

d = din;
v = vin;

dx = Lx/Nx;
dy = Ly/Ny;

%% Nord-Sud-Est-Ouest <-> 1-2-3-4
clN1 = zeros(Nx,3); 
clS2 = zeros(Nx,3);
clE3 = zeros(Ny,3);
clO4 = zeros(Ny,3);

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
  if(FluxvH(Nx+1,i,1) > 0)
    clE3(i,1) = FluxdH(Nx+1,i)/dy;
    clE3(i,2:3) = FluxvH(Nx+1,i,:)/dy;
  end;
  if(FluxvH(1,i,:) < 0)
    clO4(i,1) = FluxdH(1,i)/dy;
    clO4(i,2:3) = FluxvH(1,i,:)/dy;
  end;

  %% A COMPLETER : CL Entrantes pour les flux entrants
  FluxdH(1,i) =  FluxdH(1,i) + dy * max(0, clo(i, 1));
  FluxdH(Nx+1,i) = FluxdH(Nx+1,i) + dy * min(0, cle(i, 1));
  FluxvH(1,i,1) =   FluxvH(1,i,1) + dy * max(0, clo(i, 2));
  FluxvH(1,i,2) =   FluxvH(1,i,2) + dy * max(0, clo(i, 3));
  FluxvH(Nx+1,i,1) = FluxvH(Nx+1,i,1) + dy * min(0, cle(i, 2));
  FluxvH(Nx+1,i,2) = FluxvH(Nx+1,i,2) + dy * min(0, cle(i, 3));
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
  if(FluxvV(i, Ny+1, 2) > 0)
      clN1(i,1) = max(0, FluxdV(i, Ny+1)/dx);
      clN1(i,2:3) = FluxvV(i, Ny+1)/dx;
  end;
  if(FluxvV(i,1,2) < 0)
    clS2(i,1) = FluxdV(i, 1)/dx;
    clS2(i,2:3) = FluxvV(i,1,:)/dx;
  end;
  %% A COMPLETER : CL Entrantes pour les flux entrants
  FluxdV(i,1) = FluxdV(i, 1) + dy * max(0, cls(i, 1));
  FluxdV(i,Ny+1) = FluxdV(i, Ny+1) + dy * min(0, cln(i, 1));
  FluxvV(i,1,1) =  FluxvV(i,1,1) + dy * max(0, cls(i, 2));
  FluxvV(i,1,2) = FluxvV(i,1,2) + dy * max(0, cls(i, 3));
  FluxvV(i,Ny+1,1) = FluxvV(i,Ny+1,1) + dy * min(0, cln(i, 2));
  FluxvV(i,Ny+1,2) = FluxvV(i,Ny+1,2) + dy * min(0, cln(i, 3));
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
  qmvt(i,j,1) = qmvt(i,j,1) + dt/tau * (din(i,j) * Ug(i,j,1) - qmvt(i,j,1));
  qmvt(i,j,2) = qmvt(i,j,2) + dt/tau * (din(i,j) * Ug(i,j,2) - qmvt(i,j,2));

%% Reconstruction de la vitesse a partir de la quantite de mouvement(pour affichage)
  if (d(i,j) > 0.00000001) 
    v(i,j,1) = qmvt(i,j,1)/d(i,j);
    v(i,j,2) = qmvt(i,j,2)/d(i,j);
  else
    v(i,j,:) = 0.;
  end;
end;
end;
