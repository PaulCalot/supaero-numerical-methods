function [d,v,fl,cl11,cl21,cl31,cl41] = Schema_Lagrange(Nx,Ny,Lx,Ly,fin,cln,cls,cle,clo,dt,Ug,tau)

Dx = Lx/Nx;
Dy = Ly/Ny;

%% Nombre total de particules numeriques
Np = size(fin,1);
fl = fin;

%% Nombre de particules numeriques injectees par face du maillage pour representer les CL
kcl = 1;

%% I - On introduit les nouvelles particules venant des CL
%% A COMPLETER 
%% Nord
for i=1:Nx
  if (cln(i,1) < 0) %% On vérifie que le flux est entrant dans le domaine 
    mass = - cln(i, 1) * Dx * dt / kcl;
    qdmx = cln(i, 2) * Dx * dt / kcl; % on devrait avoir /kcl ici également je pense
    qdmy = cln(i, 3) * Dx * dt / kcl;
    for k=1:kcl
      fl(Np+k,1) = mass;    %% Masse de la particule numerique injectee
      fl(Np+k,2) = qdmx / mass; %% Vitesse en x de la particule numerique injectee
      fl(Np+k,3) = qdmy / mass; %% Vitesse en y de la particule numerique injectee
      fl(Np+k,4) = (i - 1 + k/(kcl + 1)) * Dx; %% Abscisse de la particule numerique injectee
      fl(Np+k,5) = Ny * Dy; %% Ordonnee de la particule numerique injectee
    end;
    Np = Np+kcl;
  end;
end;
%% Sud
for i=1:Nx
  if (cls(i,1) > 0) %% On vérifie que le flux est entrant dans le domaine
    mass = cls(i, 1) * Dx * dt / kcl;
    qdmx = cls(i, 2) * Dx * dt / kcl;
    qdmy = cls(i, 3) * Dx * dt / kcl;
    for k=1:kcl 
      fl(Np+k,1) = mass;
      fl(Np+k,2) = qdmx / mass; 
      fl(Np+k,3) = qdmy / mass;
      fl(Np+k,4) = (i - 1 + k/(kcl + 1)) * Dx;
      fl(Np+k,5) = 0;
    end;
    Np = Np+kcl;
  end;
end;
%% Est
for i=1:Ny
  if (cle(i,1) < 0) %% On vérifie que le flux est entrant dans le domaine 
    mass = - cle(i, 1) * Dy * dt / kcl;
    qdmx = cle(i, 2) * Dy * dt / kcl;
    qdmy = cle(i, 3) * Dy * dt / kcl;
    for k=1:kcl
      fl(Np+k,1) = mass;
      fl(Np+k,2) = qdmx / mass; 
      fl(Np+k,3) = qdmy / mass;
      fl(Np+k,4) = (Nx) * Dx;
      fl(Np+k,5) = (i - 1 + k/(kcl + 1)) * Dy;
    end;
    Np = Np+kcl;
  end;
end;
%% Ouest
for i=1:Ny
  if (clo(i,1) > 0) %% On vérifie que le flux est entrant dans le domaine 
    mass = clo(i, 1) * Dy * dt / kcl;
    qdmx = clo(i, 2) * Dy * dt / kcl;
    qdmy = clo(i, 3) * Dy * dt / kcl;
    for k=1:kcl
      fl(Np+k,1) = mass;
      fl(Np+k,2) = qdmx / mass;
      fl(Np+k,3) = qdmy / mass;
      fl(Np+k,4) = 0;
      fl(Np+k,5) = (i - 1 + k/(kcl + 1)) * Dy;
    end;  
    Np= Np+kcl;
  end;
end;

%% II - On deplace les particules pendant dt et on ajuste les effets de trainee

%% A COMPLETER
for i=1:Np
  %% Position de la particule dans le maillage
  npx = max(1, min(Nx, floor(fl(i,4)/Dx)+1 ));
  npy = max(1, min(Ny, floor(fl(i,5)/Dy)+1 ));

  %% Effets de trainee
  fl(i,2) = fl(i,2) + dt/tau * (Ug(npx, npy, 1) - fl(i, 2));
  fl(i,3) = fl(i,3) + dt/tau * (Ug(npx, npy, 2) - fl(i, 3));

  %% Transport libre
  fl(i,4) = fl(i,4) + dt * fl(i, 2); 
  fl(i,5) = fl(i,5) + dt * fl(i, 3);
end;


%% III - On localise les particules et on calcule d et v et les cl
%%       On ne garde que les particules dans le domaine
%%       On calcule les flux sortants du domaine (CLs de sortie)
LocPart = cell(Nx,Ny);
Npc = zeros(Nx,Ny);
d = zeros(Nx,Ny);
v = zeros(Nx,Ny,2);

%% Nord-Sud-Est-Ouest <-> 11-21-31-41
cl11 = zeros(Nx,3);
cl21 = zeros(Nx,3);
cl31 = zeros(Ny,3);
cl41 = zeros(Ny,3);

flt = [];
Ntmp = 0;
for i=1:Np
  npx = floor(fl(i,4)/Dx)+1; %% Localisation de la particule numerique
  npy = floor(fl(i,5)/Dy)+1;
  %% A COMPLETER
  if (npx > Nx)
    ntpy = min(Ny,max(1,npy));
    if(fl(i,2) > 0) % vx > 0
        cl31(ntpy,1) = cl31(ntpy,1) + fl(i,1) / (Dy * dt); %% Flux de densite en kg.m^(-1).s^(-1)
        cl31(ntpy,2) = cl31(ntpy,2) + fl(i,2) * fl(i, 1) / (Dy * dt); %% Flux de quantite de mouvement (en x) en kg.s^(-2)
        cl31(ntpy,3) = cl31(ntpy,3) + fl(i,3) * fl(i, 1) / (Dy * dt); %% Flux de quantite de mouvement (en y) en kg.s^(-2)
    end
  elseif (npx < 1)
    ntpy = min(Ny,max(1,npy));
    if(fl(i,2) < 0)
        cl41(ntpy,1) = cl41(ntpy,1) - fl(i,1) / (Dy * dt); 
        cl41(ntpy,2) = cl41(ntpy,2) + fl(i,2) * fl(i,1) / (Dy * dt); 
        cl41(ntpy,3) = cl41(ntpy,3) + fl(i,3) * fl(i,1) / (Dy * dt); 
    end
  elseif (npy > Ny)
    ntpx = min(Nx,max(1,npx));
    if(fl(i, 3) > 0)
        cl11(ntpx,1) = cl11(ntpx,1) + fl(i,1) / (Dx * dt);  
        cl11(ntpx,2) = cl11(ntpx,2) + fl(i,2) * fl(i,1) / (Dx * dt);
        cl11(ntpx,3) = cl11(ntpx,3) + fl(i,3) * fl(i,1) / (Dx * dt); 
    end
  elseif (npy < 1)
    ntpx = min(Nx,max(1,npx));
    if(fl(i, 3) < 0)
        cl21(ntpx,1) = cl21(ntpx,1) - fl(i,1) / (Dx * dt);   
        cl21(ntpx,2) = cl21(ntpx,2) + fl(i,2) * fl(i,1) / (Dx * dt);
        cl21(ntpx,3) = cl21(ntpx,3) + fl(i,3) * fl(i,1) / (Dx * dt); 
    end
  else %% Reconstruction densite et quantite de mouvement (pour affichage)
    d(npx,npy) = d(npx,npy)+ fl(i,1)/Dx/Dy;
    v(npx,npy,1) = v(npx,npy,1)+ fl(i,1)*fl(i,2)/Dx/Dy;
    v(npx,npy,2) = v(npx,npy,2)+ fl(i,1)*fl(i,3)/Dx/Dy;

    Ntmp = Ntmp+1;
    flt(Ntmp,:) = fl(i,:);
  end;
end;

%% Reconstruction de la vitesse a partir de la quantite de mouvement(pour affichage)
for i=1:Nx
  for j=1:Ny
	  if (d(i,j) > 0.0000000001)
            v(i,j,:) = v(i,j,:)/d(i,j);
          else 
	    v(i,j,:)=0.;
          end;
  end; 
end;

fl = flt;
