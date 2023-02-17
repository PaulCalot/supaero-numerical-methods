clear all;close all
dx=1.;dy=1;
% nbre depoints par sous domaine et nommbre de sousdomaines (dans les 
% 2 directions
nxb=32;nyb=32;px=2;py=2;idelta=4;
ntx=nxb*px;nty=nyb*py;
hx=dx/ntx;hy=dy/nty;
dpx=dx/px;dpy=dy/py;
% calcul de la solution exacte
[A,N]=PoissonInit(ntx,nty,hx,hy); 
b=zeros(N,1);u=zeros(N,1);
[b]=PoissonRHS(ntx,nty,hx,hy,0.,0.);
u=A\b';
ugrido=reshape(u,ntx+1,nty+1);
contourf(ugrido)
% resolution sur deux domaines par la methode de Schwarz
% idelta est le nombre de mailles  de recouvrement
%
    Ap=cell(px,py); bp=cell(px,py);up=cell(px,py);
for i=1:px
    for j=1:py
        xmin(i,j)?
        ymin(i,j)?;
        nx(i,j)=?;ny(i,j)=?;
       
        [AA,Np(i,j)]=PoissonInit(nx(i,j),ny(i,j),hx,hy);
        Ap{i,j}=AA;

        % calcul des second membres locaux 
         
         bp{i,j}=?;

        % initialisation de la solution dans chaque sous domaine
        up{i,j}=?;
    end
end

;epsilon=1d-4;it=0;upn=up;residu=1.;itmax=500;
% % Algorithme de Schwarz
    while(residu>epsilon&it<itmax)
        %echange avec le ou les sous-domaines voisins
         
        % calcul du residu
      uu= Recons(up,nx,ny,nxb,nyb, px,py,idelta) ;
      res= norm(uu-ugrido)/norm(ugrido)
         
        up=upn;
        it=it+1, residu
        residu(it)=res;
       % pause
    end 
       
%     reconstruction de la solution global
for i=1:px
       for j=1:py
           ugrid = reshape(up{i,j},nx(i,j)+1,ny(i,j)+1);  
%contourf(ugrid'),colorbar
%pause
       end
end
contourf(uu'),colorbar
   semilogy(residu);residu(it)/residu(it-1)
