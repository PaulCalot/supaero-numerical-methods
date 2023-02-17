clear all;close all
dx=1.;dy=1;
% nvbre depoints par sous domaine et nommbre de sousdomaines (dans les 
% 2 directions
npx=2;npy=2;
nxb=40;nyb=80;px=2;py=1;idelta=2;
iii=0;
ncas=20;
for iii=1:ncas
idelta=iii;
ntx=nxb*px;nty=nyb*py;

hx=dx/ntx;hy=dy/nty;
dpx=dx/px;dpy=dy/py;
% calcul de la solution exacte
[A,N]=PoissonInit(ntx,nty,hx,hy); 
b=zeros(N,1);u=zeros(N,1);
[b]=PoissonRHS(ntx,nty,hx,hy,0.,0.);

u=A\b';

ugrido = reshape(u,ntx+1,nty+1);  
%contourf(ugrido'),colorbar
% resolution sur deux domaines par la methode de Schwarz
% idelta est le nombre de mailles  de recouvrement
%
    Ap=cell(px,py); bp=cell(px,py);up=cell(px,py); ApL=cell(px,py); ApU=cell(px,py); 
for i=1:px
    for j=1:py
        xmin(i,j)=max((i-1)*dpx-idelta*hx,0.);
        ymin(i,j)=max((j-1)*dpy-idelta*hy,0.);
        nx(i,j)=nxb+2*idelta;ny(i,j)=nyb+2*idelta;
        if(i==1)  nx(i,j)=nx(i,j)-idelta;end
        if(i==px) nx(i,j)=nx(i,j)-idelta;end  
        if(j==1)  ny(i,j)=ny(i,j)-idelta;end
        if(j==py) ny(i,j)=ny(i,j)-idelta; end
        ii=i+(j-1)*px;
        [AA,Np(i,j)]=PoissonInit(nx(i,j),ny(i,j),hx,hy); 
        [ApL{i,j},ApU{i,j}]=lu(AA);
        Ap{i,j}=AA;
        [bb]=PoissonRHS(nx(i,j),ny(i,j),hx,hy,xmin(i,j),ymin(i,j));
         bp{i,j}=bb;
        up{i,j}=ApL{i,j}\bb';up{i,j}=ApU{i,j}\up{i,j};
    end
end

;epsilon=1d-6;it=0;upn=up;itmax=500;t=cputime;residu0=norm(u)
residu=residu0;it=1;residu(1:itmax)=0;residu(1)=1;

% % Algorithme de Schwarz
    while(residu(it)>epsilon*residu0&it<itmax)
        %echange avec le ou les sous-domaines voisins
        for ip=1:px
            for jp=1:py
               [bp{ip,jp}]=Exchange(bp{ip,jp},up,nx,ny,ip,jp,px,py,idelta) ;
             
            end
        end
 
       for ip=1:px
            for jp=1:py
                  
                bb=bp{ip,jp};
               
                  upn{ip,jp}=Ap{ip,jp}\bb';   
                 
            end
       end
       
        up=upn;
        it=it+1;
    
          
%     reconstruction de la solution global

    kk=0; clear uu;
    for i=1:px
           for j=1:py
               
           lmin=1;
           if(j>1)lmin=1+idelta;end
           lmax=ny(i,j)+1;
           if(j<py)lmax=ny(i,j)+1-idelta;end
            kmin=1;
           if(i>1)kmin=1+idelta;end
           kmax=nx(i,j)+1;
           if(i<px)kmax=nx(i,j)+1-idelta;end
          
    
           for k=kmin:kmax
             for l=lmin:lmax  
              kk=(i-1)*nxb+k-kmin+1;
              ll=(j-1)*nyb+l-lmin+1;
              uu(kk,ll)=up{i,j}(k+(l-1)*(nx(i,j)+1));
             end  
           end
           
       end
    end 
     residu(it)= norm(uu-ugrido)/norm(ugrido);
        
    end
    it
 
    npt(iii)=idelta;nit(iii)=it;


    
end
 semilogy(npt,nit);hold on