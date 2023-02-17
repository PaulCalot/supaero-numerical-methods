function [uu] = Recons(up,nx,ny,nxb,nyb,px,py,idelta) 
% reconstruction de la solution

  kk=0; uu=0;
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

 end