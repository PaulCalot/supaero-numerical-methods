
function [b] = PoissonRHS(nx,ny,hx,hy,xmin,ymin) 
% remplissage du second membre pour les conditions aux limites
for i=1:nx+1
    xi=(i-1)*hx+xmin;
   for j=1:ny+1
       
       yj=(j-1)*hy+ymin;
       k = i + (j-1)*(nx+1);  b(k)=0.;
       if(i>1&i<nx+1&j>1&j<ny+1)
            b(k)=sin(xi*pi)*sin(yj*pi);
               end
   end
end 


end