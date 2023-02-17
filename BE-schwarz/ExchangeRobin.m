
function [b] = ExchangeRobin(b,u,nx,ny,ip,jp,px,py,idelta,alpha,hx,hy) 
% remplissage du second membre pour les conditions aux limites
nxp=nx(ip,jp);nyp=ny(ip,jp);
if(jp<py) 
for i=1:nxp+1
    % envoie et reception du domaine voisin Nord: jp+1
   k = i + (nyp)*(nxp+1)   ;       
   b(k)= u {ip,jp+1}( i + (2*idelta)*(nx(ip,jp+1)+1)  )*(1+alpha*hy)-u {ip,jp+1}( i + (2*idelta-1)*(nx(ip,jp+1)+1)  );
end 
end

if(jp>1) 
for i=1:nxp+1
    % envoie et reception du domaine voisin Sud: jp-1
   k = i ;      

  b(k)= u {ip,jp-1}( i + (ny(ip,jp-1)-2*idelta)*(nx(ip,jp-1)+1)  )*(1+alpha*hy) -u{ip,jp-1}( i + (ny(ip,jp-1)-2*idelta+1)*(nx(ip,jp-1)+1)  ); 
end 
end

if(ip<px) 
for j=1:nyp+1
    % envoie et reception du domaine voisin Est: ip+1
 
     k = nxp+1 + (j-1)*(nxp+1); 
     b(k)=u{ip+1,jp}( 2*idelta+1 + (j-1)*(nx(ip+1,jp)+1))*(1+alpha*hx)-u{ip+1,jp}( 2*idelta + (j-1)*(nx(ip+1,jp)+1));
end
end


   if(ip>1)
     % envoie et reception du domaine voisinouest: ip-1
     for j=1:nyp+1

   k = 1 + (j-1)*(nxp+1);    
   b(k)=u{ip-1,jp}(nx(ip-1,jp)-2*idelta+1 + (j-1)*(nx(ip-1,jp)+1))*(1+alpha*hx)-u{ip-1,jp}(nx(ip-1,jp)-2*idelta+2 + (j-1)*(nx(ip-1,jp)+1));
     
     end
    

 end 
 
 end