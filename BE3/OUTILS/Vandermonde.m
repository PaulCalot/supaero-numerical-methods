function [VD] = Vandermonde(x,y,ordre) %Matrice de Vendermond associée aux monomes x^iy^j, i+j<=ordre

N = (ordre+1)*(ordre+2)/2;
VD = zeros(N,N);
compt = 0;
for i =0:ordre
 for j=0:i
   compt = compt +1;
   VD(:,compt) = (x.^(i-j)).*(y.^j);
 end 
end  


end 