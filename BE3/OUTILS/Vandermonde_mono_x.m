function Vec = Vandermonde_mono_x(x,y,ordre) %derivée en x

N = (ordre+1)*(ordre+2)/2;
Vec= zeros(1,N);
compt = 0;
for i =0:ordre
 for j=0:i
   compt = compt +1;
   if(i-j>0) 
    Vec(1,compt) = (i-j)*x^(i-j-1)*y^j;
   else
    Vec(1,compt) = 0;
   end  
 end 
end  


end 