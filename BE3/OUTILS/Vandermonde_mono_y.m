function Vec = Vandermonde_mono_y(x,y,ordre) %derivée en y

N = (ordre+1)*(ordre+2)/2;
Vec= zeros(1,N);
compt = 0;
for i =0:ordre
 for j=0:i
   compt = compt +1;
   if (j>0)
    Vec(1,compt) = j*x^(i-j)*y^(j-1);
   else
    Vec(1,compt) = 0;
   end
 end 
end  


end 