function Vec = Vandermonde_mono(x,y,ordre)

N = (ordre+1)*(ordre+2)/2;
Vec= zeros(1,N);
compt = 0;
for i =0:ordre
 for j=0:i
   compt = compt+1;
   Vec(1,compt) = x^(i-j)*y^j;
 end 
end  


end 