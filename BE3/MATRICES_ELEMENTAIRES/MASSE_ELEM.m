function MELEM = MASSE_ELEM(A,ordre,Nq,xq,yq,wq)

N = (ordre +1)*(ordre+2)/2; 
MELEM = zeros(N,N);
for n = 1:Nq
 for k = 1:N
   fct_basek = fct_base(k,A,xq(n),yq(n),ordre);
  for l = 1:N
   fct_basel = fct_base(l,A,xq(n),yq(n),ordre);
   MELEM(k,l) = MELEM(k,l) + wq(n)*fct_basel*fct_basek;
  end 
  end 
end 

end 