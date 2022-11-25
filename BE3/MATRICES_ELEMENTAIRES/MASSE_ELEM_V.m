function MELEM = MASSE_ELEM_V(A,ordre,Nq,xq,yq,wq)

N = (ordre +1)*(ordre+2)/2; 
MELEM = zeros(2*N,2*N);
for n = 1:Nq
 for k = 1:N
   fct_basek = fct_base(k,A,xq(n),yq(n),ordre);
  for l = 1:N
   fct_basel = fct_base(l,A,xq(n),yq(n),ordre);
   MELEM(2*(k-1)+1,2*(l-1)+1) = MELEM(2*(k-1)+1,2*(l-1)+1) + wq(n)*fct_basel*fct_basek;
   MELEM(2*(k-1)+2,2*(l-1)+2) = MELEM(2*(k-1)+2,2*(l-1)+2) + wq(n)*fct_basel*fct_basek;
  end 
  end 
end 

end 