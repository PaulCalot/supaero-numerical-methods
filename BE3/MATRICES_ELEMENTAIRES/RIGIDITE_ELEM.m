function RELEM = RIGIDITE_ELEM(A,ordre,Nq,xq,yq,wq,INV_DFT)

N = (ordre +1)*(ordre+2)/2; 
RELEM = zeros(2*N,N);

for n = 1:Nq
  for k = 1:N
     fct_basek = fct_base(k,A,xq(n),yq(n),ordre);
      
     base_x = [fct_basek;0];
     base_y = [0;fct_basek];
     
     for l = 1:N
      grad_l = grad_fct_base(l,A,xq(n),yq(n),ordre);
      gradT_l = INV_DFT*grad_l;
      curl_l =  [gradT_l(2);-gradT_l(1)];
         
      RELEM(2*(k-1)+1,l) = RELEM(2*(k-1)+1,l) + wq(n)*curl_l'*base_x;
      RELEM(2*(k-1)+2,l) = RELEM(2*(k-1)+2,l) + wq(n)*curl_l'*base_y;
      
  end 
  end 
end 
 end