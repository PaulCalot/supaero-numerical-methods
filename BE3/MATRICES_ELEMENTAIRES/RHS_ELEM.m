function JELEM = RHS_ELEM(xc,yc,r0,A,ordre,Nq,xq,yq,wq,xT,yT)


N= (ordre +1)*(ordre+2)/2; 
JELEM = zeros(2*N,1);

for n = 1:Nq
 [x,y] = eval_FT(xq(n),yq(n),xT,yT);
 curl_Js = eval_curl_point_source_space(x,y,xc,yc,r0);
 
 for k = 1:N
     
   fct_basek = fct_base(k,A,xq(n),yq(n),ordre);
  
   JELEM(2*(k-1)+1,1) = JELEM(2*(k-1)+1,1) + wq(n)*curl_Js(1,1)*fct_basek;
   JELEM(2*(k-1)+2,1) = JELEM(2*(k-1)+2,1) + wq(n)*curl_Js(2,1)*fct_basek;
   

 end 
end 




end 