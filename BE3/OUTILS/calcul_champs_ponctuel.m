function [Es, Hs ] = calcul_champs_ponctuel(ET,HT,xsref,ysref,A,ordre)

N= (ordre +1)*(ordre+2)/2; 
Es = [0; 0];
Hs = 0;
for k = 1:N
 fct_basek = fct_base(k,A,xsref,ysref,ordre);
 Es(1,1) = Es(1,1) + ET(2*(k-1)+1)*fct_basek;   
 Es(2,1) = Es(2,1) + ET(2*(k-1)+2)*fct_basek;   
 Hs      = Hs + HT(k)*fct_basek;   
end  
   
end 