function SELEME = SAUTE_ELEM(A,ordre,Ng1D,xg1D,wg1D,xT,yT,xvT,yvT,nT,xa,ya,la)

% A = Coefficients monomiaux des fonctions de base
% (Ng1D,xg1D,wg1D) = formule d'intégrations sur [0,1]
% (xT,yT) = coordonées des sommets du triangle T
% INV_DFT = DF_T^*{-1}
% (xvT,yvT) = coordonées des sommets d'un voisin de triangle T
% nT = normale unitaire sortante à l'arete considéré
% (xa,ya) = coordonnées des 2 extrémités de l'arete considéré

N = (ordre +1)*(ordre+2)/2; 
SELEME = zeros(2*N,N);
beta = 1/2;


for n = 1:Ng1D
  [x,y] = eval_Fa(xg1D(n),xa(:),ya(:));
  [r,s] = eval_inv_FT(x,y,xT,yT); % F_T^{-1}(x,y)
  [rv,sv] = eval_inv_FT(x,y,xvT,yvT); %F_vT^{-1}(x,y)
 
 
 for k =1:N
  fct_basek = fct_base(k,A,r,s,ordre);
  basek_x = [fct_basek;0];
  basek_y = [0;fct_basek];
  
  
  for l=1:N 
  fct_basel = fct_base(l,A,rv,sv,ordre);
  trace_basel = [nT(2)*fct_basel;-nT(1)*fct_basel];
 
   SELEME(2*(k-1)+1,l) = SELEME(2*(k-1)+1,l)+beta*wg1D(n)*la*trace_basel'*basek_x;
   SELEME(2*(k-1)+2,l) = SELEME(2*(k-1)+2,l)+beta*wg1D(n)*la*trace_basel'*basek_y;
 
  end 
 end 
end 
end 