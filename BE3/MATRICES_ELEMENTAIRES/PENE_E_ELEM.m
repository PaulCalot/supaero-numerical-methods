function PENE_E = PENE_E_ELEM(A,ordre,Ng1D,xg1D,wg1D,xT,yT,xvT,yvT,nT,xa,ya,la)

% A = Coefficients monomiaux des fonctions de base
% (Ng1D,xg1D,wg1D) = formule d'intégrations sur [0,1]
% (xT,yT) = coordonées des sommets du triangle T
% INV_DFT = DF_T^*{-1}
% (xvT,yvT) = coordonées des sommets d'un voisin de triangle T
% nT = normale unitaire sortante à l'arete considéré
% (xa,ya) = coordonnées des 2 extrémités de l'arete considéré

N = (ordre +1)*(ordre+2)/2; 
PENE_E = zeros(2*N,2*N);

for n = 1:Ng1D
  [x,y] = eval_Fa(xg1D(n),xa(:),ya(:));
  [r,s] = eval_inv_FT(x,y,xT,yT); % F_T^{-1}(x,y)
  [rv,sv] = eval_inv_FT(x,y,xvT,yvT); %F_vT^{-1}(x,y)
 
 
 for k =1:N
  fct_basek = fct_base(k,A,r,s,ordre);
  tracek_x = -nT(2)*fct_basek;
  tracek_y =  nT(1)*fct_basek;
  
  for l=1:N 
  fct_basel = fct_base(l,A,rv,sv,ordre);
  tracel_x = -nT(2)*fct_basel;
  tracel_y =  nT(1)*fct_basel;
 
  PENE_E(2*(k-1)+1,2*(l-1)+1) = PENE_E(2*(k-1)+1,2*(l-1)+1)+wg1D(n)*la*tracel_x*tracek_x;
  PENE_E(2*(k-1)+1,2*(l-1)+2) = PENE_E(2*(k-1)+1,2*(l-1)+2)+wg1D(n)*la*tracel_y*tracek_x;
  PENE_E(2*(k-1)+2,2*(l-1)+1) = PENE_E(2*(k-1)+2,2*(l-1)+1)+wg1D(n)*la*tracel_x*tracek_y;
  PENE_E(2*(k-1)+2,2*(l-1)+2) = PENE_E(2*(k-1)+2,2*(l-1)+2)+wg1D(n)*la*tracel_y*tracek_y;
 
 
  end 
 end 
end 
end 