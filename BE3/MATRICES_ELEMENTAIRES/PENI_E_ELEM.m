function PENI_E = PENI_E_ELEM(A,ordre,Ng1D,xg1D,wg1D,xT,yT,voisT,xa,ya,la,nT)

% A = Coefficients monomiaux des fonctions de base
% (Ng1D,xg1D,wg1D) = formule d'intégrations sur [0,1]
% (xT,yT) = coordonées des sommets du triangle T
% INV_DFT = DF^*{-1}
% Liste des voisins de T. Si voisT(l) = 0 alors l'arete l de T est sur le bord!
% (xa,ya)(1:2,l) = coordonnées des extrémités de l'arete l
% la(l) = longueur de l'arete l 
% nT(l) = normale unitaire sortante associée à l'arete l de T


N = (ordre +1)*(ordre+2)/2; 
PENI_E = zeros(2*N,2*N);



for i=1:3
 if(voisT(i)>0) 
 for n = 1:Ng1D
  [x,y] = eval_Fa(xg1D(n),xa(:,i),ya(:,i)); 
  [r,s] = eval_inv_FT(x,y,xT,yT);
 
 for k =1:N
  fct_basek = fct_base(k,A,r,s,ordre);
  tracek_x = -nT(2,i)*fct_basek;
  tracek_y =  nT(1,i)*fct_basek;
  
  for l=1:N 
   fct_basel = fct_base(l,A,r,s,ordre);
   tracel_x = -nT(2,i)*fct_basel;
   tracel_y =  nT(1,i)*fct_basel;
 
   PENI_E (2*(k-1)+1,2*(l-1)+1) = PENI_E (2*(k-1)+1,2*(l-1)+1)+wg1D(n)*la(i)*tracel_x*tracek_x;
   PENI_E (2*(k-1)+1,2*(l-1)+2) = PENI_E (2*(k-1)+1,2*(l-1)+2)+wg1D(n)*la(i)*tracel_y*tracek_x;
   PENI_E (2*(k-1)+2,2*(l-1)+1) = PENI_E (2*(k-1)+2,2*(l-1)+1)+wg1D(n)*la(i)*tracel_x*tracek_y;
   PENI_E (2*(k-1)+2,2*(l-1)+2) = PENI_E (2*(k-1)+2,2*(l-1)+2)+wg1D(n)*la(i)*tracel_y*tracek_y;
 
  end 
 end 
 end 
end
end


end 