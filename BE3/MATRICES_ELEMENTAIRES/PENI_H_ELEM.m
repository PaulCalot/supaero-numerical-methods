function PENI_H = PENI_H_ELEM(A,ordre,Ng1D,xg1D,wg1D,xT,yT,voisT,xa,ya,la,nT)

% A = Coefficients monomiaux des fonctions de base
% (Ng1D,xg1D,wg1D) = formule d'intégrations sur [0,1]
% (xT,yT) = coordonées des sommets du triangle T
% INV_DFT = DF^*{-1}
% Liste des voisins de T. Si voisT(l) = 0 alors l'arete l de T est sur le bord!
% (xa,ya)(1:2,l) = coordonnées des extrémités de l'arete l
% la(l) = longueur de l'arete l 
% nT(l) = normale unitaire sortante associée à l'arete l de T


N = (ordre +1)*(ordre+2)/2; 
PENI_H = zeros(N,N);



for i=1:3
 if(voisT(i)>0) 
 for n = 1:Ng1D
  [x,y] = eval_Fa(xg1D(n),xa(:,i),ya(:,i)); 
  [r,s] = eval_inv_FT(x,y,xT,yT);
 
 for k =1:N
  fct_basek = fct_base(k,A,r,s,ordre);
  trace_basek = [nT(2)*fct_basek;-nT(1)*fct_basek];
  
  for l=1:N 
   fct_basel = fct_base(l,A,r,s,ordre);
   trace_basel = [nT(2)*fct_basel;-nT(1)*fct_basel];
   
   PENI_H (k,l) = PENI_H (k,l)+wg1D(n)*la(i)*trace_basel'*trace_basek;
   
  end 
 end 
 end 
end
end


end 