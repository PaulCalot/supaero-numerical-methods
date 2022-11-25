function SELEMI = SAUTI_ELEM(A,ordre,Ng1D,xg1D,wg1D,xT,yT,voisT,xa,ya,la,nT)

% A = Coefficients monomiaux des fonctions de base
% (Ng1D,xg1D,wg1D) = formule d'intégrations sur [0,1]
% (xT,yT) = coordonées des sommets du triangle T
% INV_DFT = DF^*{-1}
% Liste des voisins de T. Si voisT(l) = 0 alors l'arete l de T est sur le bord!
% (xa,ya)(1:2,l) = coordonnées des extrémités de l'arete l
% la(l) = longueur de l'arete l 
% nT(l) = normale unitaire sortante associée à l'arete l de T


N = (ordre +1)*(ordre+2)/2; 
SELEMI = zeros(2*N,N);
beta = 1/2;


for i=1:3
 if(voisT(i)>0) 
 for n = 1:Ng1D
  [x,y] = eval_Fa(xg1D(n),xa(:,i),ya(:,i)); 
  [r,s] = eval_inv_FT(x,y,xT,yT);
 
 for k =1:N
  fct_basek = fct_base(k,A,r,s,ordre);
  basek_x = [fct_basek;0];
  basek_y = [0;fct_basek];
  
  
  for l=1:N 
   fct_basel = fct_base(l,A,r,s,ordre);
   trace_basel = [nT(2,i)*fct_basel;-nT(1,i)*fct_basel];
 
   SELEMI(2*(k-1)+1,l) = SELEMI(2*(k-1)+1,l)+wg1D(n)*la(i)*trace_basel'*basek_x*beta;
   SELEMI(2*(k-1)+2,l) = SELEMI(2*(k-1)+2,l)+wg1D(n)*la(i)*trace_basel'*basek_y*beta;
 
  end 
 end 
 end 
end
end


end 