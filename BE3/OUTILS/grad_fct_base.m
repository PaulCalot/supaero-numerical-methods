function grad_fct_base = grad_fct_base(i,A,x,y,ordre) 
%%valeur du gradient de la fonction de base i au point (x,y)

Vecx = Vandermonde_mono_x(x,y,ordre);
Vecy = Vandermonde_mono_y(x,y,ordre);

grad_fct_base = [Vecx*A(:,i);Vecy*A(:,i)];

end