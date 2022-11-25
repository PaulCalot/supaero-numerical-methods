function fct_base = fct_base(i,A,x,y,ordre) %%valeur de la fonction de base i au point (x,y)

Vec = Vandermonde_mono(x,y,ordre);

fct_base = Vec*A(:,i);

end