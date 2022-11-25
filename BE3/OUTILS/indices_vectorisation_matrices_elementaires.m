function [ii,jj,iiv,jjv,iir,jjr]=indices_vectorisation_matrices_elementaires(N)

%ii(1:N^2) = indices des lignes pour une matrice carrée d'une inconnue scalaire
%jj(1:N^2) = indices des colonnes pour une matrice carrée d'une inconnue scalaire
%iiv(1:4*N^2) = indices des lignes pour une matrice carrée d'une inconnue vectorielle
%jjv(1:4*N^2) = indices des colonnes pour une matrice carrée d'une inconnue vectorielle
%iir(1:2*N^2) = indices des lignes pour la matrice de rigidité
%jjr(1:2*N^2) = indices des colonnes pour la matrice de rigidité

ii = [];
jj = [];
iiv = [];
jjv = [];
iir = [];
jjr = [];

for i = 1:N
 iiv = [iiv [1:2*N] [1:2*N]];
 jjv = [jjv (2*(i-1)+1)*ones(1,2*N) (2*(i-1)+2)*ones(1,2*N)];
 ii  = [ii [1:N]];
 jj  = [jj i*ones(1,N)];
 iir = [iir [1:2*N]];
 jjr = [jjr i*ones(1,2*N)];
end


end