function [N,x,y,A] = ddl_ET_fct_base_SUR_triangle_reference(ordre)

% Nombre de points d'interpolation sur le triangle de référence
%(x,y) sont les coordonnées des points d'interpolation
%A(:,i) contient les coefficients monomiaux décrivant la i ème fonctions 

N = (ordre+1)*(ordre+2)/2; 
[x,y] = Nodes2D_ref(ordre);  %Points d'interpolations
[VD] = Vandermonde(x,y,ordre);
A = VD\diag(ones(N,1)); %Coefficients monomiaux des fonctions de base

end 