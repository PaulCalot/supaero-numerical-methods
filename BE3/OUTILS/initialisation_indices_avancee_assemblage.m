function [kk,kkv,kkr,kkpe,kkph,kk_SM] = initialisation_indices_avancee_assemblage(N)


kk  = 1:N^2;
kkv = 1:4*N^2;
kkr = 1:2*N^2;
kkpe = 1:4*N^2;
kkph = 1:N^2;
kk_SM = 1:N^2;

end