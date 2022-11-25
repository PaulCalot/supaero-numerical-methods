function [NNZ_IME,NNZ_IMH,NNZ_R,NNZ_P_E,NNZ_P_H,NNZ_SM] =NB_NON_ZEROS_MATRICES(N,mesh,Type_arete,CL_SM)

%%Calcul du nombre  de non-zeros de chaque matrice
NNZ_IME = (2*N)^2*mesh.Nbtri;
NNZ_IMH = N^2*mesh.Nbtri;
NNZ_R =  2*N^2*(mesh.Nbtri+nnz(mesh.voisin(:,:)));
NNZ_P_E =4*N^2*(mesh.Nbtri+nnz(mesh.voisin(:,:)));
NNZ_P_H =N^2*(mesh.Nbtri+nnz(mesh.voisin(:,:)));

NNZ_SM =0;
if CL_SM == 1
for T = 1:mesh.Nbtri
 n1 =mesh.Tri(1,T);
 n2 =mesh.Tri(2,T);
 n3 =mesh.Tri(3,T);
 %%Récupération des informations sur les voisins de T
 voisT = mesh.voisin(:,T);
 na = [n2 n3 n1; n3 n1 n2];
 trouve =0;
 for i=1:3
  if voisT(i)==0 && CL_SM ==1 && Type_arete(na(1,i),na(2,i))==1
   trouve =1;
  end
 end
 if trouve == 1
  NNZ_SM = NNZ_SM + N^2;
 end
 end
end
