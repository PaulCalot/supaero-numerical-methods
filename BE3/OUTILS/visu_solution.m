function visu_solution(E,H,A,ordre,mesh,mesh_ref)


N = (ordre +1)*(ordre+2)/2;

mesh_visu.Nbtri = 0;
mesh_visu.Nbpts = 0;
Exvisu = zeros(mesh.Nbtri*mesh_ref.Nbpts,1);
Eyvisu = zeros(mesh.Nbtri*mesh_ref.Nbpts,1);
Hvisu  = zeros(mesh.Nbtri*mesh_ref.Nbpts,1);

for T = 1:mesh.Nbtri
 xT = mesh.coor(1,mesh.Tri(1:3,T));
 yT = mesh.coor(2,mesh.Tri(1:3,T));
 corresp = zeros(mesh_ref.Nbpts,1);
  ET = E(2*N*(T-1)+1:2*N*(T-1)+2*N,1);
  HT = H(N*(T-1)+1:N*(T-1)+N,1);
 for l=1:mesh_ref.Nbpts
  [x,y] = eval_FT(mesh_ref.coor(1,l),mesh_ref.coor(2,l),xT,yT);
  mesh_visu.Nbpts = mesh_visu.Nbpts +1;
  corresp(l,1) = mesh_visu.Nbpts;
  mesh_visu.coor(:,mesh_visu.Nbpts) = [x;y];
 
  [Ep,Hp] = Eval_champs_ponctuels(ET,HT,mesh_ref.coor(1,l),mesh_ref.coor(2,l),ordre,A);
  Exvisu(mesh_visu.Nbpts,1) = Ep(1);
  Eyvisu(mesh_visu.Nbpts,1) = Ep(2);
  Hvisu(mesh_visu.Nbpts,1)  = Hp;
 end 

 for l=1:mesh_ref.Nbtri
  mesh_visu.Nbtri = mesh_visu.Nbtri+1;
  mesh_visu.Tri(1,mesh_visu.Nbtri) = corresp(mesh_ref.Tri(1,l));
  mesh_visu.Tri(2,mesh_visu.Nbtri) = corresp(mesh_ref.Tri(2,l));
  mesh_visu.Tri(3,mesh_visu.Nbtri) = corresp(mesh_ref.Tri(3,l));
 end 
end 


figure(1);
%legend("Composante Ex")
trimesh(mesh_visu.Tri',mesh_visu.coor(1,:)',mesh_visu.coor(2,:)',Exvisu);
figure(2);
%legend("Composante Ey")
trimesh(mesh_visu.Tri',mesh_visu.coor(1,:)',mesh_visu.coor(2,:)',Eyvisu);
figure(3);
%legend("H")
trimesh(mesh_visu.Tri',mesh_visu.coor(1,:)',mesh_visu.coor(2,:)',Hvisu);




end 