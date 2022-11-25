function frame = visu_solution_Film_E(E,comp,A,ordre,mesh,mesh_ref,Xmin,Xmax,Ymin,Ymax,Zmin,Zmax)

%comp = 1 pour Ex et comp = 2 pour Ey  
N = (ordre +1)*(ordre+2)/2;

mesh_visu.Nbtri = 0;
mesh_visu.Nbpts = 0;
if comp == 1 
 Exvisu =zeros(mesh.Nbtri*mesh_ref.Nbpts,1);
end 
if comp == 2 
 Eyvisu =zeros(mesh.Nbtri*mesh_ref.Nbpts,1);
end 
 
for T = 1:mesh.Nbtri
 xT = mesh.coor(1,mesh.Tri(1:3,T));
 yT = mesh.coor(2,mesh.Tri(1:3,T));
 corresp = zeros(mesh_ref.Nbpts,1);
  ET = E(2*N*(T-1)+1:2*N*(T-1)+2*N,1);
  HT = zeros(1:N);
 for l=1:mesh_ref.Nbpts
  [x,y] = eval_FT(mesh_ref.coor(1,l),mesh_ref.coor(2,l),xT,yT);
  mesh_visu.Nbpts = mesh_visu.Nbpts +1;
  corresp(l,1) = mesh_visu.Nbpts;
  mesh_visu.coor(:,mesh_visu.Nbpts) = [x;y];
 
  [Ep,Hp] = Eval_champs_ponctuels(ET,HT,mesh_ref.coor(1,l),mesh_ref.coor(2,l),ordre,A);
  
  if comp == 1
   Exvisu(mesh_visu.Nbpts,1) = Ep(1);
  end 
  if comp == 2
   Eyvisu(mesh_visu.Nbpts,1) = Ep(2);
  end  

 end 

 for l=1:mesh_ref.Nbtri
  mesh_visu.Nbtri = mesh_visu.Nbtri+1;
  mesh_visu.Tri(1,mesh_visu.Nbtri) = corresp(mesh_ref.Tri(1,l));
  mesh_visu.Tri(2,mesh_visu.Nbtri) = corresp(mesh_ref.Tri(2,l));
  mesh_visu.Tri(3,mesh_visu.Nbtri) = corresp(mesh_ref.Tri(3,l));
 end 
end 


if comp == 1  
 trimesh(mesh_visu.Tri',mesh_visu.coor(1,:)',mesh_visu.coor(2,:)',Exvisu,'FaceColor','interp');
end
if comp == 2  
 trimesh(mesh_visu.Tri',mesh_visu.coor(1,:)',mesh_visu.coor(2,:)',Eyvisu,'FaceColor','interp');
end

axis([Xmin Xmax Ymin Ymax Zmin Zmax]);
caxis([Zmin Zmax]);
view(0,90);

frame = getframe;

end 