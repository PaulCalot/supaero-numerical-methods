function [mesh_ref] = creation_maillage_raffine_reference(p)
%p = le taux de raffinement

x =[0:1/2^p:1];
y =[0:1/2^p:1];
%%Creation des points
n = 0;
for j = 1:2^p+2
 for i = 1:2^p+2-j
  n = n+1;
  mesh_ref.coor(1,n) = x(i);
  mesh_ref.coor(2,n) = y(j);
  corresp(i,j) = n;
 end 
end 
mesh_ref.Nbpts = n;
%%Creation des triangles
n = 0;
for j = 1:2^p+1
 for i = 1:2^p+1-j
  if((j<2^p+1) && (i<2^p+1-j)) %on construit 2 triangles
  n = n+1;
  mesh_ref.Tri(1,n) = corresp(i,j);
  mesh_ref.Tri(2,n) = corresp(i+1,j);
  mesh_ref.Tri(3,n) = corresp(i,j+1);
  n = n+1;
  mesh_ref.Tri(1,n) = corresp(i+1,j);
  mesh_ref.Tri(2,n) = corresp(i+1,j+1);
  mesh_ref.Tri(3,n) = corresp(i,j+1);
  
  else % on construit un seul triangle
  n = n+1;
  mesh_ref.Tri(1,n) = corresp(i,j);
  mesh_ref.Tri(2,n) = corresp(i+1,j);
  mesh_ref.Tri(3,n) = corresp(i,j+1);
  end
 
 end 
end 
mesh_ref.Nbtri = n;

end 