function [mesh,Type_arete] = lecture_mesh(fich_mesh)

%mesh.Nbtri = nombre de triangles
%mesh.Nbpts = nombres de noeuds
%mesh.coor(1:2,1:Nbpts) = coordonnées des noeuds
%mesh.Tri(1:3,1:Nbtri) = connectique des triangles (par rapport aux noeuds)
%mesh.voisin(l=1:3,i) = numéro du triangle voisin de i opposé au noeud l
%mesh.boundary = connectique des arêtes où est imposée la condition de
%Silver-Müller
%Type_arete(n1,n2) = le type de l'arête [n1,n2]
%Si Type_arete(n1,n2) =1 alors [n1,n2] est une arête sujjete à une
%condition de Sliver-Müller. Dans le cas contraire, Type_arete(n1,n2) =0

load(fich_mesh); 
mesh.Nbtri = size(Tri,1);
mesh.Nbpts = size(coor,1);
mesh.coor  = coor';
mesh.Tri   = Tri';
mesh.voisin = voisin';
mesh.boundary = frontiere;

Type_arete = zeros(mesh.Nbpts,mesh.Nbpts);
for i=1:size(mesh.boundary,1)
 Type_arete(mesh.boundary(i,1),mesh.boundary(i,2)) =1;
 Type_arete(mesh.boundary(i,2),mesh.boundary(i,1)) =1;
end

end 