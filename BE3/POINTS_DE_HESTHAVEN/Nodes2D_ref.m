function [x,y] = Nodes2D_ref(ordre)

[x1,y1] = Nodes2D(ordre);

[r,s] = xytors(x1,y1);

x = (r+1)/2;
y = (s+1)/2;

end 