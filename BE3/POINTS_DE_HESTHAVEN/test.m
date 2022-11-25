close all
format long
fd = inline (’sqrt ( sum ( p .^2 ,2) ) -1 ’ , ’p ’) ;
[p, t ]= distmesh2d (fd , @huniform ,0.2 ,[ -1 , -1;1 ,1] ,[]) ;