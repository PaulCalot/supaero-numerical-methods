function [A,N] = PoissonInit(nx,ny,hx,hy)
%
%  Cette routine initialise la matrice de Poisson 
% dans un rectangle avec des conditions aux limites de 
% type Dirichlet n



N = (nx+1)*(ny+1);


A = sparse(zeros(N,N));  



for j=2:ny
  for i=2:nx
    k = i + (j-1)*(nx+1);    
    A(k,k) = 2/hx^2+2/hy^2;
    A(k,k-1) =- 1/hx^2;      
    A(k,k+1) = -1/hx^2; 
    A(k,k-(nx+1)) =- 1/hy^2;
    A(k,k+(nx+1)) = -1/hy^2; 
  end;
end;
% Prise en compte des conditions de Dirichlet (0)

    for i=1: nx+1
        % condition Nord c est a dire j=ny+1

        k = i + (ny)*(nx+1);   
        
        
        A(k,:)=0;A(k,k)=1.;
       
        
         % condition Sud c est a dire j=1

        k = i;   
        
       
        A(k,:)=0;A(k,k)=1.;
       
        
      
    end     
    
   for j=1: ny+1
        % condition West c est a dire i=nx+1

        k = nx+1 + (j-1)*(nx+1);   

        A(k,:)=0;A(k,k)=1.;
       
          
         % condition Est c est a dire i=1

        k = 1+(j-1)*(nx+1);   
       
        A(k,:)=0;A(k,k)=1.;
         
   end 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
