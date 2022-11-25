function [Ep,Hp] = Eval_champs_ponctuels(E,H,r,s,ordre,A)

%E et H conntiennent les DDLS sur le triangle considér
N = (ordre +1)*(ordre+2)/2;
Ep =[0;0];
Hp = 0;
for k = 1:N
  fct_basek = fct_base(k,A,r,s,ordre);
  Ep = Ep +[E(2*(k-1)+1,1);E(2*(k-1)+2,1)]*fct_basek;
  Hp = Hp +H(k)*fct_basek; 
end 

end 