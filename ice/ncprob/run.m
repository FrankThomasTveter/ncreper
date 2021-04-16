# mathics
apt -get  install  python-dev libsqlite3-dev python-setuptools
sudo apt install python-pip
pip install mathics
#

mathics

n[4] := 3*JJ;
n[n_] := If[Mod[n,2]==1,0, 2^(-n/2)*(n!)/((n/2)!)];

rl0 = { E[a_+b_]:>E[a]+E[b],
   E[a_Constant]:>a,
   E[a_Constant^2]:>a^2,
   E[a_Constant^3]:>a^3,
   E[a_Constant^4]:>a^4,
   E[a_Constant*b_]:>a*E[b],
   E[a_Constant^2*b_]:>a^2*E[b],
   E[a_Constant^3*b_]:>a^3*E[b],
   E[a_Constant^4*b_]:>a^4*E[b],
   E[a_Integer*b_]:>a*E[b]};

rl1 = {y -> Constant[a]+ Constant[b]*u + Constant[c]*u^2 + Constant[d]*u^3};

rl2 = {E[u^nn_] :> n[nn], E[u] :> n[1],Constant[a_]:>a, E[y^nn_]:>e[nn], E[y]:>e[1], E[1]:>1, E[-1]:>-1};

s1 = ((E[y//.rl1] - E[y])//.rl0)//.rl2

s2 = ((E[ExpandAll[(y^2)//.rl1]] - E[y^2])//.rl0)//.Join[rl2]

s3 = Simplify[((E[ExpandAll[(y^3)//.rl1]] - E[y^3])//.rl0)//.Join[rl2]];

s4 = ((E[ExpandAll[(y^4)//.rl1]] - E[y^4])//.rl0)//.Join[rl2];


s1da = D[s1,a];
s1db = D[s1,b];
s1dc = D[s1,c];
s1dd = D[s1,d];

s2da = D[s2,a];
s2db = D[s2,b];
s2dc = D[s2,c];
s2dd = D[s2,d];

s3da = D[s3,a];
s3db = D[s3,b];
s3dc = D[s3,c];
s3dd = D[s3,d];

s4da = D[s4,a];
s4db = D[s4,b];
s4dc = D[s4,c];
s4dd = D[s4,d];






eqs4b = Simplify[eqs4]//.{d^n_:>dd[n]}

cost[i_]:=(((E[ExpandAll[(y^i)//.rl1]]//.rl0) - E[i])//.rl2);


			 
f[i_] := (E[ExpandAll[(((y-Constant[mu]))^i)]]//.rl0)//.rl2;

g[i_] := (E[ExpandAll[(((y-Constant[aa]-Constant[bb]*x))^i)]]//.rl0)//.{Constant[a_]:>a,E[-1]:>-1, E[1] :> 1};


f[i_] := (E[ExpandAll[(((y-Constant[mu])*Constant[gamma])^i)]]//.rl0)//.rl2;







#
# maxima
#


normalExpect(nn) :=
block( [],
       if mod(nn,2) = 1
       then return(0)
       else return(2^(-nn/2) * (nn!)/((nn/2)!))
     );
u: a+b*x+c*x^3;



