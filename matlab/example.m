N = 1000;
%Basic problem
A = rand(N);
B = rand(N);
X = rand(N);
C=A*X-X*B;

[Xc, Tc]=solver_AX_XB_C(A,B,C);

compute_norm(A,B,C,X)
compute_norm(A,B,C,Xc)

norm(X-Xc,'fro')

%Second problem
% A = rand(N)
% T = triu(rand(N))
% T(3,2) = 0.5;
% Y = rand(N)
% D=A*Y-Y*T;
% 
% [Yc]=solver_AY_YT_D(A,T,D)
% 
% compute_norm(A,T,D,Y)
% compute_norm(A,T,D,Yc)