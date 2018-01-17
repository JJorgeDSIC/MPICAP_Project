function [X,T] = solver_AX_BX_C(A,B,C)

%Decomposite B = QTQ'
[Q,T]=schur(B);
D = C*Q;
Y = solver_AY_YT_D(A,T,D);

%XQ=Y
%XQQ'=YQ'
X = Y*Q';
end

