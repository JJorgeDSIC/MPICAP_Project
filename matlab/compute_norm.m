function [value] = compute_norm(A,B,C,X)
tot=A*X-X*B-C;
value=norm(tot,'fro');
end

