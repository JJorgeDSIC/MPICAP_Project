function [Y] = solver_AY_YT_D(A,T,D)

N = size(A,1);
Y = zeros(N);
Daux = zeros(N,1);
I = eye(N);
i = 1;

while i<=N

    if i~=N & T(i+1,i) ~=0

        disp("No es triangular superior");
        
        for j=1:i
            D(:,i) = D(:,i) + Y(:,j)*T(j,i);
            D(:,i+1) = D(:,i+1) + Y(:,j)*T(j,i+1);
        end
        
        Dk = D(:,i);
        Ds = D(:,i+1);
        P = Dk/T(i+1,i);
        R = (A-T(i,i)*I)/T(i+1,i);
        Z = A*R-T(i,i+1)*I-R*T(i+1,i+1);
        W = Ds + A*P - P.*T(i+1,i+1);
        %Zx = W
        %x = Z\W
        Y(:,i) = Z\W;
        Y(:,i+1) = R*Y(:,i) - P;
        i=i+1;
        
    else 
        
        Z = A-(T(i,i)*I);
        for j=1:i
            D(:,i) = D(:,i) + Y(:,j)*T(j,i);
        end
        b = D(:,i);
        %x = Y(:i);
        %Zx = b
        %x = Z\b
        Y(:,i) = Z\b;
        
    end
    i=i+1;
end

end

