N=5;
A = zeros(5);
B = zeros(5);
X = zeros(5);
C = zeros(5);

for i=1:N
    for j=1:N
        A(i,j)=i-1;
        B(i,j)=i-1;
        X(i,j)=i-1;
        C(i,j)=i-1;
    end
end

C = A * B