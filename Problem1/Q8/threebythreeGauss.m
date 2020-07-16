function k = threebythreeGauss(B, D, J,w)

k = zeros(16,16);
s = [-sqrt(3/5) 0 sqrt(3/5)];
t = [-sqrt(3/5) 0 sqrt(3/5)];
W = [5/9 8/9 5/9];

for i=1:3
    for j = 1:3
        B_ij = double(B(s(i),t(j)));
        J_ij = double(J(s(i),t(j)));
        k = k + W(i)*W(j)*B_ij'*D*B_ij*det(J_ij);
    end
end        

k = w*k;
end