function k = twobytwoGauss(B, D, J,w)

k = zeros(12,12);
s = [-1/sqrt(3) 1/sqrt(3)];
t = [-1/sqrt(3) 1/sqrt(3)];

for i=1:2
    for j = 1:2
        k = k + (B(s(i),t(j)))'*D*B(s(i),t(j))*det(J);
        k = double(k);
    end
end        

k = w*k;
end