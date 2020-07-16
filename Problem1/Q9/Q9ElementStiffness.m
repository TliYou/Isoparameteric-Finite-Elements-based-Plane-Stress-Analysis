function k = Q9ElementStiffness(E,NU,w,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,x8,y8,x9,y9)


syms s t;

N1 = (1-s)*  (1-t)/4;
N2 = (1+s)*  (1-t)/4;
N3 = (1+s)*  (1+t)/4;
N4 = (1-s)*  (1+t)/4;
N5 = (1-s^2)*(1-t)/2;
N6 = (1+s)*  (1-t^2)/2;
N7 = (1-s^2)*(1+t)/2;
N8 = (1-s)*  (1-t^2)/2;
N9 = (1-s^2)*(1-t^2);

%Modification
N5 = N5-N9/2;
N6 = N6-N9/2;
N7 = N7-N9/2;
N8 = N8-N9/2;             %N5 to N8 Need to be modified first
N1 = N1-N5/2-N8/2-N9/4;
N2 = N2-N5/2-N6/2-N9/4;
N3 = N3-N6/2-N7/2-N9/4;
N4 = N4-N7/2-N8/2-N9/4;


% Partial Ni with respect to s or t
N1s = diff(N1,s); N2s = diff(N2,s); N3s = diff(N3,s); N4s = diff(N4,s); 
N5s = diff(N5,s); N6s = diff(N6,s); N7s = diff(N7,s); N8s = diff(N8,s);
N9s = diff(N9,s);

N1t = diff(N1,t); N2t = diff(N2,t); N3t = diff(N3,t); N4t = diff(N4,t);
N5t = diff(N5,t); N6t = diff(N6,t); N7t = diff(N7,t); N8t = diff(N8,t);
N9t = diff(N9,t);

J_f = [N1s N2s N3s N4s N5s N6s N7s N8s N9s;
       N1t N2t N3t N4t N5t N6t N7t N8t N9t];
 

J = J_f*[x1 y1; x2 y2; x3 y3; x4 y4; x5 y5; x6 y6; x7 y7; x8 y8; x9 y9];

Jinv = inv(J);
a = Jinv(1,1); b = Jinv(1,2); c = Jinv(2,1); d = Jinv(2,2);

B = [];
for p=1:9
    B_i = [a*J_f(1,p)+b*J_f(2,p)                      0;
        0                       c*J_f(1,p)+d*J_f(2,p);
        c*J_f(1,p)+d*J_f(2,p)   a*J_f(1,p)+b*J_f(2,p)];
    B = [B B_i];
end

B(s,t) = B;
J(s,t) = J;

D = (E/(1-NU*NU))*[1, NU, 0 ; NU, 1, 0 ; 0, 0, (1-NU)/2];


z = threebythreeGauss(B, D, J,w);     %Element stiffness matrix

k = z;