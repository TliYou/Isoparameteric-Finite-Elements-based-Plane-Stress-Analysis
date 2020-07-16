function k = Q6ElementStiffness(E,NU,w,x1,y1,x2,y2,x3,y3,x4,y4)

syms s t;

J_f = [-(1-t)   1-t    1+t  -(1+t);
          -(1-s) -(1+s)  (1+s)   1-s];
J = (1/4)*J_f*[x1 y1; x2 y2; x3 y3; x4 y4];

Jinv = inv(J);
a = Jinv(1,1); b = Jinv(1,2); c = Jinv(2,1); d = Jinv(2,2);

B1 = [a*(t-1)/4+b*(s-1)/4                    0 ;
      0                    c*(t-1)/4+d*(s-1)/4 ;
      c*(t-1)/4+d*(s-1)/4  a*(t-1)/4+b*(s-1)/4 ];
  
B2 = [a*(1-t)/4+b*(-1-s)/4                      0;
      0                      c*(1-t)/4+d*(-1-s)/4;
      c*(1-t)/4+d*(-1-s)/4   a*(1-t)/4+b*(-1-s)/4];
  
B3 = [a*(t+1)/4+b*(s+1)/4                       0 ; 
      0                       c*(t+1)/4+d*(s+1)/4 ;
      c*(t+1)/4+d*(s+1)/4     a*(t+1)/4+b*(s+1)/4];
  
B4 = [a*(-1-t)/4+b*(1-s)/4                      0 ;
      0                      c*(-1-t)/4+d*(1-s)/4 ;
      c*(-1-t)/4+d*(1-s)/4    a*(-1-t)/4+b*(1-s)/4];
  
B5 = [a*(-2*s) b*(-2*t);
      0               0;
      c*(-2*s)  d*(-2*t)];
  
B6 = [0             0;
      c*(-2*s)  d*(-2*t);
      a*(-2*s) b*(-2*t)];


B(s,t) = [B1 B2 B3 B4 B5 B6];

D = (E/(1-NU*NU))*[1, NU, 0 ; NU, 1, 0 ; 0, 0, (1-NU)/2];


K_all = twobytwoGauss(B, D, J,w);     % Non-condensed Element stiffness matrix

% Here K_all = [k11, k12; k21, k22]
k11 = K_all(1:8,1:8);
k12 = K_all(1:8,8:12);
k21 = K_all(8:12,1:8);
k22 = K_all(8:12,8:12);

keq = k11 - k12*(k22)^-1*k21;

k = keq;