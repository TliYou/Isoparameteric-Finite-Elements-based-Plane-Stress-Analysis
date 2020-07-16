function stress = CSTElementStresses(E,NU,xi,yi,xj,yj,xm,ym,d)
%d is the element displacement vector
A = (xi*(yj-ym) + xj*(ym-yi) + xm*(yi-yj))/2;
beta_i = yj-ym;
beta_j = ym-yi;
beta_m = yi-yj;
gamma_i = xm-xj;
gamma_j = xi-xm;
gamma_m = xj-xi;

B = (1/(2*A))*[beta_i    0       beta_j   0     beta_m       0 ;
               0       gamma_i    0     gamma_j    0     gamma_m;
               gamma_i beta_i  gamma_j  beta_j  gamma_m   beta_m];

D = (E/(1-NU*NU))*[1 NU 0 ; NU 1 0 ; 0 0 (1-NU)/2];
stress = D*B*d;