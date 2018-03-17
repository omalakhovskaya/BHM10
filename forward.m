function [A_tilde, D_tilde, d_star, G] = forward(A, D)
%FORWARD Function finds non-transformed matrices if transformed matrices
%are known
if size(A) ~= [5,5] | size(D) ~= [5,5]
display('there is a dimension errors, matrix A and D should have 5 rows and 5 columns') 
end
A_tilde = zeros(5,5);
D_tilde = zeros(5,5);
d_star = zeros(5,1);
G = eye(5);

rho = (A(4,4) - 1)/A(3,4);
sigma_e_sq = D(3,4)*A(3,4)-(A(4,4)-1)*D(3,3);
gamma = -A(5,4)*(D(3,3)*rho+D(3,3)*D(4,4)-A(3,4)*sigma_e_sq*D(4,4)-...
    D(3,3)*rho*A(3,4)*sigma_e_sq+rho*A(3,4)*sigma_e_sq^2-sigma_e_sq -...
    rho*sigma_e_sq-rho*A(3,4)*sigma_e_sq)^(-1);

D_tilde(1,1)= D(1,1);
D_tilde(2,2)= D(2,2);
D_tilde(3,3)= D(3,3);
D_tilde(3,4) = D(3,4) - rho*D(3,3);
D_tilde(4,3) = D_tilde(3,4);
D_tilde(4,4) = D(4,4) - rho*D_tilde(3,4) ;
D_tilde(3,5)= -A(3,4)*gamma*sigma_e_sq ;
D_tilde(5,3)= -A(3,4)*gamma*sigma_e_sq ;
D_tilde(4,5) = -gamma*sigma_e_sq;
D_tilde(5,4) = -gamma*sigma_e_sq;

tau = (D_tilde(3,3)*D_tilde(4,5) - D_tilde(3,4)*D_tilde(3,5))/...
    (D_tilde(3,4)^2 - D_tilde(3,3)*D_tilde(4,4));
phi = -(D_tilde(3,5)+tau*D_tilde(3,4))/D_tilde(3,3);

D_tilde(5,5) = D(5,5) - tau*(phi*D_tilde(3,4)+tau*D_tilde(4,4)+D_tilde(4,5))...
    -phi*D_tilde(3,5)-tau*D_tilde(4,5);

G(4,3)= rho; 
G(5,3) = phi;
G(5,4) = tau;

A_tilde = inv(G)*A;

d_star(1) = D(1,1);
d_star(2) = D(2,2);
%d_star(3) = D(3,3)-D(3,4)*A(3,4)^3+A(3,4)^2*(A(4,4)-1)*D(3,3);
d_star(3) = D(3,3) - A(3,4)^2*sigma_e_sq; 
d_star(4) = (D_tilde(4,4)-sigma_e_sq)*A(3,4)^2;
d_star(5) = D(5,5) - gamma^2*sigma_e_sq;








