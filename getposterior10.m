function [posterior] =getposterior10(theta,param,S,m,yyy1,yyy2,Pinv,Xtilde,mu,n,c,kappa,T1,T2,omega_tildeT)

A = [1 0 theta(1) 0 0; ...
    0 1 theta(2) 0 0; ...
    1 theta(3:4)' theta(5) 0;...
    theta(6:9)'  0;...
    theta(10:13)' 1]
D = diag([theta(14:18)])

[A_tilde, D_tilde,dstar,~] = forward_operator(A,D); 

x = [-A_tilde(1,3);-A_tilde(2,3);-A_tilde(3,2);-A_tilde(3,3);-1/A_tilde(3,4);...
    -A_tilde(4,1);-A_tilde(4,2);-A_tilde(4,3);-A_tilde(5,1);-A_tilde(5,2);...
    -A_tilde(5,3);-A_tilde(5,4)];

prior = getprior10(x,param,c);
log_priors = log(prod(prior));

kappastar=kappa+(mu*T1+T2)/2;

Ytilde = [sqrt(mu)*yyy1*A';yyy2*A';Pinv'*m'];
omega=A*S*A';
[~, taustar] = gettau3(kappa, omega, Ytilde, Xtilde,n);

up=log_priors+(mu*T1+T2)/2*log(det(A*omega_tildeT*A'))+ sum(kappa*log(kappa*diag(omega)));
down = kappastar*(n*(log(2) - log(mu*T1+T2))+ sum(log(taustar)));
posterior=up-down





