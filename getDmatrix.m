DD = function getDmatrix(mu,yyy1,yyy2, kappa, omega, S,AA,Pinv,m,Xtilde,Ytilde, n, T1, T2) 

Ytilde = [sqrt(mu)*yyy1*AA';yyy2*AA';Pinv'*m'];
         mstar = (Xtilde'*Xtilde)\(Xtilde'*Ytilde);
         omega=AA*S*AA';
         
         [~, tau_mh] = gettau3(kappa,omega,Ytilde,Xtilde,n);
         d= zeros(n,1);
         kappastar = kappa+(mu*T1+T2)/2;
         for jj = 1:n
             d(jj) = inv(gamrnd(kappastar,1/tau_mh(jj)));
         end 
         DD = diag(d);