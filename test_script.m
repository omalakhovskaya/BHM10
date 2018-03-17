
A = randn(5,5);
A(1:2,1:2)=eye(2);
A(1:2,4:5)=zeros(2,2);
A(3,1)=1;
A(1:4,5)=zeros(4,1);
A(5,5)=1;
A


D = diag(abs(randn(5,1)))

[A_tilde, D_tilde, d_star, G]= forward(A,D)

[A_new,D_new] = backward(A_tilde, D_tilde, G)

A_new - A
D_new - D