function [A, D] = backward(A_tilde, D_tilde, G)
%BACKWARD Function finds transformed matrices if non-transformed matrices
%are known
A = G*A_tilde; 
D = G*D_tilde*G';
end

