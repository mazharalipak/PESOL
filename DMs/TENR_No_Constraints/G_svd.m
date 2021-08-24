function [ G, u, v] = G_svd(J)
%tr_svd Calculate the value of G and its derivatives
%   x - variable 
%   J - Jacobian function: evaluates J(x)

 [u, G, v] = svds(J, 1,'smallest');
end

