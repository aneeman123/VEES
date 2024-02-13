function [P] = projection(V)
% function [P] = projection(V) V is column vector, P is a matrix (alisa)
sz = size(V,1)
% size test sz
P = eye(size(V,1)) - ( (V*(V')) / dot(V,V) );

return;

