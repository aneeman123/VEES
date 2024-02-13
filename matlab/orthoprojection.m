function [T] = orthoprojection(V,W)
% function [T] = orthoprojection(V,W) V is column vector, 
% W is a column vector 
% T is the result of projecting V against W ;
% you would subtract this whole thing from V to get the actual projection
sz = size(V,1)
% size test sz

T =   (( (W') * V ) / ( (W')*W ))*W; 

return;

