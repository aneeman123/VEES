function [H6] = higham6(n)
% function [H6] = higham6(n)
% returns a 6x6 of rank 5 based on gallery(5), a matlab function
% attributed to Higham.  See END OF "help gallery".
% H6 = [1 zeros(1,5); zeros(5,1) gallery(5)];

H6 = [1 zeros(1,5); zeros(5,1) gallery(5)];
return;
