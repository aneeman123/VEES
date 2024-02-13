function [q,diff] = matquality(Approx, Ideal)
% function [q,diff] = matquality(Approx, Ideal)
% Uses 1-norm (max column) to measure relative error
% of Approx vs. Ideal in units of eps.
% I.e., q = 1 is great, bigger is worse.
% q = norm(Approx - Ideal, 1) / (norm(Ideal, 1) * eps);

diff = Approx - Ideal;
q = norm(diff, 1) / (norm(Ideal, 1) * eps);
return;
