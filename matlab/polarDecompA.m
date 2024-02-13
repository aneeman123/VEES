function [R,U,T,J,SQRTJ,U0,F,Gorig,G,C,zeroEigs0] = polarDecompA(D, toler)
% function [R,U,T,J,SQRTJ,U0,F,Gorig,G,C,zeroEigs0] = polarDecompA(D, toler)
% D is a square matrix with det(D) >= 0 (not checked, see below).
% toler controls near-zero threshold.
%   Plausible values of toler range from 1 to 10000, but 0 and sqrt(1/eps)
%   are also worth experimenting with.
%   The actual cutoff for an eigenvalue of J = D^T D is toler * eps * J(n,n),
%   where J(n,n) is the maximum eigenvalue.
%
% Returns
%    R orthogonal rotation, with det(R) = +1.
%    U symmetric positive semidefinite stretch, D = R U is the polar decomp
%    T eigenvectors of U (actually of U0); column-1 is eigenvector of 0 for
%        D IF D IS SINGULAR.
%    Other diagnostic values. zeroEigs0 is how many eigenvalues were at
%    or below the cutoff.
%
% Adjustment for det(D) < 0:  Define T1 = T with column-1 negated.
% Define R1 = R * T * T1'.  Now det(R1) < 0 and D ~= R1 * U.
% In theory T1 can be T with any odd number of columns negated.
% Is it true that R1 comes out the same in all cases ????
% The answer is unclear, so this is not implemented, and we require
% det(D) >= 0 to guarantee a correct decomposition.

n = size(D,1);
sz2 = size(D,2);
if n ~= sz2 return; end;

DTD = (D')*D;
[T,J] = eig(DTD);

if J(n,n) == 0
    ALL_ZERO_EIGS = 'polarDecompA(): WARNING, ALL_ZERO_EIGS'
    R = eye(n);
    U = zeros(n);
    T = eye(n);
    J = zeros(n);
    SQRTJ = zeros(n);
    F = eye(n);
    Gorig = eye(n);
    G = eye(n);
    C = eye(n);
    zeroEigs0 = n;
    return;
    end;

% Ensure det(T) > 0 by reversing eigenvector of smallest eigenvalue,
% if necessary.
if det(T) < 0
    T(:,1) = -T(:,1);
    end;

zeroEigs = 0;
for j = 1:n
    if J(j,j) < 0
        J(j,j) = 0;
        end;
    if J(j,j) <= toler * eps * J(n,n)
        zeroEigs = zeroEigs + 1;
        end;
    end;

zeroEigs0 = zeroEigs;

% Later code cannot handle zeroEigs > 1.
if zeroEigs > 1 zeroEigs = 1; end;

SQRTJ = sqrt(J);
%DTD = UTU = U^2
U0 = (T)*SQRTJ*(T');

% Ensure U0 is symmetric
U0 = (U0 + U0') / 2;


% C = RT, ultimately we will solve for R orthogonal
% calculate F, the fake inverse of sqrt(J)
% DTF = CF

F = SQRTJ;
% case of one or more zero eigenvalue
for j = 1:zeroEigs
    if F(j,j) == 0 F(j,j) = 1.0; end;
    end;

F = inv(F);

G = D*T*F;
Gorig = G;

% Tune up G.  All columns of G not corresponding to a zero eigenvalue
% should be orthogonal.
% All columns of G not corresponding to a zero eigenvalue
% should be unit length.  Force this to be true.
%
for j = n : -1 : zeroEigs+1
    temp = G(:,j);
    for k = j+1:n
        temp = temp - orthoprojection(temp, G(:,k));
        end;
    G(:,j) = normc(temp);
    end;

% find row with smallest magnitude in G for 1
X = zeros(n,1);
for i = 1:n
    X(i) = norm(G(i,:));
    end;

% This only allows for zeroEigs = 0 or 1.
[Rmin,Imin] = min(X);
G1 = zeros(n,1);
G1(Imin) = 1.0;

% Ensure that det(final G) > 0.
if (det( [G1 G(:,2:n)] ) < 0)
    G1 = -G1;
    end;

method = 2;
if method == 1
    % Use projection matrices
    % Only coded for zeroEigs = 0 or 1, but generalizes.
    P1 = G1;
    for j = 2:n
        P1 = projection(G(:,j)) * P1;
        end;
    C1 = normc(P1);
elseif method == 2
    % Use Gram-Schmidt
    % Only coded for zeroEigs = 0 or 1, but generalizes.
    P1 = zeros(n,1);
    for j = 2:n
        P1 = P1 - orthoprojection(G1, G(:,j));
        end;
    C1 = normc(G1 + P1);
elseif method == 3
    % Use Cramer's rule instead (ignores G1).
    % This only works for zeroEigs = 0 or 1.
    % Does not generalize in any obvious way.
    C1 = zeros(n,1);
    I = eye(n);
    for i = 1:n
        C1(i) = det([I(:,i) G(:, 2:n)] );
        end;
    end;

C = [C1 G(:, 2:n)];

R = C*(T');

U2 = inv(R) * D;   % more accurate than assuming inv(R) = R'.
U = (U2 + U2') / 2;

return;
