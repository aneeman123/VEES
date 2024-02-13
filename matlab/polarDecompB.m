function [R,U,T,J,SQRTJ,F,Gorig,G,C,zeroEigs] = polarDecompB(D)
% function [R,U,T,J,SQRTJ,F,Gorig,G,C,zeroEigs] = polarDecompB(D)
% D is a square matrix , R orthogonal rotation, 
% U symmetric positive (semi?)definite stretch

n = size(D,1);
sz2 = size(D,2);
if n ~= sz2 return; end;

DTD = (D')*D;
[T,J] = eig(DTD);

if J(n,n) == 0
    allZeroEigs = 1
    return;
    end;

zeroEigs = 0;
for j = 1:n
    if J(j,j) < 0
        J(j,j) = 0;
        zeroEigs = zeroEigs + 1;
    elseif J(j,j) < eps*J(n,n)
        J(j,j) = 0;
        zeroEigs = zeroEigs + 1;
    elseif J(j,j) < sqrt(eps)*J(n,n)
        zeroEigs = zeroEigs + 1;
        end;
    end;

SQRTJ = sqrt(J);
%DTD = UTU = U^2
U = (T)*SQRTJ*(T');

% Ensure U is symmetric
U = (U + U') / 2;


% C = RT, ultimately we will solve for R orthogonal
% calculate F, the fake inverse of sqrt(J)
% DTF = CF

F = SQRTJ;
% case of one or more zero eigenvalue - do I need to invert U?
for j = 1:zeroEigs
    F(j,j) = 1.0;
    end;

F = inv(F);

G = D*T*F;
Gorig = G;

% Tune up G.  All columns of G not corresponding to a zero eigenvalue
% should be orthogonal.
% All columns of G not corresponding to a zero eigenvalue
% should be unit length.  Force this to be true.


for j = n:-1: zeroEigs+1
    temp = G(:,j);
    for k = j+1:n
        temp = temp - orthoprojection(temp, G(:,k));
        end;
    G(:,j) = normc(temp);
    end;


%P1 = G1 -  orthoprojection(G1,G(:,6)) -orthoprojection(G1,G(:,5)) ... 
%-orthoprojection(G1,G(:,4)) - orthoprojection(G1,G(:,3)) ... 
%-orthoprojection(G1,G(:,2));

% All columns of G not corresponding to a zero eigenvalue
% should be unit length.  Force this to be true.
% for j = zeroEigs+1:n
%     G(:,j) = G(:,j) / norm(G(:,j));
%     end;
G(:, zeroEigs+1:n) = normc(G(:, zeroEigs+1:n));

% find row with smallest magnitude in G for 1
X = zeros(n,1);
for i = 1:n
    X(i) = norm(G(i,:));
    end;

% This only allows for zeroEigs = 0 or 1.
[Rmin,Imin] = min(X);
G1 = zeros(n,1);
G1(Imin) = 1.0;

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
    P1 = G1;
    for j = 2:n
        P1 = P1 - orthoprojection(P1, G(:,j));
        end;
    C1 = normc(P1);
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
