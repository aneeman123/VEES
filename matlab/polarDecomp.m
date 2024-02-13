function [R,U] = polarDecomp(D)
% function [R,U] = polarDecop(D) D is a matrix , R orthogonal rotation, 
% U symmetric positive (semi?)definite stretch

DTD = (D')*D;
[T,J] = eig(DTD);

SQRTJ = sqrt(J);
%DTD = UTU = U^2
U = (T)*SQRTJ*(T');

% 

% C = RT, ultimately we will solve for R orthogonal
% calculate F, the fake inverse of sqrt(J)
% DTF = CF

F = SQRTJ;
% case of one zero eigenvalue - do I need to invert U?
if ( (F(1,1) <= 0.00001) && (F(1,1) >= 0.0) )
   F(1,1) = 1.0;
end;

F = inv(F);

G = D*T*F;
Gorig = G;
% find row with smallest magnitude in G for 1
X = [norm(G(1,:)) ; norm(G(2,:)) ; norm(G(3,:)) ; norm(G(4,:)) ; ... 
norm(G(5,:)) ; norm(G(6,:)) ] ;

[Rmin,Imin] = min(X);

G1 = [0; 0; 0; 0; 0; 0]; % column

G1(Imin) = 1.0;

P2 = G(1:6,2);

P3 = G(1:6,3) - orthoprojection(G(1:6,3),P2);

P4 = G(1:6,4) - orthoprojection(G(1:6,4),P3) - orthoprojection(G(1:6,4),P2);

P5 = G(1:6,5) - orthoprojection(G(1:6,5),P4) - orthoprojection(G(1:6,5),P3) - orthoprojection(G(1:6,5),P2);
 
P6 = G(1:6,6) - orthoprojection(G(1:6,6),P5) - orthoprojection(G(1:6,6),P4) - orthoprojection(G(1:6,6),P3) - orthoprojection(G(1:6,6),P2);

%P1 = G1 -  orthoprojection(G1,G(:,6)) -orthoprojection(G1,G(:,5)) ... 
%-orthoprojection(G1,G(:,4)) - orthoprojection(G1,G(:,3)) ... 
%-orthoprojection(G1,G(:,2));

P1 = projection(G(:,6)) *projection(G(:,5)) ...
*projection(G(:,4)) * projection(G(:,3)) ...
*projection(G(:,2))*G1;



P1 = P1/(norm(P1)); % normalize

R = [P1 P2 P3 P4 P5 P6];

R = R*(T');
