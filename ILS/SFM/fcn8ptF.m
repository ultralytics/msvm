function F = fcn8ptF(X1, X2)
% eightpoint - 8 Point algorithm for estimating the Fundamental matrix
% Computes the fundamental matrix F given a set of at least 8 image point 
% correspondaces following the 8-point algorithm described by Zisserman in p282.
% Input  - X1 -> (3xn) set of homogeneous points in image A
%        - X2 -> (3xn) set of homogeneous points in image B
% Output - F  -> (3x3) fundamental matrix
np = size(X1,2);

% Normalize the points following the centroid normalization recommended by Zisserman (p109, p282)
[X1t,T1] = normalize2Dpoints(X1,np);
[X2t,T2] = normalize2Dpoints(X2,np);

% Arrange the matrix A. Only the first two equations are needed as the last one is redundant
A= [X2t(1,:).*X1t(1,:)
    X2t(1,:).*X1t(2,:)
    X2t(1,:)
    X2t(2,:).*X1t(1,:)
    X2t(2,:).*X1t(2,:)
    X2t(2,:)
    X1t(1,:)
    X1t(2,:)
    ones(1,np)]';

%% Find the fundamental matrix F' in 2 steps
% Compute Â. Â is closest to A in frobenius form
%[U,D,V] = svd(A);
%D(end,end) = 0;
%Af = U*D*V';

% Linear solution: Determine F from the singular vector corresponding to the smaller singular value of Â
%[~,~, Va] = svd(Af);
[~,~, Va] = svd(A);

% The matrix F' is composed of the elements of the last vector of V
f = Va(:,9);

% Reorganice h to obtain H
Fp = reshape(f,3,3)';

% Constraint enforcement: Replace F by F' such that det(F´)=0 using SVD
[Uf,Df,Vf] = svd(Fp);
m = mean([Df(1,1),Df(2,2)]);
Df(1,1) = m;
Df(2,2) = m;
Df(3,3) = 0;
Ff = Uf*Df*Vf';
% Denormalization
F = T2'*Ff*T1;
end


function [X, T] = normalize2Dpoints(X,np)
% normalize2Dpoints - Normalization of 2D image points for the DLT algorithm
% Normalizes the points so that the origin is at the centroind of the
% points and the average distance to the origin is sqr(2)
% Input  - X  -> (2xn) matrix of image points (homogeneous points 3xn also work)
% Output - X -> (3xn) matrix of normalized points (including w=1)
%        - T  -> (3x3) Transformation matrix
X1 = X(1,:);  
X2 = X(2,:);

% Calculate centroid
c1=sum(X1)/np;  c2=sum(X2)/np;

% Adjust points so that the origin is at the centroid
Xn1 = X1-c1;
Xn2 = X1-c2;

% Compute the mean distance
r = sqrt(Xn1.^2+Xn2.^2);
meandistance = sum(r)/np;

% Compute the scale factor
scale = sqrt(2)/meandistance;

% Transformation matrix
T = [scale     0     -scale*c1
    0      scale     -scale*c2
    0          0             1 ];

% Recalculate points
X = T*X;
end