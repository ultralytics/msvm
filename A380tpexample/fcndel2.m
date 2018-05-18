function v = fcndel2(f,varargin)
%DEL2 Discrete Laplacian.
%   L = DEL2(U), when U is a matrix, is a discrete approximation of
%   0.25*del^2 u = (d^2u/dx^2 + d^2u/dy^2)/4.  The matrix L is the same
%   size as U, with each element equal to the difference between an 
%   element of U and the average of its four neighbors.
%
%   L = DEL2(U), when U is an N-D array, returns an approximation of
%   (del^2 u)/2/n, where n is ndims(u).
%
%   L = DEL2(U,H), where H is a scalar, uses H as the spacing between
%   points in each direction (H=1 by default).
%
%   L = DEL2(U,HX,HY), when U is 2-D, uses the spacing specified by HX
%   and HY. If HX is a scalar, it gives the spacing between points in
%   the x-direction. If HX is a vector, it must be of length SIZE(U,2)
%   and specifies the x-coordinates of the points.  Similarly, if HY
%   is a scalar, it gives the spacing between points in the
%   y-direction. If HY is a vector, it must be of length SIZE(U,1) and
%   specifies the y-coordinates of the points.
%
%   L = DEL2(U,HX,HY,HZ,...), when U is N-D, uses the spacing given by
%   HX, HY, HZ, etc. 
%
%   Class support for input U:
%      float: double, single
%
%   See also GRADIENT, DIFF.

%   Copyright 1984-2011 The MathWorks, Inc. 
%   $Revision: 5.16.4.6 $  $Date: 2011/05/17 02:22:18 $

%[err,f,ndim,loc,cflag] = parse_inputs(f,varargin);

ndim = ndims(f);
if ndim == 1
  perm = [1 2];
else
  perm = [2:ndim 1]; % Cyclic permutation
end

v = zeros(size(f),class(f));
for k = 1:ndim
   [n,p] = size(f);
   g  = zeros(size(f),class(f)); % case of singleton dimension
   
   if isempty(varargin)
       h = 1; %spacing
   else
       h = varargin{1};
   end
   
   % Take centered second differences on interior points to compute g/2
   if n > 2
      hk3 = f(2:n-1,:);
      %g(2:n-1,:) = ( (f(3:n,:)-hk3)/h  - (hk3-f(1:n-2,:))/h) / (h+h);
      %g(2:n-1,:) = ( (f(3:n,:)-hk3)  - (hk3-f(1:n-2,:)) ) / (2*h^2);
      g(2:n-1,:) = ( f(3:n,:)-2*hk3+f(1:n-2,:) ) / (2*h^2);
   end

   % Linearly extrapolate second differences from interior
   if n > 3
      g(1,:) = g(2,:)*2 - g(3,:);
      g(n,:) = -g(n-2,:) + g(n-1,:)*2;
   elseif n==3
      g(1,:) = g(2,:);
      g(n,:) = g(2,:);
   else
      g(1,:)=0;
      g(n,:)=0;
   end

   if ndim==1,
     v = v + g;
   else
     v = v + ipermute(g,[k:ndim 1:k-1]);
   end

   % Set up for next pass through the loop
   f = permute(f,perm);
end 
v = v./ndims(f);


