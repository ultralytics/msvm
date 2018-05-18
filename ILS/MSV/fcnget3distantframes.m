function [n1, n2, n3] = fcnget3distantframes(xyz)
%selects the 2 farthest frames from n1
nf = size(xyz,1);
n1=1; %first frame
f = ceil(rand(1E4,2)*nf);  %s=rng;  rng(1);  f = ceil(rand(1E4,2)*nf);  rng(s);
ax = xyz(n1,:);
bx = xyz(f(:,1),:);
cx = xyz(f(:,2),:);
r = fcnrange(ax,bx).*fcnrange(ax,cx).*fcnrange(bx,cx);
[~, i] = max(r);
x = sort(f(i,:));
n2 = x(1);
n3 = x(2);
end


