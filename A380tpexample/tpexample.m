cdata = importdata('A380.JPG');

I = single(mean(cdata,3));
[rows, cols] = size(I);
zc = zeros(1,cols);
zr = zeros(rows,1);

%ORIGINAL
%x = [zc; diff(I,1,1)];
%y = [diff(I,1,2) zr];
%d = abs(x).*abs(y);

%DEL2 PERFECT BUT SLOW :(
d = fcndel2(I,1);
figure; image(d)

%[dx, dy] = gradient(I);
%d = dx.*dy;

%SECOND DERIVATIVE MANUALLY!
x = [zc; diff(I,2,1); zc];
y = [zr diff(I,2,2) zr];
d = sqrt(x.^2 + y.^2);
d = d/max3(d) * 255;
figure; image(d)

% d=abs(d(:));
% v1 = find(d>10);
% [ds, i] = sort(d(v1,:),1,'descend');
% i = v1(i(1:300));
% [ii, jj] = ind2sub([rows cols],i);
% 
% figure; image(cdata); hold on; axis equal vis3d
% plot(jj,ii,'g+')

