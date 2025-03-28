% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function xy = camcc2pixel(cam, cc)
% el = sc(:,2);
% az = sc(:,3);
gain = cam.focalLength./cc(:,1);
x0 = cam.xppo + (cam.width+1)/2;
y0 = cam.yppo + (cam.height+1)/2;

[nr,~] = size(cc);  xy = zeros(nr,2);
xy(:,1) = cc(:,2).*gain + x0;
xy(:,2) = cc(:,3).*gain + y0;


