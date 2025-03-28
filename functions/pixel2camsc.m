% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [sc, cc] = pixel2camsc(cam, xy)
xppo = cam.xppo + (cam.width+1)/2; %x principal point offset
yppo = cam.yppo + (cam.height+1)/2; %y principal point offset

xn = cam.focalLength;
yn = xy(:,1) - xppo;
zn = xy(:,2) - yppo;
ir = 1./sqrt(xn^2 + yn.^2 + zn.^2); %inverse range

[nr, ~] = size(xy);  sc = ones(nr,3);  cc = sc;

%sc = fcnCC2SC([xn yn zn]);
%sc(:,1) = 1;
sc(:,2) = asin(-zn.*ir) * (180/pi);
sc(:,3) = atan(yn/xn) * (180/pi);

cc(:,1) = xn.*ir;
cc(:,2) = yn.*ir;
cc(:,3) = zn.*ir;
end
