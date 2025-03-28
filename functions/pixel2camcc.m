% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function cc = pixel2camcc(cam, xy)
x0 = cam.xppo + (cam.width+1)/2; %x principal point offset
y0 = cam.yppo + (cam.height+1)/2; %y principal point offset

xn = cam.focalLength;
yn = xy(:,1) - x0;
zn = xy(:,2) - y0;
ir = 1./sqrt(xn^2 + yn.^2 + zn.^2); %inverse range

[nr, ~] = size(xy);  cc = zeros(nr,3);
cc(:,1) = xn.*ir;
cc(:,2) = yn.*ir;
cc(:,3) = zn.*ir;
end
