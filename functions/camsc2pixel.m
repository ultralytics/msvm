% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function xy = camsc2pixel(cam, sc)
% el = sc(:,2);
% az = sc(:,3);
[nr,~] = size(sc);
cc = fcnSC2CCd([ones(nr,1) sc(:,2:3)]);
xy = camcc2pixel(cam, cc);
end
