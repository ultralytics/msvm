% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [cam,a] = fcncropimages(cam,a,i)

cam.apriori.lla = cam.apriori.lla(i,:);    cam.aposteriori.lla = cam.aposteriori.lla(i,:);    cam.true.lla = cam.true.lla(i,:);
cam.apriori.ecef = cam.apriori.ecef(i,:);  cam.aposteriori.ecef = cam.aposteriori.ecef(i,:);  cam.true.ecef = cam.true.ecef(i,:); 
cam.apriori.ned = cam.apriori.ned(i,:);    cam.aposteriori.ned = cam.aposteriori.ned(i,:);    cam.true.ned = cam.true.ned(i,:); 
cam.apriori.rpy = cam.apriori.rpy(i,:);    cam.aposteriori.rpy = cam.aposteriori.rpy(i,:);    cam.true.rpy = cam.true.rpy(i,:);     cam.true.t=cam.true.t(i);  cam.true.dt=[0 diff(cam.true.t)];

cam.frames = numel(i);   cam.frameID = cam.frameID(i);  

a.upx = a.upx(:,i);  
a.upy = a.upy(:,i);  
a.state = a.state(:,i);  
end

