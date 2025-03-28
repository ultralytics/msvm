% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function sc = ned2camsc(cam,i,ned,flag)
switch flag
    case 'true'
        orig = cam.true.ned(i,:);
        rpy = cam.true.rpy(i,:);
    case 'apriori'
        orig = cam.apriori.ned(i,:);
        rpy = cam.apriori.rpy(i,:);
    case 'aposteriori'
        orig = cam.aposteriori.ned(i,:);
        rpy = cam.aposteriori.rpy(i,:);
end

C_NED2CAM = fcnRPY2DCM_B2W(rpy*(pi/180)); %W2B but need transpose!
dx_ned = [ned(:,1)-orig(1), ned(:,2)-orig(2), ned(:,3)-orig(3)];
dx_cam = dx_ned*C_NED2CAM;
sc = fcnCC2SC(dx_cam);
end
