function ecef = pixel2ecef(DEM,cam,i,xy,flag)
switch flag
    case 'true'
        orig = cam.true.ecef(i,:);
        rpy = cam.true.rpy(i,:);
    case 'apriori'
        orig = cam.apriori.ecef(i,:);
        rpy = cam.apriori.rpy(i,:);
    case 'aposteriori'
        orig = cam.aposteriori.ecef(i,:);
        rpy = cam.aposteriori.rpy(i,:);
end

C_NED2ECEF = DEM.DCM_ECEF2NED'; %from DEM centerpoint
C_CAM2NED = fcnRPY2DCM_B2W(rpy*(pi/180));
C_CAM2ECEF = C_NED2ECEF*C_CAM2NED;

sc_cam = pixel2camsc(cam, double(xy));
cc_cam = fcnSC2CCd(sc_cam);
cc_ecef = cc_cam*C_CAM2ECEF'; %uvecs

dest = [cc_ecef(:,1)+orig(1), cc_ecef(:,2)+orig(2), cc_ecef(:,3)+orig(3)]; %target

ecef = fcnrayDEMintersection(DEM,orig,dest); %ecef
end

