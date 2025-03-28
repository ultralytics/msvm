% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function ned = pixel2neduvec(cam,i,xy,flag)
switch flag
    case 'true'
        rpy = cam.true.rpy(i,:);
    case 'apriori'
        rpy = cam.apriori.rpy(i,:);
    case 'aposteriori'
        rpy = cam.aposteriori.rpy(i,:);
end
C_CAM2NED = fcnRPY2DCM_B2W(rpy*d2r);

sc_cam = pixel2camsc(cam, double(xy));  cc_cam = fcnSC2CCd(sc_cam);
%cc_cam = pixel2camcc(cam,double(xy));
ned = cc_cam*C_CAM2NED'; %uvecs


end

