% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function lla = pixel2lla(DEM,cam,i,xy,flag)
ecef = pixel2ecef(DEM,cam,i,xy,flag);
lla = ecef2lla(ecef);
end

