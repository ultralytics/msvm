% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [xy] = ecef2pixel(cam,i,ecef,flag)
ned = ecef2ned(cam.DEM,ecef);
xy = ned2pixel(cam,i,ned,flag);
end


