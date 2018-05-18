function [xy] = lla2pixel(cam,i,lla,flag)
ecef = lla2ecef(lla);
ned = ecef2ned(cam.DEM,ecef);
xy = ned2pixel(cam,i,ned,flag);
end