clc; clear all

syms roll pitch yaw dx dy dz tpx tpy tpz ox oy oz cam_focalLength cam_width cam_height cam_xppo cam_yppo real


dx = tpx - ox;
dy = tpy - oy;
dz = tpz - oz;

C_NED2CAM = fcnRPY2DCM_B2W([roll pitch yaw]); %W2B but need transpose!
dx_ned = [dx, dy, dz];
dx_cam = dx_ned*C_NED2CAM;

%sc = fcnCC2SC(dx_cam);
cc = dx_cam;
% sc(:,1) = sqrt(cc(:,1).^2 + cc(:,2).^2 + cc(:,3).^2); %faster than sum
% sc(:,2) = asin(-cc(:,3)./sc(:,1));
% sc(:,3) = atan(cc(:,2)/cc(:,1));

% xy = camcc2pixel(cam, dx_cam);
gain = cam_focalLength./cc(:,1);
px = cc(:,2).*gain + ((cam_width+1)/2 + cam_xppo);  % +x principal point offset
py = cc(:,3).*gain + ((cam_height+1)/2 + cam_yppo);  % +y principal point offset


dpxdrpy = [     diff(px,roll)  diff(px,pitch)  diff(px,yaw)
                diff(py,roll)  diff(py,pitch)  diff(py,yaw)     ];

dpxdtpxyz = [   diff(px,tpx)  diff(px,tpy)  diff(px,tpz)
                diff(py,tpx)  diff(py,tpy)  diff(py,tpz)        ];

dpxdK = [       diff(px,cam_xppo)  diff(px,cam_yppo)  diff(px,cam_focalLength) 
                diff(py,cam_xppo)  diff(py,cam_yppo)  diff(py,cam_focalLength)        ];



% ddpxdrpy = [    diff(diff(px,roll),roll)  diff(diff(px,pitch),pitch)  diff(diff(px,yaw),yaw)
%                 diff(diff(py,roll),roll)  diff(diff(py,pitch),pitch)  diff(diff(py,yaw),yaw)];
% 
% ddpxdtpxyz = [  diff(diff(px,tpx),tpx)  diff(diff(px,tpy),tpy)  diff(diff(px,tpz),tpz)
%                 diff(diff(py,tpx),tpx)  diff(diff(py,tpy),tpy)  diff(diff(py,tpz),tpz)];
            
            
            