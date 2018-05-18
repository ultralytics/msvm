function K = fcnK(nx,ny,xppo,yppo,f)
%Fundamental Matrix F to Essential Matrix E
%http://www.umiacs.umd.edu/~ramani/cmsc828d/lecture27.pdf
%nx = number of horizontal pixels
%ny = number of vertical pixels
%fovh = horizontal fov (deg)
%fovv = vertical fov (deg)
d2r = pi/180;
 
%fx = (nx/2)/tan(fovh/2*d2r); 
%fy = (ny/2)/tan(fovv/2*d2r); 
 
kx = 1; %x scaling factor
ky = 1; %y scaling factor
 
fx =  f*kx; %x focal length (pixels)
fy =  -f*ky; %y focal length (pixels)
 
x0 = xppo+(nx+1)/2; %image center x pixel
y0 = yppo+(ny+1)/2; %image center y pixel
 
s = 0; %skew parameter
 
K = [fx   s    x0
     0   fy    y0
     0    0     1];  %Calibration Matrix
end



% function K = fcnK(nx,ny,f,fovh,fovv)
% %Fundamental Matrix F to Essential Matrix E
% %http://www.umiacs.umd.edu/~ramani/cmsc828d/lecture27.pdf
% %nx = number of horizontal pixels
% %ny = number of vertical pixels
% %fovh = horizontal fov (deg)
% %fovv = vertical fov (deg)
% d2r = pi/180;
%  
% fx = (nx/2)/tan(fovh/2*d2r); 
% %fy = (ny/2)/tan(fovv/2*d2r); 
%  
% 
% kx = 1; %x scaling factor
% ky = 1; %y scaling factor
%  
% fx =  f*kx; %x focal length (pixels)
% fy =  -f*ky; %y focal length (pixels)
%  
% x0 = (nx+1)/2; %image center x pixel
% y0 = (ny+1)/2; %image center y pixel
%  
% s = 0; %skew parameter
%  
% K = [fx   s    x0
%      0   fy    y0
%      0    0     1];  %Calibration Matrix
% end
