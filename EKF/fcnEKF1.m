% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [cam, ekf, xhat, P, residual] = fcnEKF1(cam, ekf, index, xhat, P, dt, z, R, f)
%[x y z, dx dy dz, r p y, dr dp dy, xppo yppo focalLength,  tp1x tp1y tp1z, ...  tpnz   tpny   tpnz]', (15+3n x 1)
% 1 2 3  4   5  6  7 8 9  10 11 12    13   14          15     16   17   18, ... 10+3n  11+3n  12+3n]
nx = numel(xhat);  
Ix = eye(nx);
m2km = 1/1000;

%0. Initialize State xhat and Covariance Matrix P
ekf = fcnT2PPhiQ(cam,ekf,dt);

%1. Define Measurement Noise Matrix R

%2. Define Process Noise Matrix Q
om = speye(nx,nx);
Q = om * ((0*m2km)^2 * dt); %tie point transition uncertainty (m)
Q(1:12,1:12) = ekf.Q12; %process covariance
%Q(index.K,index.K) = eye(3)*(1E-16)^2 * dt;
% %3. Define State Transition Matrix Phi
Phi = om;  Phi(1:12,1:12) = ekf.Phi12;

%4. Propagate State xhat
xhat = Phi*xhat;

%5. Propagate Covariance Matrix P
P = Phi*P*Phi' + Q;

%R = ekf.zsigma^2;
if nargin>6 %if measurements present, update.
    for updateiterations = 1 %iterations
        %6. Get Estimated Measurement zhat for time t, zhat = H*xhat
        cam = fcnx2cam(cam,xhat,index,f);
        xy = ned2pixel(cam,f,xhat(index.tp),'aposteriori');  zhat=xy';  %xy = ned2pixel(cam,oned,tpned);
        
        %7. Get Linearized Measurement Matrix H for time t, %   zhat = H*xhat, (2x1 = 2x6*6x1) (nz x 1) = (nz x nx)*(nx x 1)
        H = fcnH(cam,xhat,index,f);
        
        %8. Compute Residual (Error = True Measurements - Estimated Measurements)
        residual = z - zhat(:); %pixels
        
        %9. Compute Measurement Prediction Covariance S
        PH = P*H';
        S = H*PH + R;
        
        %10. Compute Kalman Gain Matrix K
        K = PH/S;
        
        %11. Update xhat
        xhat = xhat + K*residual;
        
        %12. Update Covariance Matrix P
        IKH = Ix - K*H;
        P = IKH*P;
        %P = IKH*P*IKH' + K*R*K'; %less sensitive to round off errors
    end
end

%Update aposteriori estimate for frame i
cam = fcnx2cam(cam,xhat,index,f);

%Update aposteriori estimate for frame i+1 (for help with next frame's tie point matching/propagation)
if f < cam.frames;  cam = fcnx2cam(cam,Phi*xhat,index,f+1);  end
end

function cam = fcnx2cam(cam,x,index,f)
cam.aposteriori.ned(f,:) = cam.apriori.ned(f,:) + x(index.ned)'; %m
cam.aposteriori.rpy(f,:) = cam.apriori.rpy(f,:) + x(index.rpy)'; %deg
cam.xppo = x(index.K(1));  cam.yppo = x(index.K(2));  cam.focalLength = x(index.K(3));
end


function H = fcnH(cam,xhat,index,i)
d2r = pi/180;
nx = numel(xhat);
ntp = size(index.tp,1);
nz = ntp*2;
H = zeros(nz,nx);

tpnedhat = xhat(index.tp);
cfl=xhat(index.K(3));

o = cam.aposteriori.ned(i,:);  ox = o(1);  oy = o(2);  oz = o(3);
rpy = cam.aposteriori.rpy(i,:)*d2r;
roll=rpy(1);   cr=cos(roll);   sr=sin(roll);
pitch=rpy(2);  cp=cos(pitch);  sp=sin(pitch);
yaw=rpy(3);    cy=cos(yaw);    sy=sin(yaw);
for j=1:ntp
    k=(1:2)+(j-1)*2;
    dx=ox-tpnedhat(j,1);  dy=oy-tpnedhat(j,2);  dz=oz-tpnedhat(j,3);  k1 = (cp*cy*dx - sp*dz + cp*sy*dy)^2;
    
    dpxdtpxyz = [   (cfl*(cr*sy - cy*sp*sr))/(cp*cy*dx - sp*dz + cp*sy*dy) + (cfl*cp*cy*((cr*cy + sp*sr*sy)*dy - (cr*sy - cy*sp*sr)*dx + cp*sr*dz))/k1, (cfl*cp*sy*((cr*cy + sp*sr*sy)*dy - (cr*sy - cy*sp*sr)*dx + cp*sr*dz))/k1 - (cfl*(cr*cy + sp*sr*sy))/(cp*cy*dx - sp*dz + cp*sy*dy), - (cfl*cp*sr)/(cp*cy*dx - sp*dz + cp*sy*dy) - (cfl*sp*((cr*cy + sp*sr*sy)*dy - (cr*sy - cy*sp*sr)*dx + cp*sr*dz))/k1
        (cfl*cp*cy*((sr*sy + cr*cy*sp)*dx - (cy*sr - cr*sp*sy)*dy + cp*cr*dz))/k1 - (cfl*(sr*sy + cr*cy*sp))/(cp*cy*dx - sp*dz + cp*sy*dy), (cfl*(cy*sr - cr*sp*sy))/(cp*cy*dx - sp*dz + cp*sy*dy) + (cfl*cp*sy*((sr*sy + cr*cy*sp)*dx - (cy*sr - cr*sp*sy)*dy + cp*cr*dz))/k1, - (cfl*cp*cr)/(cp*cy*dx - sp*dz + cp*sy*dy) - (cfl*sp*((sr*sy + cr*cy*sp)*dx - (cy*sr - cr*sp*sy)*dy + cp*cr*dz))/k1];
    dpxdrpy = d2r*[    (cfl*((sr*sy + cr*cy*sp)*dx - (cy*sr - cr*sp*sy)*dy + cp*cr*dz))/(cp*cy*dx - sp*dz + cp*sy*dy), (cfl*(cp*cy*sr*dx - sp*sr*dz + cp*sr*sy*dy))/(cp*cy*dx - sp*dz + cp*sy*dy) + (cfl*((cr*cy + sp*sr*sy)*dy - (cr*sy - cy*sp*sr)*dx + cp*sr*dz)*(cp*dz + cy*sp*dx + sp*sy*dy))/k1, - (cfl*((cr*cy + sp*sr*sy)*dx + (cr*sy - cy*sp*sr)*dy))/(cp*cy*dx - sp*dz + cp*sy*dy) - (cfl*(cp*cy*dy - cp*sy*dx)*((cr*cy + sp*sr*sy)*dy - (cr*sy - cy*sp*sr)*dx + cp*sr*dz))/k1
        -(cfl*((cr*cy + sp*sr*sy)*dy - (cr*sy - cy*sp*sr)*dx + cp*sr*dz))/(cp*cy*dx - sp*dz + cp*sy*dy), (cfl*(cp*cr*cy*dx - cr*sp*dz + cp*cr*sy*dy))/(cp*cy*dx - sp*dz + cp*sy*dy) + (cfl*((sr*sy + cr*cy*sp)*dx - (cy*sr - cr*sp*sy)*dy + cp*cr*dz)*(cp*dz + cy*sp*dx + sp*sy*dy))/k1,   (cfl*((cy*sr - cr*sp*sy)*dx + (sr*sy + cr*cy*sp)*dy))/(cp*cy*dx - sp*dz + cp*sy*dy) - (cfl*(cp*cy*dy - cp*sy*dx)*((sr*sy + cr*cy*sp)*dx - (cy*sr - cr*sp*sy)*dy + cp*cr*dz))/k1];
    dpxdK = [ 1,  0,  ((cr*cy + sp*sr*sy)*dy - (cr*sy - cy*sp*sr)*dx + cp*sr*dz)/(cp*cy*dx - sp*dz + cp*sy*dy)
        0,  1,  ((sr*sy + cr*cy*sp)*dx - (cy*sr - cr*sp*sy)*dy + cp*cr*dz)/(cp*cy*dx - sp*dz + cp*sy*dy)];
    
    H(k,index.tp(j,:)) = dpxdtpxyz;
    H(k,index.ned) = -dpxdtpxyz;  %camera positions
    H(k,index.rpy) = dpxdrpy;
    H(k,index.K)= dpxdK;  %K = [xppo, yppo, focalLength];
end

H = sparse(H);
end