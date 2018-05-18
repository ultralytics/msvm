function rpy = solveAzPt_givenRoll(cam,a,fi) %fi = frame index
tpi=1;
x = cam.aposteriori.rpy(fi,:)'*(pi/180);  roll = x(1);
%sc = pixel2camsc(cam, fz);
%find el az, return DCM
%fhat1 = ned2pixel(cam,fi,xned,'aposteriori');
%C = fcnRPY2DCM_B2W(rpy);
x = x(2:3);

fz = pixel2camsc(cam, [a.upx(tpi,fi) a.upy(tpi,fi)]); fz=fz(:,2:3)'*(pi/180);
dx_ned = cam.tpnedhat(tpi,:)-cam.aposteriori.ned(fi,:);  TAx=dx_ned(:,1); TAy=dx_ned(:,2); TAz=dx_ned(:,3);
x = fcnelaz(dx_ned)'; %updated!!
rcc = sqrt(sum(dx_ned.^2));
diff1 = 9999;
iter = 0;
while diff1>1E-16 && iter<30
    iter = iter+1;
    r = roll;  sr=sin(r); cr=cos(r);
    p = x(1);  sp=sin(p); cp=cos(p);
    y = x(2);  sy=sin(y); cy=cos(y);
    k1 = (TAx*(sr*sy + cr*cy*sp) - TAy*(cy*sr - cr*sp*sy) + TAz*cp*cr);
    k2 = (TAy*(cr*cy + sp*sr*sy) - TAx*(cr*sy - cy*sp*sr) + TAz*cp*sr);
    k3 = (TAx*cp*cy - TAz*sp + TAy*cp*sy);
    k4 = (k3.^2 + k1.^2 + k2.^2);
        
    %C_NED2CAM = fcnRPY2DCM_B2W([roll x']); %W2B but need transpose!
    %cc = dx_ned*C_NED2CAM;
    %cc = fcnrotateB2W(roll,p,y,dx_ned);
    cc = [dx_ned(:,1).*(cp.*cy)+dx_ned(:,2).*(cp.*sy )+dx_ned(:,3).*(-sp),   dx_ned(:,1).*(sr.*sp.*cy-cr.*sy)+dx_ned(:,2).*(sr.*sp.*sy+cr.*cy)+dx_ned(:,3).*(sr.*cp),   dx_ned(:,1).*(cr.*sp.*cy+sr.*sy)+dx_ned(:,2).*(cr.*sp.*sy-sr.*cy)+dx_ned(:,3).*(cr.*cp)];
    %rcc = sqrt(cc(:,1).^2 + cc(:,2).^2 + cc(:,3).^2);
    fhat1 = [asin(-cc(:,3)./rcc);  atan2(cc(:,2),cc(:,1))]; %elaz
    
    dazdpy = [((TAx*cp*cy*sr - TAz*sp*sr + TAy*cp*sr*sy)./k3 + (k2.*(TAz*cp + TAx*cy*sp + TAy*sp*sy))./k3.^2)./(k2.^2./k3.^2 + 1), -((TAx*(cr*cy + sp*sr*sy) + TAy*(cr*sy - cy*sp*sr))./k3 + ((TAy*cp*cy - TAx*cp*sy).*k2)./k3.^2)./(k2.^2./k3.^2 + 1)];
    deldpy = [-((TAx*cp*cr*cy - TAz*cr*sp + TAy*cp*cr*sy)./k4.^.5 - (k1.*(2*k2.*(TAx*cp*cy*sr - TAz*sp*sr + TAy*cp*sr*sy) - 2*k3.*(TAz*cp + TAx*cy*sp + TAy*sp*sy) + 2*k1.*(TAx*cp*cr*cy - TAz*cr*sp + TAy*cp*cr*sy)))./(2*k4.^(3./2)))./(1 - k1.^2./k4).^.5, -((TAx*(cy*sr - cr*sp*sy) + TAy*(sr*sy + cr*cy*sp))./k4.^.5 - (k1.*(2*(TAx*(cy*sr - cr*sp*sy) + TAy*(sr*sy + cr*cy*sp)).*k1 + 2*(TAy*cp*cy - TAx*cp*sy).*k3 - 2*(TAx*(cr*cy + sp*sr*sy) + TAy*(cr*sy - cy*sp*sr)).*k2))./(2*k4.^(3./2)))./(1 - k1.^2./k4).^.5];
    B = [deldpy
        dazdpy];
    
%     dx = 1E-8;
%     sc = fcnCC2SCr(dx_ned*fcnRPY2DCM_B2W(x+[0 0 0]'));   fhat = sc(:,2:3); f0=fhat(:);
%     sc = fcnCC2SCr(dx_ned*fcnRPY2DCM_B2W(x+[dx 0 0]'));  fhat = sc(:,2:3); f1=fhat(:);
%     sc = fcnCC2SCr(dx_ned*fcnRPY2DCM_B2W(x+[0 dx 0]'));  fhat = sc(:,2:3); f2=fhat(:);
%     sc = fcnCC2SCr(dx_ned*fcnRPY2DCM_B2W(x+[0 0 dx]'));  fhat = sc(:,2:3); f3=fhat(:);
%     B = [f1-f0 f2-f0 f3-f0]/dx;
    
    f = fz - fhat1;
    Bt = B';
    N = Bt*B;
    x = x+N\(Bt*f);
    diff1 = sum(abs(f))/2;%
end
rpy = [roll x'];
%C_NED2CAM = fcnRPY2DCM_B2W(rpy);
%dM = M - Mtrue % comes out as zeros to prove that it works
end

function x2 = fcnrotateB2W(r,p,y,x)
sr=sin(r); sp=sin(p); sy=sin(y);
cr=cos(r); cp=cos(p); cy=cos(y);
x2 = [x(:,1).*(cp.*cy)+x(:,2).*(cp.*sy )+x(:,3).*(-sp),   x(:,1).*(sr.*sp.*cy-cr.*sy)+x(:,2).*(sr.*sp.*sy+cr.*cy)+x(:,3).*(sr.*cp),   x(:,1).*(cr.*sp.*cy+sr.*sy)+x(:,2).*(cr.*sp.*sy-sr.*cy)+x(:,3).*(cr.*cp)];
end

function x2 = fcnrotateW2B(r,p,y,x)
sr=sin(r);  sp=sin(p);  sy=sin(y);
cr=cos(r);  cp=cos(p);  cy=cos(y);
x2 = zeros(size(x));
x2(:,1)=x(:,1).*(cp.*cy)+x(:,2).*(sr.*sp.*cy-cr.*sy)+x(:,3).*(cr.*sp.*cy+sr.*sy);
x2(:,2)=x(:,1).*(cp.*sy)+x(:,2).*(sr.*sp.*sy+cr.*cy)+x(:,3).*(cr.*sp.*sy-sr.*cy);
x2(:,3)=x(:,1).*(-sp)+x(:,2).*(sr.*cp)+x(:,3).*(cr.*cp);
end

% function [M, Mtrue,Mned] = solveAzPt_givenRoll(cam,a,fi)
% myXpix = a.upx(1,fi);
% myYpix = a.upy(1,fi);
% 
% fov = 5*pi/180; % FOV in one direction
% d = 1280; % number of pixels in array in one direction
% 
% fl = d/(2*tan(fov/2));
% 
% % generate synthetic data
% XL = fcnNED2ENU*cam.aposteriori.ned(fi,:)'; % sensor position
% X = fcnNED2ENU*cam.true.focus.ned(fi,:)'; % aimpoint on ground associated with center of image
% up = [0; 0; 1]; % up vector
% zu = (XL-X)/norm(XL-X);
% xu = cross(up,zu)/norm(cross(up,zu));
% yu = cross(zu,xu);
% Mnr = [xu yu zu]'; % rotation matrix with no roll (nr)
% r = cam.aposteriori.rpy(fi,1)*pi/180; % true roll angle
% ca = cos(r);
% sa = sin(r);
% Mr = [ca sa 0
%     -sa ca 0
%     0 0 1];
% Mtrue = Mr*Mnr; % true direction cosine matrix (DCM)
% 
% % run image-to-ground for an arbitrary tie point p
% Xp = zeros(3,1); % allocate a matrix to fill later
% xp = myXpix - cam.width/2;
% yp = cam.height/2 - myYpix;
% Xp(3) = -cam.tpnedhat(1,3); 
% %Xp(3) = -a.ipned(1,3);
% M = Mtrue;
% uvw = M'*[xp; yp; -fl];
% Xp(1) = XL(1) + (Xp(3)-XL(3))*uvw(1)/uvw(3);
% Xp(2) = XL(2) + (Xp(3)-XL(3))*uvw(2)/uvw(3);
% Xp = fcnNED2ENU*cam.tpnedhat(1,:)';
% 
% % finished generating the synthetic data
% % following is the algorithm that you need to apply
% % given XL, YL, ZL, Xp, Yp, Zp, and the true roll angle, compute the DCM
% % first, compute the rotation matrix with no roll
% zu = (XL-Xp)/norm(XL-Xp);
% xu = cross(up,zu)/norm(cross(up,zu));
% yu = cross(zu,xu);
% Mnr = [xu yu zu]'; % rotation matrix with no roll (nr)
% % extract azimuth and pitch angles
% a = atan2(Mnr(1,2),Mnr(1,1));
% b = atan(Mnr(2,3)/Mnr(3,3));
% M = Mr*Mnr; % initial guess at DCM
% 
% % set up Newton iteration to refine values of a and b as a function of the
% % collinearity equations.  a and b are used to construct the updated DCM
% sum_abs_deltas = 999;
% iter = 0;
% while sum_abs_deltas>1e-12 && iter<30
%     iter = iter + 1;
%     UVW = M*(Xp-XL);
%     U = UVW(1);
%     V = UVW(2);
%     W = UVW(3);
%     f(1,1) = -(xp + fl*U/W);
%     f(2,1) = -(yp + fl*V/W);
%     if iter == 1
%         sa = sin(a);
%         ca = cos(a);
%         Ma = [ca sa 0
%             -sa ca 0
%             0 0 1];
%         sb = sin(b);
%         cb = cos(b);
%         Mb = [1 0 0
%             0 cb sb
%             0 -sb cb];
%     end
%     dMda = Mr*Mb*[-sa ca 0; -ca -sa 0; 0 0 0];
%     dUVW = dMda*(Xp-XL);
%     dU = dUVW(1);
%     dV = dUVW(2);
%     dW = dUVW(3);
%     B(1,1) = fl*(W*dU-U*dW)/W^2;
%     B(2,1) = fl*(W*dV-V*dW)/W^2;
%     dMdb = Mr*[0 0 0; 0 -sb cb; 0 -cb -sb]*Ma;
%     dUVW = dMdb*(Xp-XL);
%     dU = dUVW(1);
%     dV = dUVW(2);
%     dW = dUVW(3);
%     B(1,2) = fl*(W*dU-U*dW)/W^2;
%     B(2,2) = fl*(W*dV-V*dW)/W^2;
%     delta = B\f;
%     a = a + delta(1);
%     b = b + delta(2);
%     sum_abs_deltas = abs(delta(1)) + abs(delta(2));
%     sa = sin(a);
%     ca = cos(a);
%     Ma = [ca sa 0
%         -sa ca 0
%         0 0 1];
%     sb = sin(b);
%     cb = cos(b);
%     Mb = [1 0 0
%         0 cb sb
%         0 -sb cb];
%     M = Mr*Mb*Ma;
% end
% %M
% %dM = M - Mtrue % comes out as zeros to prove that it works
% Mned = [-M(3,2) M(1,2) -M(2,2)
%         -M(3,1) M(1,1) -M(2,1)
%          M(3,3) M(1,3)  M(2,3)];
% end

