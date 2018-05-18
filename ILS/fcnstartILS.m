function [ils, cam] = fcnstartILS(cam,ekf,a)
startclock1 = clock;
d2r = pi/180;  r2d=180/pi;
nframes = cam.frames;
cfl=cam.focalLength;

ntp = size(a.upx,1); %number of tie points
activetp = a.state==1;    vp2(1,:,:) = activetp;  vp2(2,:,:) = activetp;  %valid pixel indices for pixelresiduals
ntpe = sum3(activetp);%number of tie point equations
mtpl = mean(sum(activetp,2)); %mean tie point life

nz = ntpe*2 + 12*nframes; %number of measurments total
nx = ntp*3 + 12*nframes; %number of parameters
nf = ntpe*2 + 12*nframes; %number of equations
cam.apriori.ned = cam.apriori.ned(1:nframes,:);   cam.aposteriori.ned = cam.apriori.ned;  cam.true.ned = cam.true.ned(1:nframes,:);   
cam.apriori.rpy = cam.apriori.rpy(1:nframes,:);   cam.aposteriori.rpy = cam.apriori.rpy;  cam.true.rpy = cam.true.rpy(1:nframes,:);

% %RUN SFM INITIAL GUESS ----------------------------------------------------
% [cam, metric] = fcnSFM(cam,a,30);
% figure; plot(cam.apriori.rpy(:,2)); hold on; plot(cam.true.rpy(:,2),'g')
% fprintf('aposteriori = %.2fdeg, Metric=%.2fkm\n',fcnrms(unwrap(cam.apriori.rpy(1:nframes,:)*d2r-cam.true.rpy(1:nframes,:)*d2r)*r2d),metric)
% %RUN SFM INITIAL GUESS ----------------------------------------------------

% %RUN SFM INITIAL GUESS ----------------------------------------------------
% rpy = zeros(nframes,3);
% v1 = round(linspace(1,nframes,cam.sfmframes));
% for i=v1
%     [cam, metric] = fcnSFM(cam,a,i);
%     sc = cam.aposteriori.rpy*d2r; sc(:,1)=1;
%     rpy = rpy + fcnSC2CC(sc);  
% end
% rpy=fcnCC2SC(rpy/numel(v1)); rpy(:,1)=0;
% x = fcnSC2CC([ones(nframes,1) rpy(:,2:3)*d2r]);
% frames = 1:nframes; ha=fig; hq = quiver3(ha,cam.true.ned(frames,1),cam.true.ned(frames,2),cam.true.ned(frames,3),x(:,1),x(:,2),x(:,3),10,'g');
% box(ha,'on'); axis(ha,'equal','vis3d'); fcn3label(ha); set(ha,'zdir','Reverse');
% hold on; fcnplot3(a.ipned(1:ntp,:),'r.');
% fcnplotrpy(hrpy(3),cam,'1. aposteriori rpy using only SFM')
% cam.aposteriori.rpy = rpy;
% %RUN SFM INITIAL GUESS ----------------------------------------------------

%RUN MSV INITIAL GUESS ----------------------------------------------------
%hrpy=fig(2,2); fcnplotrpy(hrpy(1),cam,'1. apriori rpy')
startclock2 = clock;
[cam, ils.initRMSE] = fcnMSV(cam,a); %fcnplotrpy(hrpy(2),cam,'2. apriori rpy after MSV')
ils.msv.tpnedhat = cam.tpnedhat;
ils.t.init = etime(clock,startclock2);
%RUN MSV INITIAL GUESS ----------------------------------------------------

msvonly = false;
if msvonly
    ils.xhat            = zeros(nframes*12+ntp*3,12);
    ils.xhat12          = zeros(nframes,12);
    ils.xhatsigma12     = zeros(nframes,12);
    ils.xhattp          = zeros(ntp,3);
    ils.xhattpsigma     = zeros(ntp,3);
    ils.initguess       = zeros(nframes,12);
    ils.successflag     = true;
    ils.t.ils           = 0;
    return
end

if ~isfield(a,'iis') %include only specific features in error statistics
   a.iis = true(ntp,1); 
end

fprintf('ILS:  %.0f frames at %.2f fps, %.0f unique tie points (%.1f frames mean lifespan), %.0f measurements (z), %.0f parameters (x), %.0f equations (f)\n',nframes,cam.fps,ntp,mtpl,nz,nx,nf)
fprintf('A = %.0f x %.0f\nB = %.0f x %.0f\nR = %.0f x %.0f\nf = %.0f x %.0f\n\n',nf,nz,nf,nx,nz,nz,nf,1)

%GET UPDATED PHI Q AND P --------------------------------------------------
ekf = fcnT2PPhiQ(cam,ekf,mean(cam.true.dt(2:end))); %fake values passed to ILS, separate from true Phi Q and P
residualrms=zeros(1E3,1);
cam.apriori.rpy = cam.aposteriori.rpy;  xhat12=zeros(nframes,12);
ils.initguess = xhat12;  ils.initguess(:,[7 9 11])=cam.apriori.rpy;
xhat = [reshape(cam.tpnedhat',[ntp*3 1]); reshape(xhat12',[nframes*12 1])];

fprintf('Initial apriori values:\n')
if cam.syntheticVideoFlag;  fprintf('Tie Point Position RMSE (m) = %.2f\nSensor Position RMSE (m) = %.2f\nSensor Angle RMSE (deg RPY) = %.3f\n',1E3*fcnrms(a.ipned(a.iis,:) - cam.tpnedhat(a.iis,:)),1E3*fcnrms(cam.aposteriori.ned-cam.true.ned),fcnrms(fcndrpyd(cam.true.rpy,cam.aposteriori.rpy))); end

I12 = eye(12,12);
ils.successflag = true;
B0 = spalloc(nf,nx,round(nf*nx*.03));
BWe0 = spalloc(nx,nf,round(nf*nx*.03));
R2 = zeros(12*nframes, 12*nframes);
B2 = zeros(12*nframes, nx);
for kk=1:20
    f = zeros(nf,1);
    B = B0;
    BWe = BWe0;
    %R = zeros(nz,nz);
    
    pixelresiduals = zeros([2 ntp nframes]);
    fi = 1;
    zi = 1;
    for i = 1:nframes
        vt = find(activetp(:,i)); %valid tp number indices
        vti3 = [vt*3-2, vt*3-1, vt*3];  vti=vti3';  vti=vti(:); %valid tp ned indices
        
        nvt = numel(vt);  nvf = nvt*2+12;  nvz = nvt*2+12;  nvx = nvt*3+12;
        vf = fi:fi+nvf-1; %valid f indices
        vz = zi:zi+nvz-1; %valid z indices
        vx = [vti' ntp*3+12*(i-1)+(1:12)]; %valid x indices
        
        z1 = [a.upx(vt,i), a.upy(vt,i)]; z1t=z1';
        tpnedhat = xhat(vti3);
        
        %BUILD f --------------------------------------------------------------
        zhat1 = ned2pixel(cam,i,tpnedhat,'aposteriori');  zhat1t=zhat1';  %zhat xy
        f1 = z1t-zhat1t;
        pixelresiduals(:,vt,i) = f1;
        
        %get f2 (frame transitions)
        indcurrent = ntp*3+12*(i-1)+1 : ntp*3+12*(i-1)+12;
        if i>1
            indprevious = indcurrent - 12;
            ekf = fcnT2PPhiQ(cam,ekf,cam.true.dt(i)); %fake values passed to ILS, separate from true Phi Q and P
            f2 = 0 - (xhat(indcurrent) - ekf.Phi12 * xhat(indprevious)); %Phi12
        else
            f2 = -xhat(indcurrent);
        end
        
        f(vf) = [f1(:); f2];
        fi = fi+nvf;
        
        %BUILD A, B and R -----------------------------------------------------
        zvarvec = ones(nvz,1)*ekf.zsigma.^2;
        R1 = diag(zvarvec); %A1 small (no zeros)
        B1 = zeros(nvf, nvx); 
        
        ox = cam.aposteriori.ned(i,1);
        oy = cam.aposteriori.ned(i,2);
        oz = cam.aposteriori.ned(i,3);
        
        rpy = cam.aposteriori.rpy(i,:)*d2r;
        roll=rpy(1);   cr=cos(roll);   sr=sin(roll);
        pitch=rpy(2);  cp=cos(pitch);  sp=sin(pitch);
        yaw=rpy(3);    cy=cos(yaw);    sy=sin(yaw);
        
        for j=1:nvt
            k=(1:2)+(j-1)*2; %fi
            dx=ox-tpnedhat(j,1);  dy=oy-tpnedhat(j,2);  dz=oz-tpnedhat(j,3);  k1 = (cp*cy*dx - sp*dz + cp*sy*dy)^2;

            dpxdtpxyz =  [   (cfl*(cr*sy - cy*sp*sr))/(cp*cy*dx - sp*dz + cp*sy*dy) + (cfl*cp*cy*((cr*cy + sp*sr*sy)*dy - (cr*sy - cy*sp*sr)*dx + cp*sr*dz))/k1, (cfl*cp*sy*((cr*cy + sp*sr*sy)*dy - (cr*sy - cy*sp*sr)*dx + cp*sr*dz))/k1 - (cfl*(cr*cy + sp*sr*sy))/(cp*cy*dx - sp*dz + cp*sy*dy), - (cfl*cp*sr)/(cp*cy*dx - sp*dz + cp*sy*dy) - (cfl*sp*((cr*cy + sp*sr*sy)*dy - (cr*sy - cy*sp*sr)*dx + cp*sr*dz))/k1
                (cfl*cp*cy*((sr*sy + cr*cy*sp)*dx - (cy*sr - cr*sp*sy)*dy + cp*cr*dz))/k1 - (cfl*(sr*sy + cr*cy*sp))/(cp*cy*dx - sp*dz + cp*sy*dy), (cfl*(cy*sr - cr*sp*sy))/(cp*cy*dx - sp*dz + cp*sy*dy) + (cfl*cp*sy*((sr*sy + cr*cy*sp)*dx - (cy*sr - cr*sp*sy)*dy + cp*cr*dz))/k1, - (cfl*cp*cr)/(cp*cy*dx - sp*dz + cp*sy*dy) - (cfl*sp*((sr*sy + cr*cy*sp)*dx - (cy*sr - cr*sp*sy)*dy + cp*cr*dz))/k1];
            
            dpxdrpy = (pi/180)*[  (cfl*((sr*sy + cr*cy*sp)*dx - (cy*sr - cr*sp*sy)*dy + cp*cr*dz))/(cp*cy*dx - sp*dz + cp*sy*dy), (cfl*(cp*cy*sr*dx - sp*sr*dz + cp*sr*sy*dy))/(cp*cy*dx - sp*dz + cp*sy*dy) + (cfl*((cr*cy + sp*sr*sy)*dy - (cr*sy - cy*sp*sr)*dx + cp*sr*dz)*(cp*dz + cy*sp*dx + sp*sy*dy))/k1, - (cfl*((cr*cy + sp*sr*sy)*dx + (cr*sy - cy*sp*sr)*dy))/(cp*cy*dx - sp*dz + cp*sy*dy) - (cfl*(cp*cy*dy - cp*sy*dx)*((cr*cy + sp*sr*sy)*dy - (cr*sy - cy*sp*sr)*dx + cp*sr*dz))/k1
                -(cfl*((cr*cy + sp*sr*sy)*dy - (cr*sy - cy*sp*sr)*dx + cp*sr*dz))/(cp*cy*dx - sp*dz + cp*sy*dy), (cfl*(cp*cr*cy*dx - cr*sp*dz + cp*cr*sy*dy))/(cp*cy*dx - sp*dz + cp*sy*dy) + (cfl*((sr*sy + cr*cy*sp)*dx - (cy*sr - cr*sp*sy)*dy + cp*cr*dz)*(cp*dz + cy*sp*dx + sp*sy*dy))/k1,   (cfl*((cy*sr - cr*sp*sy)*dx + (sr*sy + cr*cy*sp)*dy))/(cp*cy*dx - sp*dz + cp*sy*dy) - (cfl*(cp*cy*dy - cp*sy*dx)*((sr*sy + cr*cy*sp)*dx - (cy*sr - cr*sp*sy)*dy + cp*cr*dz))/k1];
            
%             tp = tpnedhat(j,:);  pdx = .000001;  flag = 'aposteriori';
%             dpxdtpxyz = [     (ned2pixel(cam,i,tp+[pdx 0 0],flag)-ned2pixel(cam,i,tp,flag))./pdx
%                               (ned2pixel(cam,i,tp+[0 pdx 0],flag)-ned2pixel(cam,i,tp,flag))./pdx
%                               (ned2pixel(cam,i,tp+[0 0 pdx],flag)-ned2pixel(cam,i,tp,flag))./pdx      ]';  %true numerical derivative

            B1(k,j*3-[2 1 0]) = dpxdtpxyz;
            B1(k,nvt*3+[1 3 5]) = -dpxdtpxyz;
            B1(k,nvt*3+[7 9 11]) = dpxdrpy;
        end
        B1(nvt*2+(1:12),nvt*3+(1:12)) = I12; %current frame B1
        
        B(vf,vx) = B1;     
        vz1 = 1:nvt*2;
        vz2 = (nvt*2+1):nvz;
        if i==1
            R1(vz2,vz2) = ekf.P0_12;
            R2(1:12,1:12) = ekf.P0_12;
        else
            R1(vz2,vz2) = ekf.Q12;
            R2((i-1)*12+(1:12),(i-1)*12+(1:12)) = ekf.Q12;
            B(vf(vz2),vx(end-11:end)-12) = -ekf.Phi12;
        end
        %R(vz,vz) = R1;
        B2((i-1)*12+(1:12),:) = B(vf(vz2),:);       

        %BWe(:,vf) = B(vf,:)' / R1;
        %BWe(:,vf(vz1)) = bsxfun(@times, B(vf(vz1),:), 1./zvarvec(vz1))';
        BWe(:,vf(vz1)) = B(vf(vz1),:)' * (1./ekf.zsigma.^2);
        BWe(:,vf(vz2)) = B(vf(vz2),:)' / R1(vz2,vz2);
        zi = zi + nvz;
    end
    %j = repmat((1:nvz)',[1 nframes]) + repmat(0:nvz:nf-1,[nvz 1]);  j=j(1:nvt*2,:);  j=j(:);  BWe(:,j) = B(j,:)' * (1./ekf.zsigma.^2);

    
    residualrms(kk) = fcnrms(pixelresiduals(vp2));
    fprintf('Pixel Residuals RMS = %.4f\n\n',residualrms(kk))

    if kk==1 %Pixel Residuals
        [~,hpr] = fig;
        y = pixelresiduals;  y = sqrt(y(1,:,:).^2 + y(2,:,:).^2);  y = sum(y,3)./sum(activetp'); %#ok<UDIM>
        bar(1:ntp, y,1, 'r','edgecolor',[.7 .7 .7]); alpha(0.5); hold on;
        plot([1 ntp], mean(y)*[1 1],'r','linewidth',2); axis tight
        plot([1 ntp], mean(y)*[1 1]+1*std(y)*[1 1],'--r','linewidth',1);
        plot([1 ntp], mean(y)*[1 1]-1*std(y)*[1 1],'--r','linewidth',1);
        xlabel('Tie Point ID'); ylabel('Mean Residual (pixels)'); title(sprintf('mean tie point residuals for %s across %.0f frames',cam.filename,nframes))
        
        [~,hpv] = fig;
        z1 = [a.upx(vt,i), a.upy(vt,i)]; f2 = (z1 - zhat1)*10;
        plot(z1(:,1),z1(:,2),'r.'); hold on;
        hpvl(1)=quiver(z1(:,1),z1(:,2),f2(:,1),f2(:,2),0,'r.');
    end
    
    %BWe = B'/R;
    N = BWe*B; %sigma = inv(N);  %BWe = B'*W;
    xhat = xhat + N\(BWe*f  - B2'*(R2\B2)*xhat*0);
    
    tpnedhat = reshape(xhat(1:ntp*3,1),[3 ntp])';
    cam.aposteriori.ned = cam.apriori.ned + xhat(  ntp*3 + [12*(0:nframes-1)'+1 12*(0:nframes-1)'+3 12*(0:nframes-1)'+5 ]  );
    cam.aposteriori.rpy = cam.apriori.rpy + xhat(  ntp*3 + [12*(0:nframes-1)'+7 12*(0:nframes-1)'+9 12*(0:nframes-1)'+11]  );
    
    fprintf('Iteration %.0f results:\n',kk)
    if cam.syntheticVideoFlag;  fprintf('Tie Point Position RMSE (m) = %.2f\nSensor Position RMSE (m) = %.2f\nSensor Angle RMSE (deg RPY) = %.3f\n',1E3*fcnrms(a.ipned(a.iis,:) - tpnedhat(a.iis,:)),1E3*fcnrms(cam.aposteriori.ned-cam.true.ned),fcnrms(fcndrpyd(cam.true.rpy,cam.aposteriori.rpy)));  end
    
    if kk>1
        df = abs(residualrms(kk-1)-residualrms(kk))/residualrms(kk-1); %fractional change
        if df<.001 || residualrms(kk)>1000;  break;  end  %converged or diverged!
        if kk>2
            if residualrms(kk)>residualrms(kk-1) && residualrms(kk-1)>residualrms(kk-2);  ils.successflag=false;  break; end
        end
    end
end
y = pixelresiduals;  y = sqrt(y(1,:,:).^2 + y(2,:,:).^2);  y = sum(y,3)./sum(activetp');
ils.pixelresiduals = y;

%Pixel Residuals
figure(hpr)
bar(1:ntp, y,1, 'b','edgecolor',[.7 .7 .7]); alpha(0.5); hold on;
plot([1 ntp], mean(y)*[1 1],'b','linewidth',2); axis tight
plot([1 ntp], mean(y)*[1 1]+1*std(y)*[1 1],'--b','linewidth',1);
plot([1 ntp], mean(y)*[1 1]-1*std(y)*[1 1],'--b','linewidth',1);
xlabel('Tie Point ID'); ylabel('Mean Residual (pixels)');
legend('{\ita priori} residuals','{\ita priori} mean','{\ita priori} +1\sigma','{\ita priori} -1\sigma',...
    '{\ita posteriori} residuals','{\ita posteriori} mean','{\ita posteriori} +1\sigma','{\ita posteriori} -1\sigma','location','best');

figure(hpv)
z1 = [a.upx(vt,i), a.upy(vt,i)]; f2 = (z1 - zhat1)*10;
hpvl(2)=quiver(z1(:,1),z1(:,2),f2(:,1),f2(:,2),0,'b.');
axis([1 cam.width 1 cam.height]); legend(hpvl,'apriori','aposteriori'); xlabel('x (pixels)'); ylabel('y (pixels)'); set(gca,'ydir','reverse')


ils.Sigma = inv(full(N));  %Covariance Matrix = N^-1
diagSigma = diag(ils.Sigma);
ils.xhat = xhat;
ils.xhat12 = reshape(xhat(ntp*3+1:end),[12 nframes])';
ils.xhat12Sigma = ils.Sigma(ntp*3+1:end,ntp*3+1:end);
ils.xhatsigma12 = reshape(sqrt(diagSigma(ntp*3+1:end)), [12 nframes])';
j = reshape(1:ntp*3,[3 ntp])';
ils.xhattp = reshape(xhat(1:ntp*3),[3 ntp])';
ils.xhattpsigma = reshape(sqrt(diagSigma(1:ntp*3)), [3 ntp])';
ils.residuals = y;
ils.t.ils = etime(clock,startclock1) - ils.t.init;

m1 = 0:12:(12*(nframes-1));
ils.index.ned = ones(nframes,1)*[1 3 5] + m1'*[1 1 1] + ntp*3;
ils.index.rpy = ones(nframes,1)*[7 9 11] + m1'*[1 1 1] + ntp*3;

if ~ils.successflag || residualrms(kk)>30 || kk==10 || ~isreal(ils.xhat) || ~isreal(diagSigma)
    ils.successflag = false;    
    fprintf('\nWARNING: ILS FAILURE after %.fs.\n',ils.t.ils)
else
    fprintf('\nILS completed successfully in %.fs. Covariance Matrix Uncertainties:\n',ils.t.ils)
    fprintf('Tie Point Position 1sigma (m NED) = [%.3g  %.3g  %.3g]\nSensor Position 1sigma (m NED) = [%.3g  %.3g  %.3g]\nSensor Angle 1sigma (deg RPY) = [%.3g  %.3g  %.3g]\n', ...
        1E3*sqrt(mean(ils.xhattpsigma.^2)),1E3*sqrt(mean(ils.xhatsigma12(:,[1 3 5]).^2)),mean(ils.xhatsigma12(:,[7 9 11])));
end

%fcncompareILS(cam,ekf,ils,a)

% ha=fig;
% x=1:cam.frames;
% y=ekf.x([7 9 11],x)';
% plot(ha,x,y(:,1),'r',x,y(:,2),'g',x,y(:,3),'b','linewidth',2,'linesmoothing','on'); hold on;
% y=ils.xhat12(x,[7 9 11]);
% plot(ha,x,y(:,1),'r',x,y(:,2),'g',x,y(:,3),'b','linewidth',1,'linesmoothing','on'); hold on;
% str = sprintf('%.2fdeg RMSE',fcnrms(unwrap(cam.aposteriori.rpy*d2r-cam.true.rpy*d2r)*r2d));
% title(ha,'rpy'); legend(ha,'roll','pitch','yaw'); fcn3label(ha,'frames','angle (deg RPY)')

%fcnplotrpy(hrpy(4),cam,'1. aposteriori rpy after ILS'); fcnfontsize(8)

% h = fig(2,3);
% axes(h(1)); pcolor(double(B~=0)); shading flat; set(gca,'ydir','reverse'); axis equal tight; xlabel('Z'), ylabel('F'); title('B matrix'); colormap(1-gray); colorbar
% axes(h(2)); pcolor(double(Bt~=0)); shading flat; set(gca,'ydir','reverse'); axis equal tight; xlabel('X'), ylabel('F'); title('Bt matrix'); colorbar
% axes(h(3)); pcolor(double(R~=0)); shading flat; set(gca,'ydir','reverse'); axis equal tight; xlabel('F'), ylabel('F'); title('R matrix'); colorbar
% axes(h(4)); pcolor(double(W~=0)); shading flat; set(gca,'ydir','reverse'); axis equal tight; xlabel('F'), ylabel('F'); title('W matrix'); colorbar
% axes(h(5)); pcolor(double(N~=0)); shading flat; set(gca,'ydir','reverse'); axis equal tight; xlabel('F'), ylabel('F'); title('N matrix'); colorbar
% axes(h(6)); pcolor(double(inv(N)~=0)); shading flat; set(gca,'ydir','reverse'); axis equal tight; xlabel('F'), ylabel('F'); title('N^{-1} matrix'); colorbar

% h = fig(2,3);
% axes(h(1)); pcolor(double(B)); shading flat; set(gca,'ydir','reverse'); axis equal tight; xlabel('Z'), ylabel('F'); title('A matrix'); colormap(jet); colorbar
% axes(h(2)); pcolor(double(B)); shading flat; set(gca,'ydir','reverse'); axis equal tight; xlabel('X'), ylabel('F'); title('B matrix'); colorbar
% axes(h(3)); pcolor(double(R)); shading flat; set(gca,'ydir','reverse'); axis equal tight; xlabel('F'), ylabel('F'); title('R matrix'); colorbar
% axes(h(4)); pcolor(double(W)); shading flat; set(gca,'ydir','reverse'); axis equal tight; xlabel('F'), ylabel('F'); title('Re matrix'); colorbar
% axes(h(5)); pcolor(double(N)); shading flat; set(gca,'ydir','reverse'); axis equal tight; xlabel('F'), ylabel('F'); title('N matrix'); colorbar
% axes(h(6)); pcolor(double(inv(N))); shading flat; set(gca,'ydir','reverse'); axis equal tight; xlabel('F'), ylabel('F'); title('N^{-1} matrix'); colorbar

end

function [] = fcnplotrpy(ha,cam,titlestr)
x=1:cam.frames;
z=cam.aposteriori.rpy;
plot(ha,x,z(:,1),'r',x,z(:,2),'g',x,z(:,3),'b','linesmoothing','on'); hold(ha,'on')
y=cam.true.rpy;
plot(ha,x,y(:,1),'r',x,y(:,2),'g',x,y(:,3),'b','linewidth',2,'linesmoothing','on');
str = sprintf('%s, %.2fdeg RMSE',titlestr,fcnrms(unwrap(cam.aposteriori.rpy*d2r-cam.true.rpy*d2r)*r2d));
title(ha,str); legend(ha,'roll','pitch','yaw'); fcn3label(ha,'frames','angle (deg)')
end
