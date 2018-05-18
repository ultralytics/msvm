function [ils, cam] = fcnLM(cam,ekf,a)
startclock1 = clock;
nframes = cam.frames;
cam.apriori.ned = cam.apriori.ned(1:nframes,:);   cam.aposteriori.ned = cam.apriori.ned;  cam.true.ned = cam.true.ned(1:nframes,:);   
cam.apriori.rpy = cam.apriori.rpy(1:nframes,:);   cam.aposteriori.rpy = cam.apriori.rpy;  cam.true.rpy = cam.true.rpy(1:nframes,:);

%RUN MSV INITIAL GUESS ----------------------------------------------------
startclock2 = clock;
%[cam, ils.initRMSE] = fcnMSV(cam,a);
ils.t.init = etime(clock,startclock2);
cam.apriori.rpy = cam.aposteriori.rpy;  tpnedhat=cam.tpnedhat;

%GET SHORT TRACKS ---------------------------------------------------------
%i = find(sum(a.state==1,2) < (nframes*.9));   a = fcncropTP(cam,a,i);  cam.tpnedhat = cam.tpnedhat(i,:);

%DEFINE MATRICES ----------------------------------------------------------
ntp = size(a.upx,1); %number of tie points
activetp = a.state==1;  clear vp2;  vp2(1,:,:) = activetp;  vp2(2,:,:) = activetp;  %valid pixel indices for pixelresiduals
ntpe = sum3(activetp);%number of tie point equations
mtpl = mean(sum(activetp,2)); %mean tie point life
ncp = 6; %number of camera parameters per frame
ncpK = 3; %number of camera parameters (xppo, yppo, focalLength)
nz = ntpe*2 + nframes*3   + ncpK; %number of measurments total
nx = ntp*3  + nframes*ncp + ncpK; %number of parameters
nf = ntpe*2 + nframes*3   + ncpK; %number of equations
xhat = [tpnedhat(:); cam.apriori.rpy(:); cam.apriori.ned(:); cam.xppo; cam.yppo; cam.focalLength];
index.tp = reshape(1:ntp*3,[ntp 3]);
index.rpy = reshape(ntp*3+(1:nframes*3),[nframes 3]);
index.ned = reshape(ntp*3+nframes*3+(1:nframes*3),[nframes 3]);
index.K = (ntp*3+ncp*nframes) + [1 2 3];

if numel(ekf.zsigma)==1;
    a.R = ones(size(a.state))*ekf.zsigma.^2;
%     age = cumsum(a.state'==1)';  %figure; pcolor(age); shading flat;
%     mage = mean(age(activetp)); %mean age
%     a.R = (age/mage * ekf.zsigma).^2 + 1;
else
    a.R = ekf.zsigma.^2 * ones(size(a.state(1,:)));
end

fprintf('ILS:  %.0f frames at %.2f fps, %.0f unique tie points (%.1f frames mean lifespan), %.0f measurements (z), %.0f parameters (x), %.0f equations (f)\n',nframes,cam.fps,ntp,mtpl,nz,nx,nf)
fprintf('A = %.0f x %.0f\nB = %.0f x %.0f\nR = %.0f x %.0f\nf = %.0f x %.0f\n',nf,nz,nf,nx,nz,nz,nf,1)
fprintf('              Pixel Residuals RMS   Tie-Point Position RMSE (m)      Sensor Position RMSE (m)       Sensor Angle RMSE (deg)\n');  
if cam.syntheticVideoFlag; fprintf('Apriori:                         '); fprintf('%30.2f%30.2f%30.3f\n',1E3*fcnrms(a.ipned(a.iis,:)-cam.tpnedhat(a.iis,:)),1E3*fcnrms(cam.aposteriori.ned-cam.true.ned),fcnrms(fcndrpyd(cam.true.rpy,cam.aposteriori.rpy))); end

ils.successflag = true;
maxi = 10;  
RMSE=zeros(maxi,1);
for i=1:10
    [f, pixelresiduals, BWe, B]=geterror(cam,a,xhat,ekf,nf,nx,ntp,activetp,nframes,index,ncp,ncpK,2);

%     if i==1 %Pixel Residuals
%         [hpr,~] = fig(1,2,2); 
%         axes(hpr(1));   y = pixelresiduals;  y = sqrt(y(1,:,:).^2 + y(2,:,:).^2);  y = sum(y,3)./tpl; %#ok<*LAXES,UDIM>
%         bar(1:ntp, y,1, 'r','edgecolor',[.7 .7 .7]); alpha(0.5); hold on;
%         plot([1 ntp], mean(y)*[1 1],'r','linewidth',2); axis tight
%         plot([1 ntp], mean(y)*[1 1]+1*std(y)*[1 1],'--r','linewidth',1);
%         plot([1 ntp], mean(y)*[1 1]-1*std(y)*[1 1],'--r','linewidth',1);
%         xlabel('Tie Point ID'); ylabel('Mean Residual (pixels)'); title(sprintf('mean tie point residuals for %s across %.0f frames',cam.filename,nframes))
%         axes(hpr(2));  i=1;  vt=find(activetp(:,i));  tpnedhat=xhat(index.tp(vt,:));  zhat1=ned2pixel(cam,i,tpnedhat,'aposteriori');  z1=[a.upx(vt,i) a.upy(vt,i)];  f2=(z1-zhat1)*10;
%         plot(z1(:,1),z1(:,2),'r.'); hold on;
%         hpvl(1)=quiver(z1(:,1),z1(:,2),f2(:,1),f2(:,2),0,'r.');
%     end
    
    %UPDATE ---------------------------------------------------------------
    H = full(BWe*B); 
    BWef = BWe*f;
    
%     %LINE SEARCH DELTA
%     g = logspace(-1,1,10);  ng=numel(g);  gy=zeros(1,ng);
%     for j=1:ng
%         %gy(j) = fcnrms(geterror(cam,a,xhat-(H + g(j)*eye(nx))\BWef,ekf,nf,nx,ntp,activetp,nframes,index,ncp,ncpK,1));
%         gy(j) = fcnrms(geterror(cam,a,xhat-g(j)*H\BWef,ekf,nf,nx,ntp,activetp,nframes,index,ncp,ncpK,1));
%     end
%     [~,j] = min(gy);  %figure; plot(g,gy);
%     %delta = -(H + g(j)*eye(nx))\BWef;
%     delta = -g(j)*H\BWef;
% 
%    %LM DELTA
%     if i==1
%         lambda = 1E-3 * mean(diag(H));
%     else
%         lambda = lambda/10;
%     end
%     delta = -(H + lambda*eye(nx))\BWef;
%     ftest=fcnrms(geterror(cam,a,xhat+delta,ekf,nf,nx,ntp,activetp,nframes,index,ncp,ncpK,1));
%     while ftest>fcnrms(f)
%         lambda = lambda*10;
%         delta = -(H + lambda*eye(nx))\BWef;
%         ftest=fcnrms(geterror(cam,a,xhat+delta,ekf,nf,nx,ntp,activetp,nframes,index,ncp,ncpK,1));
%     end

    %NORMAL DELTA
    delta = -H\BWef;

    %UPDATE
    xhat = xhat + delta;
    
    tpnedhat=xhat(index.tp);
    cam.aposteriori.ned = xhat(index.ned);
    cam.aposteriori.rpy = xhat(index.rpy);
    cam.xppo = xhat(index.K(1));  cam.yppo = xhat(index.K(2));  cam.focalLength = xhat(index.K(3));
    
    RMSE(i) = fcnrms(pixelresiduals(vp2));
    fprintf('Iteration %.0f: %20.4f',i,RMSE(i));
    if cam.syntheticVideoFlag;  fprintf('%30.2f%30.2f%30.3f',1E3*fcnrms(a.ipned(a.iis,:)-tpnedhat(a.iis,:)),1E3*fcnrms(cam.aposteriori.ned-cam.true.ned),fcnrms(fcndrpyd(cam.true.rpy,cam.aposteriori.rpy))); end
    fprintf('\n')
    
    if i>1
        %%CHECK FOR CONVERGENCE OR DIVERGENCE
        df = abs(RMSE(i-1)-RMSE(i))/RMSE(i-1); %fractional change
        if df<.0001 || RMSE(i)>1000;  break;  end  %converged or diverged!
        %if i>2 && RMSE(i)>RMSE(i-1) && RMSE(i-1)>RMSE(i-2);  ils.successflag=false;  break; end
        
%         %REMOVE OUTLIERS
%         x = pixelresiduals;  x = sqrt(x(1,:,:).^2 + x(2,:,:).^2);  x = sum(x,3)./sum(activetp'); %#ok<UDIM>     
%         j = a.R(:,1)==ekf.zsigma^2;
%         mu = mean(x(j));  sigma = std(x(j));  j = x>(mu+2*sigma);  a.R(j,:) = 100^2;
%         %fcnValidateResults(cam,a);
    end
    
    
end
fprintf('Estimated Camera Parameters:%10.1f%10.1f%10.1f',cam.xppo,cam.yppo,cam.focalLength)
fprintf('\n')

% %Pixel Residuals
% y = pixelresiduals;  y = sqrt(y(1,:,:).^2 + y(2,:,:).^2);  y = sum(y,3)./sum(activetp'); %#ok<UDIM>
% axes(hpr(1));   bar(1:ntp, y,1, 'b','edgecolor',[.7 .7 .7]); alpha(0.5); hold on;
% plot([1 ntp], mean(y)*[1 1],'b','linewidth',2); axis tight
% plot([1 ntp], mean(y)*[1 1]+1*std(y)*[1 1],'--b','linewidth',1);
% plot([1 ntp], mean(y)*[1 1]-1*std(y)*[1 1],'--b','linewidth',1);
% xlabel('Tie Point ID'); ylabel('Mean Residual (pixels)');
% legend('{\ita priori} residuals','{\ita priori} mean','{\ita priori} +1\sigma','{\ita priori} -1\sigma',...
%     '{\ita posteriori} residuals','{\ita posteriori} mean','{\ita posteriori} +1\sigma','{\ita posteriori} -1\sigma','location','best');
% axes(hpr(2));  i=1;  vt=find(activetp(:,i));  tpnedhat=xhat(index.tp(vt,:));  zhat1=ned2pixel(cam,i,tpnedhat,'aposteriori');  z1=[a.upx(vt,i) a.upy(vt,i)];  f2=(z1-zhat1)*10;
% hpvl(2)=quiver(z1(:,1),z1(:,2),f2(:,1),f2(:,2),0,'b.');
% axis([1 cam.width 1 cam.height]); legend(hpvl,'apriori','aposteriori'); xlabel('x (pixels)'); ylabel('y (pixels)'); set(gca,'ydir','reverse')

C = inv(full(H));
diagSigma = sqrt(diag(C)); %Covariance Matrix = N^-1
ils.xhat = xhat;
ils.xhatKsigma = diagSigma(index.K);
ils.xhatnedsigma = diagSigma(index.ned);
ils.xhatrpysigma = diagSigma(index.rpy);
ils.xhattp = xhat(index.tp);
ils.xhattpsigma = diagSigma(index.tp);
ils.xhattpC = C(index.tp,index.tp);
ils.xhatC = C;
ils.index = index;

% i = [ils.index.rpy(end,:) ils.index.ned(end,:) ils.index.K];
% C = ils.xhatC(i,i);
% D = cov2corr(C);
% for i=1:numel(C(:,1)); corrC(i,i) = sigmas(i); end;

y = pixelresiduals; y = sqrt(y(1,:,:).^2 + y(2,:,:).^2);  y = sum(y,3)./sum(activetp'); %#ok<UDIM>
ils.pixelresiduals = y;
ils.t.ils = etime(clock,startclock1) - ils.t.init;

if ~ils.successflag || RMSE(i)>30 || i==maxi || ~isreal(ils.xhat) %|| ~isreal(diagSigma)
    ils.successflag = false;    
    fprintf('\nWARNING: ILS FAILURE after %.fs.\n',ils.t.ils)
else
    fprintf('\nILS completed successfully in %.fs\n',ils.t.ils)
end
fprintf('Covariance Matrix Uncertainties:\n')
fprintf('Tie Point Position 1sigma (m NED) = [%.3g  %.3g  %.3g]\nSensor Position 1sigma (m NED) = [%.3g  %.3g  %.3g]\nSensor Angle 1sigma (deg) = [%.3g  %.3g  %.3g]\n', ...
    1E3*sqrt(mean(ils.xhattpsigma.^2)),1E3*sqrt(mean(ils.xhatnedsigma.^2)),sqrt(mean(ils.xhatrpysigma.^2)));
end

function [f, pixelresiduals, BWe, B]=geterror(cam,a,xhat,ekf,nf,nx,ntp,activetp,nframes,index,ncp,ncpK,mode)
%mode = 1 (f only), 2 (f, B, and BWe)
cam.aposteriori.ned = xhat(index.ned);
cam.aposteriori.rpy = xhat(index.rpy);
cam.xppo = xhat(index.K(1));  cam.yppo = xhat(index.K(2));  cam.focalLength = xhat(index.K(3));    cfl=cam.focalLength;
RK = diag(cam.sigmas.K.^2); %[xppo yppo focalLength], pixels
d2r = pi/180;  r2d=180/pi;

vfC = zeros(nframes,3); %vector of Camera position functions
f = zeros(nf,1);  %R = zeros(nz,nz);
%R = zeros(nf,nf);
B = spalloc(nf,nx,round(nf*nx*.03));
BWe = spalloc(nx,nf,round(nf*nx*.03));
pixelresiduals = zeros([2 ntp nframes]);

fi = 1;
zi = 1;
for i = 1:nframes
    vt = find(activetp(:,i)); %valid tp number indices
    nvt = numel(vt);  nvf = nvt*2+3;  nvz = nvt*2+3;  nvx = nvt*3+ncp+ncpK;  vfC(i,:) = (fi+nvt*2) + [0 1 2];
    
    vf = fi:fi+nvf-1; %valid f indices
    vx = [reshape(index.tp(vt,:)',[nvt*3 1])' index.rpy(i,:)  index.ned(i,:)  index.K]; %valid x indices
    
    z1 = [a.upx(vt,i), a.upy(vt,i)]; z1t=z1';
    tpnedhat = xhat(index.tp(vt,:));
    
    %BUILD f --------------------------------------------------------------
    zhat1 = ned2pixel(cam,i,tpnedhat,'aposteriori');  zhat1t=zhat1';  %zhat xy
    f1 = zhat1t - z1t;
    pixelresiduals(:,vt,i) = f1;
    
    f2 = xhat(index.ned(i,:)) - cam.apriori.ned(i,:)'; %position correction
    f(vf) = [f1(:); f2];
    fi = fi+nvf;
    
    %BUILD A, B and R -----------------------------------------------------
    if mode==2
        %rvec = (ones(nvz,1)*ekf.zsigma.^2);  R1 = diag(rvec);
        rvec = repmat(a.R(vt,i)',[2 1]);  R1 = diag([rvec(:)' 0 0 0]);
        B1 = zeros(nvf, nvx);
        
        ox = cam.aposteriori.ned(i,1);
        oy = cam.aposteriori.ned(i,2);
        oz = cam.aposteriori.ned(i,3);
        rpy = cam.aposteriori.rpy(i,:)*d2r;
        roll=rpy(1);   cr=cos(roll);   sr=sin(roll);
        pitch=rpy(2);  cp=cos(pitch);  sp=sin(pitch);
        yaw=rpy(3);    cy=cos(yaw);    sy=sin(yaw);
        for j=1:nvt
            k=(1:2)+(j-1)*2;
            dx=ox-tpnedhat(j,1);  dy=oy-tpnedhat(j,2);  dz=oz-tpnedhat(j,3);  
            k1 = cp*cy*dx - sp*dz + cp*sy*dy;
            k1s = k1^2;
            k2 = cr*sy - cy*sp*sr;
            k3 = cr*cy + sp*sr*sy;
            k4 = sr*sy + cr*cy*sp;
            k5 = cy*sr - cr*sp*sy;
            k6 = cp*sr*dz;
            k7 = cp*cr*dz;

            dpxdtpxyz = [   (cfl*k2)/k1 + (cfl*cp*cy*(k3*dy - k2*dx + k6))/k1s, (cfl*cp*sy*(k3*dy - k2*dx + k6))/k1s - (cfl*k3)/k1, - (cfl*cp*sr)/k1 - (cfl*sp*(k3*dy - k2*dx + k6))/k1s
                (cfl*cp*cy*(k4*dx - k5*dy + k7))/k1s - (cfl*k4)/k1, (cfl*k5)/k1 + (cfl*cp*sy*(k4*dx - k5*dy + k7))/k1s, - (cfl*cp*cr)/k1 - (cfl*sp*(k4*dx - k5*dy + k7))/k1s];
            dpxdrpy = d2r*[    (cfl*(k4*dx - k5*dy + k7))/k1, (cfl*(cp*cy*sr*dx - sp*sr*dz + cp*sr*sy*dy))/k1 + (cfl*(k3*dy - k2*dx + k6)*(cp*dz + cy*sp*dx + sp*sy*dy))/k1s, - (cfl*(k3*dx + k2*dy))/k1 - (cfl*(cp*cy*dy - cp*sy*dx)*(k3*dy - k2*dx + k6))/k1s
                -(cfl*(k3*dy - k2*dx + k6))/k1, (cfl*(cp*cr*cy*dx - cr*sp*dz + cp*cr*sy*dy))/k1 + (cfl*(k4*dx - k5*dy + k7)*(cp*dz + cy*sp*dx + sp*sy*dy))/k1s,   (cfl*(k5*dx + k4*dy))/k1 - (cfl*(cp*cy*dy - cp*sy*dx)*(k4*dx - k5*dy + k7))/k1s];
            dpxdK = [ 1,  0,  (k3*dy - k2*dx + k6)/k1
                0,  1,  (k4*dx - k5*dy + k7)/k1];
            
            B1(k,j*3-[2 1 0]) = dpxdtpxyz;  % tiepoint pos
            B1(k,nvt*3+[1 2 3]) = dpxdrpy;  % camera rpy
            B1(k,nvt*3+[4 5 6]) = -dpxdtpxyz;  %camera pos
            B1(k,nvt*3+ncp+[1 2 3])= dpxdK;  %K = [xppo, yppo, focalLength];
        end
        B1(nvt*2+(1:3),nvt*3+[4 5 6]) = eye(3); %current frame B1
        R1(end-2:end,end-2:end) = diag(([cam.ilsxy1s cam.ilsxy1s cam.ilsz1s]/1E3).^2);
        B(vf,vx) = B1;
        %R(vf,vf) = R1;
        BWe(vx,vf) = sparse(B1') / sparse(R1);
    end
    zi = zi + nvz;
end

%CAMERA PARAMETERS --------------------------------------------------------
vf3 = fi+[0 1 2];  fi=fi+3;  zi=zi+3;
f(vf3) = xhat(index.K) - [cam.apriori.K.xppo; cam.apriori.K.yppo; cam.apriori.K.focalLength];  if mode==1; return; end
B(vf3,index.K) = eye(3);
%R(vf3,vf3) = RK;
BWe(:,vf3) = B(vf3,:)' / RK;

%GPS CORRELATIONS ---------------------------------------------------------
s = [cam.ilsxy1s cam.ilsxy1s cam.ilsz1s]/1E3;
for i=1:3 
   j = vfC(:,i); %function indices
   R2 = fcnGPScov(cam, s(i));
   BWe(:,j) = B(j,:)' / R2;
end

end

function cov = fcnGPScov(cam, sigma)
t = cam.true.t;  n = numel(t);  ov = ones(1,n);
T = 1800; %s
dtm = abs(ov'*t - t'*ov); %dt matrix
corr = exp(-dtm/T);
cov = corr*sigma^2;  %cov = corr2cov(sigma*ov, corr)
end
