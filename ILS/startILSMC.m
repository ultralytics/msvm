clc; clear; close all; opengl('OpenGLLineSmoothingBug',1); startclock1 = clock;  load test20tp.mat
cam.optim.options1 = optimset('GradObj','off','display','off','largescale','off','MaxFunEvals',2000);
cam.optim.options2 = optimset('GradObj','off','display','off','MaxFunEvals',2000,'HessUpdate','bfgs','Algorithm','active-set');
cam.optim.options3 = optimset('Algorithm','levenberg-marquardt','display','off');
cam.optim.options4 = optimset('Algorithm','levenberg-marquardt','display','off','Jacobian','on','DerivativeCheck','off','gradobj','on');

%tsv = [4:2:30];  xstr = 'Angle 1\sigma (s)';  fname = 'TS - Angle 1sigma 5degroll 1cube FPW TPCENTER 9FRAMES';
%tsv = 1:2:21;  xstr = 'Angle Time Constant (s)';  fname = 'TS - Angle Time Constant 5degroll 1cube FPW';
%tsv = round([10 20 30 40 50]);  xstr = 'Number of MSV Grid Points';  fname = 'TS - MSV Grid Points 5degroll 1cube FPW SQUAREAD work';
%tsv = round([2 4 6 8 10 12 14 16 18 20]);  %xstr = 'ILS-fed {\itapriori} rpy 1\sigma (deg) (\tau=.1s)'; fname = 'TS - .1s rpy time constant';
%tsv = round([50]);  xstr = 'Number of Frames (Jocher-minsepr)';  fname = 'Thirty Frames';
%tsv = round([5]);  xstr = 'Camera Uncertainty (m) (Hank-minsepangle)';  fname = 'Thirty Frames';
%tsv = round(20:20:200);  xstr = 'TP minsepr2 (all 40 frames  used)';  fname = 'minsepr';
tsv = round([1:10]);  xstr = '{\bfTrue} Pixel 1\sigma (assumed=3pix)';  fname = 'Thirty Frames';

nts = numel(tsv);
ts.sigma12mean = zeros(nts,12);   ts.sigma12std = zeros(nts,12);   ts.xhat12mean = ts.sigma12mean;   ts.xhat12std = ts.sigma12std;   ts.ilssuccessfraction = zeros(1,nts);  ts.timeperMC = zeros(nts,2);
ts.sigma3mean = zeros(nts,3);     ts.sigma3std = zeros(nts,3);     ts.xhat3mean = ts.sigma3mean;     ts.xhat3std = ts.sigma3std;     ts.initRMSEmean = zeros(nts,6);        ts.initRMSEstd = ts.initRMSEmean;

nmc = 1000; %number MCs
cam.msvgridpoints = 30;
cam.fps = .5; %must divide equally into 30fps
cam.fovh = 5; %deg
cam.xyz1s = 2; %5m true
cam.rpy1s = 20; %2deg
cam.roll1s = 5; %deg
cam.focuspointwander = 30; %m 1 sigma
cam.fovv=cam.fovh*(cam.height/cam.width);
cam.width=cam.ssge(3);  cam.height=cam.ssge(4);
cam.true.focus.pixel = camsc2pixel(cam, [1 0 0]);
cam.ilsrpy1s = 10;  %5deg, new 15deg
cam.ilsTrpy = 5; %10s, new 2s
cam.sfmframes = 10;
cam.true.LLAstart = [40.4170, -3.7030, 0]; %Madrid Sol
cam0 = fcnaircraftpath(DEM,cam);  a0=a;

if nmc>30;  matlabpool close force local;  matlabpool open 4;  end
for i=1:nts
    nf = 80; %tsv(i); %number frames
    ntp = 20; %number of tie points
    mc_xhat12 = zeros(nf,12,nmc);   mc_xhatsigma12 = mc_xhat12;  vs = false(1,nmc);
    mc_xhattp = zeros(ntp,3,nmc);   mc_xhattpsigma = mc_xhattp;  mc_initRMSE = zeros(nmc,6);  mc_time=zeros(nmc,2);
    startclock2 = clock;
    
    parfor j=1:nmc
        a=a0;
        
        %CAMERA TRAJECTORY ----------------------------------------------------
        [cam, ekf] = fcnaircraftpath(DEM,cam0);
        
        %GET TIE POINTS -------------------------------------------------------
        ekf.zsigma = 3; %1pixels
        tp1 = [rand(ntp,1)*cam.width,  rand(ntp,1)*cam.height];  ov=ones(ntp,1);
        a.ipned = (pixel2ned(DEM,cam,10,tp1,'true') + randn(ntp,3)*.000).*[ov ov ov];
        a.upx=zeros(ntp,cam.frames);  a.upy=a.upx;
        for k=1:cam.frames
            xy = ned2pixel(cam,k,a.ipned,'true');
            a.upx(:,k) = xy(:,1);
            a.upy(:,k) = xy(:,2);
        end
        %a.upx=a.upx+randn(size(a.upx))*ekf.zsigma;  a.upy=a.upy+randn(size(a.upx))*ekf.zsigma; %add noise
        a.upx=a.upx+randn(size(a.upx))*tsv(i);  a.upy=a.upy+randn(size(a.upx))*tsv(i); %add noise
        
        %ILS ------------------------------------------------------------------        
        cam.frames = nf;
        %fcnMSV(cam,a); drawnow
        %fcnMSVhank(cam,a); drawnow
        ils = fcnstartILS(cam,ekf,a);
        
        %SAVE -----------------------------------------------------------------
        ia = [7 9 11];
        xtrue = ekf.x(:,1:nf)';  xtrue(:,ia) = cam.true.rpy(1:nf,:) - ils.initguess(1:nf,ia);
        d = xtrue - ils.xhat12;  d(:,ia) = fcndrpyd(ils.xhat12(:,ia),xtrue(:,ia));

        mc_xhat12(:,:,j) = d;
        mc_xhatsigma12(:,:,j) = ils.xhatsigma12;
        mc_xhattp(:,:,j) = a.ipned - ils.xhattp;
        mc_xhattpsigma(:,:,j) = ils.xhattpsigma;
        mc_initRMSE(j,:) = ils.initRMSE;
        mc_time(j,:) = [ils.t.init ils.t.ils];
        vs(j) = ils.successflag;
        %v1 = mc_initRMSE(:,1)~=0; disp([mean3(mc_initRMSE(v1,1:3)), mean3(mc_initRMSE(v1,4:6))]);
        fprintf('Finished MC %.0f/%.0f for tsvi %.0f/%.0f (tsv(i)=%g)\n',j,nmc,i,nts,tsv(i))
    end
    scale = [ones(1,6)*1E3 ones(1,6)];
    ts.sigma12mean(i,:) = mean(mean(mc_xhatsigma12(:,:,vs),3),1).*scale; %ILS covariance matrix
    ts.xhat12std(i,:) = mean(std(mc_xhat12(:,:,vs),[],3),1).*scale; %ILS error
    ts.sigma3mean(i,:) = mean(mean(mc_xhattpsigma(:,:,vs),3),1)*1000;
    ts.xhat3std(i,:) = mean(std(mc_xhattp(:,:,vs),[],3),1)*1000;
    
    ts.sigma12std(i,:) = mean(std(mc_xhatsigma12(:,:,vs),[],3),1).*scale;
    ts.xhat12mean(i,:) = mean(mean(mc_xhat12(:,:,vs),3),1).*scale;
    ts.sigma3std(i,:) = mean(std(mc_xhattpsigma(:,:,vs),[],3),1)*1000; %m
    ts.xhat3mean(i,:) = mean(mean(mc_xhattp(:,:,vs),3),1)*1000;
    
    ts.initRMSEmean(i,:) = mean(mc_initRMSE(vs,:));
    ts.initRMSEstd(i,:) = std(mc_initRMSE(vs,:));
    ts.timeperMC(i,:) = mean(mc_time(vs,:)); %seconds per run
    ts.ilssuccessfraction(i) = sum(vs)/nmc;
end
if matlabpool('size')>0;  matlabpool close;  end
elapsed = etime(clock,startclock1);
fprintf('\n\n%.0f MC Runs Completed in %.1fhrs (%.0fs), %.3fs per run.\n',nmc*nts,elapsed/3600,elapsed,elapsed/nmc/nts)
%save(fname,'ts','tsv','nmc','xstr')
if i==1;  v1 = find(fcnrange(mc_initRMSE(:,1:3))<30000 & mc_initRMSE(:,1)~=0); fprintf('%.2fm, %.2fdeg RMSE\n',mean3(mc_initRMSE(v1,1:3)),mean3(mc_initRMSE(v1,4:6)));  break;  end

hf = fcnplotTS(ts,tsv,nmc,xstr,cam0);  %hgsave(hf, fname)

