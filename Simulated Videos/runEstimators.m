function runEstimators()
clc; close all; ils=[]; %#ok<NASGU>
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',0));

%LOAD VIDEO ---------------------------------------------------------------
[file, path] = uigetfile({'*.mat','MAT Files (*.mat)';'*.*','All files (*.*)'},'Select video file','MultiSelect','on');  
if isnumeric(path); return; end; pathnow=what; addpath(pathnow.path); addpath(fcnpathm1(path)); addpath(path);  cd(path);
if ~iscell(file); x=file; clear file; file{1} = x; end
e = [];  for i=1:numel(file);  e1=run1file(file{i},path,i);  e=[e; e1]; end %#ok<AGROW>

% dx = e;
% r.x = dx;
% r.mean = mean(dx);
% r.std = std(dx);
% r.rms = [fcnrms(dx(:,1)) fcnrms(dx(:,2)) fcnrms(dx(:,3))];
% r
% fcnrms(dx)

% fig;
% c = {'.r','.g','.b','.c','.m'};
% str = {'SGC55','SGC69','S4','S5','S6','S9'};
% for i=1:size(dx,1)/6
%     j = (1:6) + (i-1)*6;
%     h(i)=fcnplot3(dx(j,:),c{i}); set(h(i),'markersize',20);
%     for k=1:6
%        ht1=text(dx(j(k),1),dx(j(k),2),dx(j(k),3),['  ' str{k}]); set(ht1,'Color',get(h(i),'Color'),'FontSize',8);
%     end
% end
% box on; fcn3label('x (m)','y (m)','z (m)'); set(gca,'zdir','reverse')
% legend(h,'1-100','100-200','200-300','300-400','400-500')
% plot3(0,0,0,'b.'); fcnplotline([0 0 0]); axis tight equal vis3d
end

function err=run1file(filename,pathname,ii)
fprintf('Beginning estimation for ''%s'' ...\n',filename);  load([pathname filename]);
ils=[]; err=[]; ils.xhattp=[]; cam.tpnedhat=[];
c = {'r','g','b','c','m'};

%OPTIMIZER OPTIONS ----------------------------------------------------
cam.optim.options0 = optimset('display','iter');
cam.optim.options1 = optimset('GradObj','off','display','off','largescale','off','MaxFunEvals',3000);
cam.optim.options2 = optimset('GradObj','off','display','off','MaxFunEvals',2000,'HessUpdate','bfgs','Algorithm','active-set');
cam.optim.options3 = optimset('Algorithm','levenberg-marquardt','display','off','TolX',1E-24,'Tolfun',1E-24,'MaxFunEvals',2E5,'MaxIter',1000);
cam.optim.options4 = optimset('Algorithm','levenberg-marquardt','display','off','Jacobian','on','DerivativeCheck','off','gradobj','on','TolX',1E-24,'Tolfun',1E-24,'MaxFunEvals',2E5,'MaxIter',1000);

%SCENARIO SPECIFIC --------------------------------------------------------
cam.xppo = 0;  cam.yppo = 0;
if regexpi(cam.filename,'synthetic') %GE SYNTHETIC CAMERA TRAJECTORY
    cam.focalLength = (cam.width/2) / tan(cam.fovh/2*d2r);   %pixels (0.00942 = pixel pitch in mm/pix)
    cam.sigmas.K = [3 3 3]; %[xppo yppo focalLength], pixels
elseif regexpi(cam.filename,'OSU') %OSU INDIVIDUAL
    cam.focalLength = 135/0.00942;  cam.widthOriginal = 2672;  cam.heightOriginal = 4016;
    cam.fovh = atand(cam.widthOriginal/2/cam.focalLength)*2; %deg
    cam.focalLength = (cam.width/2) / tan(cam.fovh/2*d2r); %pixels
    cam.sigmas.K = [5 10 50]; %[xppo yppo focalLength], pixels
elseif ~isempty(regexpi(cam.filename,'camera0')) %PURDUE INDIVIDUAL
    cam.focalLength = 100/0.0074;  cam.widthOriginal = 3248;  cam.heightOriginal = 4872;
    cam.fovh = atand(cam.widthOriginal/2/cam.focalLength)*2; %deg
    cam.focalLength = (cam.width/2) / tan(cam.fovh/2*d2r); %pixels
    cam.sigmas.K = [10 10 100]; %[xppo yppo focalLength], pixels
elseif ~isempty(regexpi(cam.filename,'purdue')) || ~isempty(regexpi(cam.filename,'2000x2000')) %PURDUE STITCHED
    cam.focalLength = 105/0.0074; %at 2k x 2k in pixels (0.0074 = pixel pitch in mm/pix)
    cam.xppo = 0.3/.0074;  cam.yppo = 21.5/.0074;
    cam.sigmas.K = [1000 1000 1000]/100; %[xppo yppo focalLength], pixels
elseif ~isempty(regexpi(cam.filename,'RAWvideo'))  %SHANKS THERMAL
    cam.fovh = 30; %deg
    cam.focalLength = (cam.width/2) / tan(cam.fovh/2*d2r);   %pixels (0.00942 = pixel pitch in mm/pix)
    cam.sigmas.K = [1 1 1]; %[xppo yppo focalLength], pixels
else
    fprintf('\nERROR: PLEASE SPECIFY SCENARIO SPECIFICS IN runEstimators.m\n\n')
    return
end
cam.apriori.K.xppo=cam.xppo;  cam.apriori.K.yppo=cam.yppo;  cam.apriori.K.focalLength=cam.focalLength;

%SMOOTH GPS TAGS ----------------------------------------------------------
%figure; subplot(211); plot(cam.true.dt,'.-'); title(sprintf('\\Deltat  ''%s''',cam.filename)); ylabel('s'); xlabel('frame'); axis tight
%subplot(212); plot(fcnrange(diff(cam.apriori.ned)),'.-'); title('\Deltax'); ylabel('m'); xlabel('frame'); axis tight

%COMMON -------------------------------------------------------------------
cam.msvgridpoints = 20;
ekf.zsigma = 3; %pixel uncertainty
cam.rpy1s = 20; %#ok<*PFBNS> %2deg
cam.gpsxy1s = 5; %5m true
cam.gpsz1s = 10; %10m true
cam.ilsxy1s = 2; %5m  ils assumed
cam.ilsz1s = 2; %10m ils assumed
cam.ilsrpy1s = 3;  %5deg, new 15deg
cam.ilsTrpy = .5; %10s, new 2s
cam.google.startdate = now;

cam.fovh = atand(cam.width/2/cam.focalLength)*2; %deg
cam.fovv = atand(cam.height/2/cam.focalLength)*2; %deg
cam.true.focus.pixel = camsc2pixel(cam, [1 0 0]);


%SELECT TIE POINTS --------------------------------------------------------
i = 1:numel(a.iis);
%mu = mean(a.score);  sigma = std(a.score);  i = find(a.score>mu-2*sigma);
ntp = min(100,numel(a.iis));  i = 1:ntp;  [~, j] = sort(a.score,'descend');  i=j(i);
%n=sum(a.iis);  i = ceil(rand(1,500)*n);
%i = find(a.state(:,1)==1 & a.state(:,end)==1); %persistent tie points
a = fcncropTP(cam,a,i);


%SELECT FRAMES ------------------------------------------------------------
%i = unique(round( linspace(1,cam.frames,min(10,cam.frames)) ));
%i = [1:500];
%[cam,a] = fcncropimages(cam,a,i);


% %MAKE PERFECT MEASUREMENTS ------------------------------------------------
% [ntp, nf] = size(a.upx);
% tpned = a.ipned + randn(size(a.ipned))*.005;
% for i=1:nf
%     f =  ned2pixel(cam,cam.frameID(i),tpned,'true') + randn(ntp,2)*1;
%     a.upx(:,i) = f(:,1);
%     a.upy(:,i) = f(:,2);
% end

%MSV ----------------------------------------------------------------------
tic; cam = fcnMSV(cam,a); toc; %fcnValidateResults(cam,a);

% %SFM ----------------------------------------------------------------------
%ha=fig(5,2);
%[~, hs] = fcnSFM(cam,a,'horn sequential');   %popoutsubplot(hs(2),ha(1)); popoutsubplot(hs(2),ha(2)); popoutsubplot(hs(1),ha(3)); popoutsubplot(hs(3),ha(4)); close(gcf); drawnow
[camsfm, hs] = fcnSFM(cam,a,'horn non-sequential');  %popoutsubplot(hs(1),ha(5)); popoutsubplot(hs(3),ha(6)); close(gcf); drawnow
%[~, hs] = fcnSFM(cam,a,'zisserman sequential');  popoutsubplot(hs(1),ha(7)); popoutsubplot(hs(3),ha(8)); close(gcf); drawnow
%[~, hs] = fcnSFM(cam,a,'zisserman non-sequential');  %popoutsubplot(hs(1),ha(9)); popoutsubplot(hs(3),ha(10)); close(gcf)
%view(ha(1),[-111 22]); legend(ha(1),'off'); axis(ha(1),'tight'); box(ha(1),'on'); set(ha([3 5 7 9]),'CameraViewAngle',8.75)
%fcnfontsize(10)
%export_fig(gcf,'-q95','-r150','-a4','TS1pix.jpg','-native')

% ha = fig(1,2,2);
% popoutsubplot(hs(1),ha(1)); popoutsubplot(hs(3),ha(2));  view(ha(2),90,90); 
% fcnfontsize(14); fcnmarkersize(ha(2),12); fcnmarkersize(ha(1),50); %fcnlinewidth(ha(2),2)
% export_fig(gcf,'-q90','-r150','-a4','exported.jpg','-native')
% 
% ha = fig(2,1,1.1,1);
% popoutsubplot(hs(2),ha(1)); popoutsubplot(hs(3),ha(2));  view(ha(1),90,90);  view(ha(2),90,90); 
% axis(ha,'off'); legend(ha(1),'off'); legend(ha(2),'off'); fcnfontsize(14)
% export_fig(gcf,'-q90','-r150','-a4','exported.jpg','-native')

%NLS ----------------------------------------------------------------------
[ils, cam] = fcnLM(cam,ekf,a);  if ~cam.syntheticVideoFlag; fcnValidateResults(cam,a); end

%j=ils.index.ned(1,:); C0 = ils.xhatC(j,j)*1E6
%j=ils.index.tp(1,:); C1 = ils.xhatC(j,j)*1E6

%[ils, cam] = fcnstartILS(cam,ekf,a);  fcnValidateResults(cam,a);
%[camsfm, hs] = fcnSFM(cam,a,'horn non-sequential');  %popoutsubplot(hs(1),ha(5)); popoutsubplot(hs(3),ha(6)); close(gcf); drawnow

%j=[1 2 3]; C3 = ils.Sigma(j,j)*1E6
%C4 = ils.xhat12Sigma([1 3 5],[1 3 5])*1E6


%EKF ----------------------------------------------------------------------
%[cam] = fcnEKF(cam,ekf,a);  fcnValidateResults(cam,a);


%STATISTICS ---------------------------------------------------------------
% n=sum(a.iis);
% a0 = a;
% cam0 = cam;
% matlabpool close force local
% matlabpool 4
% parfor j=1:1000
%     i = ceil(rand(1,500)*n);  b = fcncropTP(cam0,a0,i);  cam1=cam0;  cam1.tpnedhat=cam0.tpnedhat(i,:);
%     [ils, cam2] = fcnLM(cam1,ekf,b);
%     rmse(j)=fcnrms(fcnrange(fcnValidateResults(cam2,b)));
%     pr(j) = mean(ils.pixelresiduals);
% end
% matlabpool close force local
% 
% if ii==1; fig; end 
% x = [rmse' pr'];
% C = cov(x);  mu = mean(x);
% [ellipse, cross] = fcnerrorellipse( C , mu, .683);
% plot(rmse,pr,[c{ii} '.']); hold on;
% plot(ellipse.x,ellipse.y,[c{ii} '-']);
% plot(cross.x,cross.y,[c{ii} '-']);
% xlabel('Tie Point RMSE (m)'); ylabel('Mean Residuals (pixels)'); title(sprintf('RMSE = %.1f+/-%.1fm, Residuals = %.2f+/-%.2f pixels',mu(1),sqrt(C(1,1)),mu(2),sqrt(C(2,2)))); 
% % h=findobj(gca,'Marker','.');
% % legend(h([5 4 3 2 1]),'Segment 1','Segment 2','Segment 3','Segment 4','Segment 5');
% % set(h,'MarkerSize',4);

% % for j=1:3
% %      x=ils.pixelresiduals';  mu = mean(x);  sigma = std(x);  i = find(x<(mu+2*sigma));  a = fcncropTP(cam,a,i);  %i = find(x<(mu+2*sigma) | (a.state(:,1)==1 & a.state(:,end)==1));
% %      [ils, cam] = fcnLM(cam,ekf,a);  fcnValidateResults(cam,a);
% % end


%GOOGLE EARTH PLOTTING ----------------------------------------------------
[cam] = fcngetcamcornersLLA(cam,cam.tpnedhat,ils.xhattp);

%SAVE ---------------------------------------------------------------------
save([pathname cam.filename ' MSV+ILS.mat'],'a','cam','ekf','ils')
end
