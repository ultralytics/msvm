% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function C0 = fcnMIGLS(cam,a)
clc; close all
if nargin==0
    [file, path] = uigetfile({'*.mat','MAT Files (*.mat)';'*.*','All files (*.*)'},'Select video file','MultiSelect','on');  
    if isnumeric(path); return; end; pathnow=what; addpath(pathnow.path); addpath(fcnpathm1(path)); addpath(path);  cd(path);
    load([path file]);
end

a0=a;
cam0=cam;
c = {'r','g','b','c','m','y','k','w'};

%GET PIXEL VALUES
[vf(1), vf(2), vf(3)] = fcnget3distantframes(cam.aposteriori.ned);
nf=numel(vf);

[filename, pathname] = uigetfile({'*.avi;*.mpeg;*.mp4;*.wmv','Video Files (*.avi, *.mpeg, *.mp4, *.wmv)';'*.*','All files (*.*)'},'Select video file');  if filename==0; return; end;  pathnow=what; addpath(pathnow.path); addpath(fcnpathm1(pathname)); addpath(pathname);
vfr=VideoReader(filename);
hf1=figure;
hf2=figure;
mx=zeros(1,nf);  my=mx;  ntp = 10;
for i=1:nf
    j = cam.frameID(vf(i));
    f1.rgb = single(read(vfr,j))/255;  sf = size(f1.rgb);
    figure(hf1); cla; imshow(f1.rgb); title(sprintf('Frame %.0f, Pick a Point',j)); axis on
    
    ws = 50;
    for k=1:ntp
        [mxt, myt, button] = ginput(1);  mxt=round(mxt);  myt=round(myt);
        if button~=1; ntp=k-1;  break;  end
        xv = max(mxt-ws,1) : min(mxt+ws,sf(1));
        yv = max(myt-ws,1) : min(myt+ws,sf(1));
        [~, hf3] = fig;  imshow(f1.rgb(yv,xv,:),'xdata',xv,'ydata',yv); axis on;
        [mx(k,i), my(k,i)] = ginput(1);  close(hf3)
    end
    
%     figure(hf2); subplot(1,nf,i); imshow(f1.rgb); title([sprintf('Frame %.0f Points',j) sprintf('\n%.1f,  %.1f',mx(:,i),my(:,i))]); axis on; hold on
%     for k=1:ntp
%         plot(mx(k,i),my(k,i),[c{k} '+'],'Markersize',20)
%     end
end
close(hf1); drawnow
a.upx = mx;  a.upy = my;  a.state = ones(ntp,nf);

%VALIDATION POINTS
% load measuredLLA.mat
% vf = [1 50 100]; nf = 3;
% mx = reshape(gt.xy(1:18,1),[6 3]);
% my = reshape(gt.xy(1:18,2),[6 3]);  ntp=6;
% a.upx = mx;  a.upy = my;  a.state = ones(ntp,nf);


%MIG MSV
cam1 = cam;  cam1.frames = nf;
cam1.aposteriori.rpy = cam.aposteriori.rpy(vf,:);
cam1.aposteriori.ned = cam.aposteriori.ned(vf,:);
a.state = ones(ntp,nf);
[C0, msvr] = fcnMIGMSV(cam1,a);


%RUN LS
a=a0; cam=cam0; i=numel(a.score)+(1:ntp);
cam.tpnedhat(i,:) = C0;
a.upx(i,vf) = mx;
a.upy(i,vf) = my;
a.state(i,:) = 0;  a.state(i,vf) = 1;
a.score(i) = 1;
a.iis(i) = 1;

ekf.zsigma = ones(i(end),1)*3;
ekf.zsigma(i) = 1;
[ils, cam] = fcnLM(cam,ekf,a);
err = fcnValidateResults(cam,a);

%PLOT MATLAB
str=[];
for j=i    
    mu = ils.xhattp(j,:);
    C = ils.xhatC(ils.index.tp(j,:),ils.index.tp(j,:));
    sigma = ils.xhattpsigma(j,:)*1000;

    [ellipse, cross] = fcnerrorellipse(C,mu,.9);
    % fig; surf(ellipse.x,ellipse.y,ellipse.z,'FaceColor',[0 0 0],'FaceAlpha',.5,'FaceLighting','gouraud','EdgeColor','none');%,'AmbientStrength',1E6);
    % box on; axis equal vis3d; hold on
    % plot3(cross.x,cross.y,cross.z,'b-','linewidth',1.2);
    
    %PLOT GOOGLE EARTH
    s = size(ellipse.x);
    lla = lla2llag(ned2lla(cam.DEM,[ellipse.x(:) ellipse.y(:) ellipse.z(:)]));
    x = reshape(lla(:,1),s);
    y = reshape(lla(:,2),s);
    z = reshape(lla(:,3),s);
    str1 = ge_plot3(y,x,z,'lineColor',[dec2hex(round(.4*255)) fcnhexcolor([0 1 0])],'lineWidth',1,'altitudeMode','absolute','name','ellipse');
    
    s = size(cross.x);
    lla = lla2llag(ned2lla(cam.DEM,[cross.x(:) cross.y(:) cross.z(:)]));
    x = reshape(lla(:,1),s);
    y = reshape(lla(:,2),s);
    z = reshape(lla(:,3),s);
    str2 = ge_plot3(y,x,z,'lineColor',[dec2hex(round(.8*255)) fcnhexcolor([0 1 0])],'lineWidth',2,'altitudeMode','absolute','name','cross');
    
    str = [str ge_folder(sprintf('Tie Point %d',j),[str1 str2])]; %#ok<AGROW>
    
    fprintf('\nManually Selected Point TP%.0f:\nTie Point Position 1sigma (m NED) = [%5.1f %5.1f %5.1f]',j,sigma)
    fprintf('\nCovariance Matrix (m NED) = [%5.1f %5.1f %5.1f',C(1,:)*1E6)
    fprintf('\n                             %5.1f %5.1f %5.1f',C(2,:)*1E6)
    fprintf('\n                             %5.1f %5.1f %5.1f ]\n',C(3,:)*1E6)
    
end
cam.google.kml.manualPoints = ge_folder('Manual Points',str);

%SAVE TO KML FILE ---------------------------------------------------------
fields = fieldnames(cam.google.kml);  str=[];
for i=1:numel(fields);  str = [str ' cam.google.kml.' fields{i}];  end;  str = ['[' str ']'];
ge_output([cam.filename ' MIGLS.kmz'],eval(str));
end

