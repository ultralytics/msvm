% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [] = fcnimportvideo()
clc; close all; cam=[]; a=[]; fcnaddmfilepath; reset(RandStream.getGlobalStream);
flags.video = true;
flags.plots = true;

%LOAD VIDEO --------------------------------------------------------------- %example: 'GE Synthetic Video - Madrid Buildings 73f.avi'
adir = dir;
i = cell2mat(arrayfun(@(x) (any(regexpi(x.name,'.mp4')) || any(regexpi(x.name,'.avi')) || any(regexpi(x.name,'.mov'))) && ~any(regexpi(x.name,'Features')),adir,'UniformOutput', false));
if sum(i)==1 %only one possible video file in this dir
    fname = adir(i).name;
    pname = [cd filesep];
else
    if any(fcnlastfileloaded)
        [fname, pname] = uigetfile(fcnlastfileloaded);
    else
        [fname, pname] = uigetfile({'*.avi;*.mpeg;*.mp4;*.mov','Video Files (*.avi, *.mpeg, *.mp4, *.mov)'},'Select video file','MultiSelect','on');
    end
    if fname==0; return; end;  pathnow=what; addpath(pathnow.path); addpath(fcnpathm1(pname)); addpath(pname);  cd(pname);
end
pf=[pname fname];  fcnlastfileloaded(pf);
fprintf('Beginning feature tracking for ''%s'' ...\n',fname);

%filename = 'GE Synthetic - Vesuvius 1280x720 45FOV.mp4';
%pathname = fcnpathm1(which(filename));

%LOAD METADATA ------------------------------------------------------------ %example: 'GE Synthetic Video - Madrid Buildings 73f metadata.mat'
ma = [pname fname(1:end-4) '.mat']; %prediction a
mb = [pname 'llat.mat']; %prediction b
a=dir(pname); a={a.name}; mc=a{contains(a,fname(1:end-4)) & contains(a,'metadata.mat')};
if exist(ma,'file')
    matfname=ma; matpname=pname;
elseif exist(mb,'file')
    matfname=mb; matpname=pname;
elseif exist(mc,'file')
    matfname=mc; matpname=pname;
else %user picks a matfile
    [matfname, matpname] = uigetfile({'*.mat',pname},['Select metadata *.mat file for ''' fname '''']);  if matfname==0; return; end;
end
fprintf('Attaching metadata from ''%s'' ...\n',matfname);  load([matpname matfname]);  startclock=clock;
V.vfr = VideoReader(fname);  nf=V.vfr.Duration*V.vfr.FrameRate;   cam.frames=nf;  cam.width=V.vfr.Width;  cam.height=V.vfr.Height;  cam.fps=V.vfr.FrameRate;  cam.xppo=0;  cam.yppo=0;
%V.vvfr = vision.VideoFileReader(fname,'VideoOutputDataType','inherit');

%cam.frameID = 1:3:500;% nf; %purdue
cam.frameID =   1:1:min(105,nf); %synthetic
%cam.frameID =  1:min(69,nf); %OSU 1
%cam.frameID =  90:1:min(174,nf); %OSU 2

if  ~isempty(regexpi(fname,'synthetic')) %GE Synthetic Video
    cam.syntheticVideoFlag = true;
    cam.focalLength = (cam.width/2) / tan(cam.fovh/2*d2r);  %pixels (0.00942 = pixel pitch in mm/pix)
    cam.fovv = atan(cam.height/2/cam.focalLength)*r2d*2;  %deg
    cam.fps = mean(1./cam.true.dt);
    lla = [cam.true.lla cam.true.t'];
    llaID = lla(cam.frameID,1:3);  nf = numel(llaID(:,1));  cam.frames=nf;  zva = zeros(nf,3);  cam.DEM = getDEM(cam,llaID);
elseif ~isfield(cam,'syntheticVideoFlag') || ~cam.syntheticVideoFlag %Real Video
    cam.syntheticVideoFlag = false;
    llaID = lla(cam.frameID,1:3);  nf = numel(llaID(:,1));  cam.frames=nf;  zva = zeros(nf,3);  cam.DEM = getDEM(cam,llaID); %#ok<*NODEF>

    ecef = lla2ecef(llaID);
    ned = ecef2ned(cam.DEM,ecef);
    cam.apriori.lla = llaID;   cam.aposteriori.lla = llaID;   cam.true.lla = zva;
    cam.apriori.ecef = ecef;   cam.aposteriori.ecef = ecef;   cam.true.ecef = zva;
    cam.apriori.ned = ned;     cam.aposteriori.ned = ned;     cam.true.ned = zva;
    cam.apriori.rpy = zva;     cam.aposteriori.rpy = zva;     cam.true.rpy = zva;
    cam.true.t = lla(cam.frameID,4)'-lla(cam.frameID(1),4);   cam.true.dt = [0 diff(cam.true.t)];   cam.fps = 1/mean(cam.true.dt(2:end));
end
cam.google.startdate = now;
cam.filename = fname(1:end-4);
cam.pathname = pname;
%INITIALIZE VIDEO FILE READER ---------------------------------------------
ft.rgb=[];  ft.gray=[];  ft.SURF.points=[];  ft.SURF.features=[];  ft.SURF.vpts=[];  ft.valid=[];  ft.points=[];  ft.tform=[];  ft.score=[];  ft.state=[];  frames(1:cam.frames)=ft;  f1=ft;
[f1.rgb, f1.gray] = fcngetVideoReaderframe(V,cam.frameID(1));
[     ~, f2.gray] = fcngetVideoReaderframe(V,cam.frameID(2));  %find a frame about 3 second away

%SELECT INITIAL POINTS ----------------------------------------------------
bs = 11;  bs = max(round(cam.width/90),bs);  if ~fcnodd(bs); bs=bs+1; end;           cam.KLT.bs = bs; %11 block size of feature templates
ns = 13;  ns = max(round(cam.width/20),ns);  if ~fcnodd(ns); ns=ns+1; end;           cam.KLT.ns = ns; %9 neighborhood size (exclusion area around tiepoints)
edge = ceil(bs*1.5);                                    cam.KLT.edge = edge;
ROI = [edge edge cam.width-edge*2 cam.height-edge*2];   cam.KLT.ROI = ROI;
p = detectFeatures(f1.gray,2000,ROI);
%[features,valid_points] = extractFeatures(f1.gray,points,'BlockSize',bs);
%close; y=normxcorr2(reshape(double(features(100,:)),bs,bs),double(f1.gray)); imagesc(y.^3)

%INITIALIZE KLT TRACKERS --------------------------------------------------
cam.KLT.PTa = vision.PointTracker('MaxBidirectionalError',0.9,'NumPyramidLevels',1,'BlockSize',[bs bs],'MaxIterations',50);
cam.KLT.PTb = vision.PointTracker('MaxBidirectionalError',  2,'NumPyramidLevels',5,'BlockSize',[bs bs],'MaxIterations',50);
cam.KLT.PTc = vision.PointTracker('MaxBidirectionalError',  2,'NumPyramidLevels',5,'BlockSize',[bs bs],'MaxIterations',50);

%REMOVE BURN INS FROM INITIAL POINTS --------------------------------------
initialize(cam.KLT.PTc, p, f1.gray);  [p2, valid] = step(cam.KLT.PTc, f2.gray); p=p(valid,:);  p2=p2(valid,:);
[~, ~, inliers] = fcnimagexform(p,p2);  p=p(inliers,:);  %fig;  p2=p2(inliers,:);  imshow(f2.gray);  plot(p(:,1),p(:,2),'r.',p2(:,1),p2(:,2),'g.',[p(:,1) p2(:,1)]',[p(:,2) p2(:,2)]','y-')
initialize(cam.KLT.PTa, p, f1.gray);  %[p, valid] = step(cam.KLT.PTa, f1.gray);

%INITIALIZE VISION OBJECTS ------------------------------------------------
halphablend = vision.AlphaBlender('Operation','Binary mask','MaskSource','Input port');
mi = vision.MarkerInserter('Shape','Circle','Fill',true,'FillColor', 'Custom','Size',4,'Antialiasing',false);
rmi = clone(mi);  rmi.CustomFillColor = uint8([1   0   0]*255);
gmi = clone(mi);  gmi.CustomFillColor = uint8([.5  1  .5]*255);
bmi = clone(mi);  bmi.CustomFillColor = uint8([0   0   1]*255);
omi = clone(mi);  omi.CustomFillColor = uint8([1 .65   0]*255);
gli = vision.ShapeInserter('Shape','Lines','BorderColor', 'Custom', 'CustomBorderColor', uint8([.5  1  .5]*255), 'Antialiasing',false); %#ok<NASGU>

%START FRAME BY FRAME TRACKING --------------------------------------------
valid = true(size(p,1),1);
np=sum(valid);
tformchain = eye(3);
a.upx = nan(np,nf);      a.upx(:,1)=p(:,1);
a.upy = nan(np,nf);      a.upy(:,1)=p(:,2);
a.state = zeros(np,nf,'int8');  a.state(:,1)=1;
f1.points = p;  f1.valid = valid;  f1.state = a.state(:,1);  f1.score = zeros(np,1);

if flags.video
    %vid=VideoWriter([cam.filename ' Features'],'MPEG-4'); vid.FrameRate=30; vid.Quality=100; open(vid); %'Motion JPEG AVI'
    vvid = vision.VideoFileWriter([cam.filename ' Features.mp4'],'FrameRate',30,'FileFormat','MPEG4','Quality',100);
    f1base = f1.rgb;   f1ab = step(gmi, f1base, uint16(p));
end
[X,Y] = meshgrid(1:cam.width,1:cam.height);  cam.KLT.Ixy = single([X(:) Y(:) ones(numel(X),1)]);  clear p valid X Y
for i=1:numel(cam.frameID)
    ntls=0;  nels=0;
    if i>1
        if i>cam.frames;  i=i-1;  break;  end
        %f1=ft;  f1.valid = false(size(f0.valid));  f1.points = f0.points;  f1.score = f0.score;  f1.state = a.state(:,i-1);  a.state(:,i) = a.state(:,i-1);
        f1=f0;  f1.SURF=ft.SURF;  f1.valid(:)=false;  f1.tform=[];  a.state(:,i)=f1.state;
        
        %GET NEW RGB AND GRAY IMAGES
        [f1.rgb, f1.gray] = fcngetVideoReaderframe(V, cam.frameID(i));
        
        %BACK ONE FRAME
        f1=fcnKLT1(f0,f1,cam);
        
        %BACK TWO FRAME
        tprt = 0.90; %tie point retention threshold (anything above we don't go back and look)
        try
            if i>2  &&  sum(f0.valid)/sum(frames(i-2).valid) < tprt
                [f1]=fcnKLT2(frames,f0,f1,cam,a,i);
                nels = sum(f1.state==1 & f0.state==2); %number edge losses saved
                ntls = sum(f1.state==1 & f0.state==3); %number of tracking losses saved
            end
        catch
            'BACK TWO FRAMES NOT WORKING!'
        end
        
        %UPDATE STATES
        a.state(:,i) = f1.state;
        v1 = f0.valid | f1.valid;  a.upx(v1,i)=f1.points(v1,1);  a.upy(v1,i)=f1.points(v1,2);
        f1 = fcnDropTiePoints(f1,a,i,3); %drop after 3 unsuccessful frame
    end
    
    %WRITE VIDEO
    if flags.video
        if i>1
            videoWriteMode = 'overlayTracksDirectly'; %{'overlayTracksDirectly','overlaySyncedToFrame1'}
            switch videoWriteMode
                case 'overlayTracksDirectly'  %PLOT TRACKS DIRECTLY OVER RAW VIDEO
                    f1ab = step(gmi, f1.rgb, uint16(f1.points(f1.valid, :)));
                    i1 = max(1,i-10); %tie point tail
                    lines = zeros(sum(f1.valid),(i-i1+1)*2);  lines(:,1:2:end) = a.upx(f1.valid,i1:i);  lines(:,2:2:end) = a.upy(f1.valid,i1:i);  lines = lines(lines(:,1)~=0,:);
                    f1ab = step(gli, f1ab, uint16(lines));
                case 'overlaySyncedToFrame1'  %PLOT ALL TRACKS SYNCED TO FRAME 1
                    tformchain = f1.tform*tformchain;  tform = tformchain;  itform = tform^-1; %[tform, itform] = fcnimagexform(f1.points(f1.valid,:),frames(1).points(f1.valid,:));
                    fnw1 = fcnimwarp(f1.rgb,cam.KLT.Ixy,tform,itform);  xw = uint16(fcnxformpoints(f1.points,tform));
                    
                    mask = max(fnw1,[],3)~=0;
                    f1base = step(halphablend, f1base, fnw1, mask); %alpha blend base
                    
                    f1ab = step(gmi, f1base, xw(f1.valid,:));
                    f1ab = step(bmi, f1ab, xw(f1.state==2,:));
                    f1ab = step(omi, f1ab, xw(f1.state==3,:));
                    if i>2; f1ab = step(rmi, f1ab, xw(f1.state==0 & (f0.state==3 | f0.state==2), :)); end
                    %mask = repmat(~mask,[1 1 3]);  f1ab(mask) = 0.3 + 0.7*f1ab(mask);
            end
        end
        %writeVideo(vvid,f1ab);
        step(vvid, f1ab);
    end
    
    %PRINT FRAME RESULTS
    np  = sum(f1.state==1); %number alive
    nel = sum(f1.state==2); %number edge losses
    ntl = sum(f1.state==3); %number tracking losses
    fprintf('frame %3.0f/%3.0f, %3.0f surviving features (%3.0f edge losses, %3.0f saved, %3.0f tracking losses, %3.0f saved)\n',cam.frameID(i),cam.frameID(nf),np,nel,nels,ntl,ntls);
    
    %ADD MORE TIE POINTS AS NECESSARY
    tpt = 1250; %tie point threshold (any count below and we refresh)
    if i>1  &&  i<nf  &&  np<tpt  &&  (sum(f1.valid)/sum(f0.valid)>tprt || sum(f0.valid)<tpt)
        nadd = (tpt-np) + tpt; %number to add
        %nadd = 300;
        [a, f1] = fcnAddTiePoints(f1, a, i, nadd, cam);  
        np=sum(f1.state==1);  
    end
    
    %SAVE FRAME INFO
    if i>1 && i<nf; f1.rgb=[]; end
    frames(i) = f1;
    if i>5;  frames(i-4).gray = []; end %taking too much memory!
    f0 = f1;
    
    %EXIT
    if np==0;  i=max(i-1,1);  fprintf('---WARNING!--- Premature Termination, No Surviving Features Left!!\n'); break;  end
end
vf = 1:i;  cam.frameID = cam.frameID(vf);  a.score = f1.score./sum(a.state==1,2);  cam.KLT=[];
if i~=cam.frames
    a.upx = a.upx(:,vf);  a.upy = a.upy(:,vf);  a.state = a.state(:,vf);  cam.frames=i;  frames = frames(vf);
end
fprintf('Finished Feature Tracking for ''%s'' in %.1fs\n\n',fname,etime(clock,startclock))

%PLOT RESULTS -------------------------------------------------------------
if flags.plots
    fcnplotFeatureTracks(a, cam,frames(1).rgb, fcngetVideoReaderframe(V,i) ); drawnow
    hf = findobj(0,'Type','Figure');
    for i=1:numel(hf)
        %export_fig(hf(i),'-q95','-r200','-a1',sprintf('%s%s figure%g.jpg',pathname,cam.filename,i));
    end
end
cam = fcngetcamcornersLLA(cam);

%CUT SHORT TRACKS ---------------------------------------------------------
%f1v = f1.valid;  a.upx = a.upx(f1v,vf);  a.upy = a.upy(f1v,vf);  a.state = a.state(f1v,vf);  a.score = a.score(f1v);  %a.iis = all(a.state,2); %include in statistics
f1v = sum(a.state==1,2)>5;  %at least 3 valid frames
a.upx = a.upx(f1v,:);  a.upy = a.upy(f1v,:);  a.state = a.state(f1v,:);  a.score = a.score(f1v);  a.iis = ones(sum(f1v),1);


%ADD TRUE TP LOCATIONS ----------------------------------------------------
if cam.syntheticVideoFlag %if this is a GE Synthetic video and we have cam.true data
    
    ntp = size(a.state,1);
    a.ipned = zeros(ntp,3);
    for i=1:ntp
        j = find(a.state(i,:)==1,1,'first');
        a.ipned(i,:) = pixel2ecef(cam.DEM,cam,j,[a.upx(i,j) a.upy(i,j)],'true');
    end
    a.ipned = ecef2ned(cam.DEM, a.ipned);
    %a.ipned = ecef2ned(cam.DEM, pixel2ecef(cam.DEM,cam,1,[a.upx(:,1) a.upy(:,1)],'true'));
    
    i = cam.frameID; cam.frames = numel(i);
    if ~isfield(cam.apriori,'lla');     cam.apriori.lla =       ecef2lla(cam.apriori.ecef);     end
    if ~isfield(cam.aposteriori,'lla'); cam.aposteriori.lla =   ecef2lla(cam.aposteriori.ecef); end
    if ~isfield(cam.true,'lla');        cam.true.lla =          ecef2lla(cam.true.ecef);        end
    cam.apriori.lla = cam.apriori.lla(i,:);    cam.aposteriori.lla = cam.aposteriori.lla(i,:);    cam.true.lla = cam.true.lla(i,:);
    cam.apriori.ecef = cam.apriori.ecef(i,:);  cam.aposteriori.ecef = cam.aposteriori.ecef(i,:);  cam.true.ecef = cam.true.ecef(i,:);
    cam.apriori.ned = cam.apriori.ned(i,:);    cam.aposteriori.ned = cam.aposteriori.ned(i,:);    cam.true.ned = cam.true.ned(i,:);
    cam.apriori.rpy = cam.apriori.rpy(i,:);    cam.aposteriori.rpy = cam.aposteriori.rpy(i,:);    cam.true.rpy = cam.true.rpy(i,:);
    cam.true.t = lla(i,4)'-lla(i(1),4);   cam.true.dt = [0 diff(cam.true.t)];   cam.fps = 1/mean(cam.true.dt(2:end));
end

%FINISH -------------------------------------------------------------------
savefilename = [pname cam.filename ' metadata.mat'];
save(savefilename,'a','cam','lla');
evalin('base','clear')
end


function f1 = fcnDropTiePoints(f1,a,i,n)
%cut tie points after not spotting them in n frames
if i>n
    j = (i-n+1) : (i-1);
    i = f1.state~=1 & all(a.state(:,j)~=1,2);
    f1.state(i)=0;
    f1.valid(i)=0;
end
end


function [a, f1] = fcnAddTiePoints(f1, a, i, n, cam)
[x,p] = detectFeatures(f1.gray,1E4,cam.KLT.ROI);

ep = f1.points(f1.valid,:); %existing points
rs = min( bsxfun(@minus,ep(:,1),x(:,1)').^2 + bsxfun(@minus,ep(:,2),x(:,2)').^2 ); %min range squared

p = p(rs>(cam.KLT.ns)^2);  p = p.selectStrongest(n); p=p.Location;

fprintf('%d new points added of %d requested\n',numel(p)/2,n)
j = numel(f1.state) + (1:size(p,1));
f1.valid(j) = true;
f1.points(j,:) = p;
f1.score(j) = 0;
f1.state(j) = 1;
a.state(j,i) = 1;
a.upx(j,i) = p(:,1);  a.upy(j,i) = p(:,2);
end


function [x,p] = detectFeatures(I,n,ROI)
%https://la.mathworks.com/help/vision/local-feature-extraction.html
if nargin<3; ROI = [1 1 size(I)]; end
p = detectBRISKFeatures(I,'ROI',ROI);
%p = detectKAZEFeatures(I,'ROI',ROI); %NOT AVAILABLE in 2016B!
if nargin>1 && ~isempty(n);  p=p.selectStrongest(n);  end
x = p.Location;
end


function [f1]=fcnKLT1(f0,f1,cam)
PTa = cam.KLT.PTa;  PTb = cam.KLT.PTb;  f0v = find(f0.valid);

x0 = f0.points(f0v,:);  n = 500 - numel(f0v);
xr = [x0; cam.KLT.edge+[rand(n,1)*(cam.width-2*cam.KLT.edge), rand(n,1)*(cam.height-2*cam.KLT.edge)]];
release(PTb);  initialize(PTb, xr, f0.gray);  step(PTa, f0.gray);  setPoints(PTb, xr);  [x1, v1] = step(PTb, f1.gray);

if sum(v1)>200
    [tform, itform] = fcnimagexform(x1(v1,:),xr(v1,:));
else
    [tform, itform, f1, f0] = fcnSURFxform(f1,f0); %tform = f0 to f1
end
f0_1gray = fcnimwarp(f0.gray,cam.KLT.Ixy,itform,tform); %frame 0 aligned to frame 1
f1.points = fcnxformpoints(f0.points,itform);

nep = ~fcnFEP(cam,f1.points(f0v,:)); %not edge point
f0vnep = f0v(nep);
if ~isempty(f0vnep)
    x0_1 = f1.points(f0vnep,:);
    release(PTa);  initialize(PTa, x0_1, f0_1gray);  step(PTa, f0_1gray);  setPoints(PTa, x0_1);  [x1, v1, scores] = step(PTa, f1.gray);
    f1.points(f0vnep,:) = x1;
    f1.valid(f0vnep) = v1;
    f1.score(f0vnep) = f1.score(f0vnep) + scores;
end
f1.tform = tform;
edgepoints = false(size(f0.valid));  edgepoints(f0v) = ~nep;

f1.state(f1.valid) = 1;
f1.state(edgepoints) = 2;
f1.state(~f1.valid & f0.valid & ~edgepoints) = 3; %tracking loss

% fig; imshow(f0.gray); plot(x0(:,1),x0(:,2),'go',x0(v1,1),x0(v1,2),'ro','Markersize',10);
% fig; imshow(f0_1gray);  plot(x0_1(:,1),x0_1(:,2),'go','Markersize',10);
% fig; imshow(f1.gray);   plot(x1(v1,1),x1(v1,2),'go','Markersize',10);
% fig; imshowpair(f1.gray,f0_1gray,'falsecolor','colorchannels',[2 1 2]);   plot(x0_1(:,1),x0_1(:,2),'r.',x1(v1,1),x1(v1,2),'go','Markersize',10);
% fig; imshowpair(f1.gray,f0_1gray);   plot(x0_1(:,1),x0_1(:,2),'r.',x1(v1,1),x1(v1,2),'go','Markersize',10);
% c=f1.score(f0vnep(v1));  scatter(x1(v1,1),x1(v1,2),40,c,'filled'); set(gca,'clim',fcnminmax(c)+[0 1E-6]); colorbar;
end


function [f1]=fcnKLT2(frames,f0,f1,cam,a,i)
ep = f0.state==2; %edge points
lp = f0.state==3; %lost points
f0 = frames(i-2);
f0.valid = (lp | ep) & a.state(:,i-2)==1;
PTa = cam.KLT.PTa;  PTb = cam.KLT.PTb;  f0v = find(f0.valid);
if numel(f0v)==0; return; end

%x0 = f0.points(f0v,:);  n = 500 - numel(f0v);
%xr = [x0; cam.KLT.edge+[rand(n,1)*(cam.width-2*cam.KLT.edge), rand(n,1)*(cam.height-2*cam.KLT.edge)]];
%release(PTb);  initialize(PTb, xr, f0.gray);  step(PTa, f0.gray);  setPoints(PTb, xr);  [x1, v1] = step(PTb, f1.gray);

v1 = f1.valid;
if sum(v1)>100
    [tform, itform] = fcnimagexform(f1.points(v1,:),f0.points(v1,:));
else
    [tform, itform, f1, f0] = fcnSURFxform(f1,f0); %tform = f0 to f1
end
%C = fcntformchain(f1,frames,i-2,i);

f0_1gray = fcnimwarp(f0.gray,cam.KLT.Ixy,itform,tform); %frame 0 aligned to frame 1

tf = projective2d(itform);
Rin = imref2d(size(f0.gray));
%Rin.XWorldLimits = Rin.XWorldLimits-mean(Rin.XWorldLimits);
%Rin.YWorldLimits = Rin.YWorldLimits-mean(Rin.YWorldLimits);
out = imwarp(f0.gray,Rin,tf);



f1.points(f0v,:) = fcnxformpoints(f0.points(f0v,:),itform);

nep = ~fcnFEP(cam,f1.points(f0v,:)); %not edge point
f0vnep = f0v(nep);  if isempty(f0vnep); return; end
x0_1 = f1.points(f0vnep,:);
release(PTa);  initialize(PTa, x0_1, f0_1gray);  step(PTa, f0_1gray);  setPoints(PTa, x0_1);  [x1, v1, scores] = step(PTa, f1.gray);  if sum(v1)==0; return; end

f1.points(f0vnep,:) = x1;
f1.valid(f0vnep) = v1;
f1.score(f0vnep) = f1.score(f0vnep) + scores;
%f1.tform = tform;
%edgepoints = false(size(f0.valid));  edgepoints(f0v) = ~nep;

f1.state(f1.valid) = 1;
%f1.state(edgepoints) = 2;
%f1.state(~f1.valid & f0.state==3) = 3; %tracking loss
%fprintf('%.0f saved, ',sum(v1))
end


function C = fcntformchain(f1,frames,i,j)
frames(j).tform = f1.tform;
C = eye(3);
for k = (i+1):j
    C = C*frames(k).tform;
end
end


function edgepoints = fcnFEP(cam,x,valid) %Find Edge Points
edge = cam.KLT.edge;
if nargin==3
    edgepoints = false(size(valid));
    edgepoints(valid) = x(:,1)<edge | x(:,2)<edge | x(:,1)>(cam.width-edge) | x(:,2)>(cam.height-edge);
else
    edgepoints = x(:,1)<edge | x(:,2)<edge | x(:,1)>(cam.width-edge) | x(:,2)>(cam.height-edge);
end
end


function y = fcnxformpoints(x,tform)
y = [x ones(numel(x)/2,1)]*tform(:,1:2);
end


function [tform, itform, f1, f0] = fcnSURFxform(f1,f0) %f1=current frame, f0=previous frame
ndsp = 1000; %number desired surf points
st = 2500; %SURF threshold

for i=0:5
    if isempty(f0.SURF.points)
        f0.SURF.points = detectSURFFeatures(f0.gray,'MetricThreshold',st/2^i);
        [f0.SURF.features, f0.SURF.vpts] = extractFeatures(f0.gray, f0.SURF.points.selectStrongest(ndsp),'SURFSize',64);
    end
    f1.SURF.points = detectSURFFeatures(f1.gray,'MetricThreshold',st/2^i);
    [f1.SURF.features, f1.SURF.vpts] = extractFeatures(f1.gray, f1.SURF.points.selectStrongest(ndsp),'SURFSize',64);
    if f0.SURF.points.Count>250 && f1.SURF.points.Count>250; break; end
end
index_pairs = fcnmatchSURF(f1.SURF.features,f0.SURF.features,.3); %index_pairs = matchFeatures(features1, features0,'Prenormalized',true,'Method','NearestNeighborRatio','MaxRatio',.3);    % SURF feature vectors are already normalized.
mp0 = f0.SURF.vpts(index_pairs(:, 2)); %matched points in frame 0
mp1 = f1.SURF.vpts(index_pairs(:, 1)); %matched points in frame 1
[tform, itform] = fcnimagexform(mp1.Location, mp0.Location);
% close all; fig;  plot(mp0.Location(:,1),mp0.Location(:,2),'b*');  plot(mp1.Location(:,1),mp1.Location(:,2),'r*')
% plot([mp0.Location(:,1) mp1.Location(:,1)]',[mp0.Location(:,2) mp1.Location(:,2)]','g')
% x = [mp0.Location ones(mp0.Count,1)]*itform; x=x(:,1:2);  plot(x(:,1),x(:,2),'bo'); drawnow %blue circles should be centered on red astericks
end

function [tform, itform, inliers] = fcnimagexform(xin,xout)
%xin = [nx2] pixels in source image
%xout = [nx2] pixels in output image
%tform = [3x3] converts from input to output image [xout 1] = [xin 1]*tform;

n = 4;
srl = 2.8;
np0 = 0;  inliers = true(size(xin,1),1);
for i=1:n
    np = sum(inliers);  if np==np0; break; end
    ov = ones(np,1);
    z = [xout(inliers,:) ov];  H = [xin(inliers,:) ov];
    tform = (H'*H)\H'*z; %LLS
    
    %TIE VECTOR METRICS
    %dx = xout-xin;
    %dxr = sqrt(dx(:,1).^2 + dx(:,2).^2); %dx range
    %dxa = fcnatan2(dx(:,2),dx(:,1)); %dx angle
    %[~, inliersa] = fcnsigmarejection(dxr,srl,3,'onlyIndices');
    %[~, inliersb] = fcnsigmarejection(dxa,srl,3,'onlyIndices');
    
    
    %RESIDUAL METRICS
    if i<n
        res = z - H*tform; %xyz residual
        r = sqrt(sum(res.^2,2)); %range residual
        [~, inliersc] = fcnsigmarejection(r,srl,3,'onlyIndices');
        inliers(inliers) = inliersc;
    end
    np0=np;
end
itform = tform^-1;  %itform = (z'*z)\z'*H; %LLS
end


function J = fcnimwarp(I,Ixy,itform,tform)
p1 = Ixy*tform(:,1);  p2 = Ixy*tform(:,2);  %faster than p=Ixy*tform(:,1:2)

%J2=uint8(reshape(interp2(single(I),p1,p2),size(I)));

if isinteger(I)
    J = interp2mexchar(I,p1,p2);
else
    J = interp2mexsingle(I,p1,p2);
end
J = reshape(J, size(I));  %fig; imshowpair(I1,I2);
end


function [rgb, gray] = fcngetVideoReaderframe(V,j)
%rgb = step(V.vvfr); %method 1
%rgb = read(V.vfr,j); %method 2
V.vfr.CurrentTime = (j-1)/V.vfr.FrameRate; 
rgb = readFrame(V.vfr); %method 3
gray = rgb2gray(rgb);  %    gray = rgb2graymex(rgb); %mex is slower unfortunately
gray = matchintensity(gray);
%gray = imadjust(gray);
end

function x = matchintensity(x)
%x0 = x;
pdf = raylpdf(linspace(0,1,256),.3);  
i = x>1 & x<254;  
x(i) = histeq(x(i),pdf);  %fig; imshowpair(x0,x,'montage')
%x = adapthisteq(x0,'clipLimit',0.02,'Distribution','rayleigh','NBins',256,'NumTiles',[4 4]);
end