function fcnCreateSyntheticVideo()
%CREATES SYNTHETIC GE VIDEO, CIRCULAR ORBIT OR BASED ON EXISTING GPS
clc; close all; tic
%reset(RandStream.getGlobalStream);
ss = get(0,'ScreenSize');
cam.ssge = [8 ss(4)-720-77 1280 720];  %1080x720
cam.syntheticVideoFlag = true;

useLLATandLLAFOCUSflag = 0; %i.e. use real world GPS points, i.e. OSU
if useLLATandLLAFOCUSflag
    load('llafocus'); %#ok<UNRCH>
    load('llat');  nf = numel(lla(:,1));
    cam.DEM = getDEM(cam,lla(1:nf,1:3),'geoid');
else
    r = 2; %km radius from llafocus
    h = 1; %km height above llafocus
    nf = 100; %number of frames
    fps = 1.25; %fps must be evenly divisible by 30!
    %llafocus = [40.4170, -3.7030, 0]; %Madrid Puerta Del Sol
    %llafocus = [40.4214, -3.7531, 0]; %Madrid Casa De Campo
    %llafocus = [40.4229, -3.7551, 0]; %Madrid Casa De Campo with Villa
    %llafocus = [40.418825 -3.697203 690]; %Gran Via (lla) focus point. altitude in meters
    %llafocus = [40.821059 14.426130 1200]; %Vesuvius
    %llafocus = [40.749300 14.484731 10]; %Pompei
    llafocus = [38.889, -77.035 170]; %MANUALLY INPUT DERIVE FROM GOOGLE EARTH
    scenarioName = 'DC';
    
    %MAKE LLAT and LLAFOCUS -----------------------------------------------    
    cam.DEM = getDEM(cam,llafocus,'geoid');
    nedfocus = lla2ned(cam.DEM,llafocus); %km
    az = linspace(-pi,pi,nf)'/4 - pi/2;
    nedcam = [r*cos(az) r*sin(az) repmat(-h+nedfocus(3),nf,1)];
    lla = [ned2lla(cam.DEM,nedcam), (1:1:nf)'./fps]; %llat
    llafocus = repmat(llafocus,[nf 1]);
end

%INITIALIZE VARIABLES -----------------------------------------------------
lla(:,1:3) = lla2llag(lla(:,1:3));  cam.frames=nf;  zva = zeros(nf,3);
cam.apriori.lla = zva;   cam.aposteriori.lla = zva;   cam.true.lla = lla(1:nf,1:3);
cam.apriori.ecef = zva;  cam.aposteriori.ecef = zva;  cam.true.ecef = lla2ecef(cam.true.lla);
cam.apriori.ned = zva;   cam.aposteriori.ned = zva;   cam.true.ned =  ecef2ned(cam.DEM,cam.true.ecef);
cam.apriori.rpy = zva;   cam.aposteriori.rpy = zva;   cam.true.rpy = zva;
cam.true.t = lla(:,4)'-lla(1,4);   cam.true.dt = [0 diff(cam.true.t)];   cam.fps = 1/mean(cam.true.dt(2:end));
cam.syntheticVideoFlag = true;

%FOCUS POINT --------------------------------------------------------------
llafocus(:,3) = cam.DEM.F(llafocus(:,2),llafocus(:,1));
cam.true.focus.lla = llafocus;  %sol, madrid
cam.true.focus.ecef = lla2ecef(cam.true.focus.lla);
cam.true.focus.ned = ecef2ned(cam.DEM,cam.true.focus.ecef);

%CAMERA TRAJECTORY --------------------------------------------------------
cam.width=cam.ssge(3);  cam.height=cam.ssge(4);
cam.fovh = 8; %deg (15
cam.focalLength = (cam.width/2) / tan(cam.fovh/2*d2r);   %pixels (0.00942 = pixel pitch in mm/pix)
cam.fovv = atan(cam.height/2/cam.focalLength)*r2d*2; %deg
cam.xppo = 0;  cam.yppo = 0;
cam.true.focus.pixel = camsc2pixel(cam, [1 0 0]);

cam.gpsxy1s = 5; %5m true
cam.gpsz1s = 10; %5m true
cam.rpy1s = 2; %deg
cam.roll1s = 0; %deg, MUST BE 0deg true roll for all GE Synthetic videos!!
cam.focuspointwander = 20; %m 1 sigma

[cam, ekf] = fcnaddSyntheticStochasticError(cam.DEM,cam);

%INCREMENT FRAMES ---------------------------------------------------------
he = actxserver('googleearth.ApplicationGE'); %Create a COM server running Google Earth

cam.filename = sprintf('GE Synthetic - %s %.0fx%.0f %.0fFOV',scenarioName,cam.width,cam.height,cam.fovh);
vid=VideoWriter(fcnincrementfname(cam.filename),'MPEG-4');  vid.FrameRate=30;  vid.Quality=100;  open(vid);  import java.awt.Robot;  mouse=Robot;
for i=1:cam.frames
    r=cam.true.focus.sc(i,1)*1000;
    el=cam.true.focus.sc(i,2)+90;
    az=cam.true.focus.sc(i,3);
    
    %UPDATE GE
    he.SetCameraParams(cam.true.focus.lla(i,1),cam.true.focus.lla(i,2),0,1,r,el,az,5);  %h.SetCameraParams(lat,lng,alt,altMode,range(m),tilt(deg),heading(deg),speed);
    
    %GET TIE POINTS & TOP LEFT PLOT
    pause(.05);  fcnFinishStreamProgress(he); 
    
    %SCREEN CAPTURE
    if mod(i,100)==0;  mouse.mouseMove(round(rand*100)+cam.ssge(3)+100, round(rand*100)+500); end %move mouse!
    writeVideo(vid, getscreen(cam.ssge))
    fprintf('frame %.0f/%.0f\n',i,cam.frames);
end
close(vid)

save(cam.filename,'cam')
end

