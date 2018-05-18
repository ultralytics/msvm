clc
close all
clear all

addpath([cd '\GEfiles'])

%SCREEN SIZE TO RECORD
ss = get(0,'ScreenSize');
ssge = [255 40 1680-255 1050-40-75];
ssge = [255+400 40+200 900 600];

%INITIAL CONDITIONS
alt     = 1000; %m
altMode = 1;
range   = 40; %m
tilt    = 75; %off nadir, 90 = horizon
heading = -160;
speed   = .5;
lng0     = 69.179;
lat0     = 34.523;

h       = actxserver ('googleearth.ApplicationGE');%  Create a COM server running Google Earth
aviobj  = avifile('large.avi');
aviobj.quality = 10;
aviobj.compression = 'none';
aviobj.fps = 30;

%CIRCULAR WAYPOINTS
np = 30;
r = .05; %lat
az = linspace(0,360,np);
dlng = r*sind(az);
dlat = r*cosd(az);
heading = linspace(-30,30,np);
tilt = linspace(90,60,np);
range = linspace(20,50,np);


for i = 1:np
    fprintf('Frame %.0f\n',i)
    lat = lat0 + dlat(i);
    lng = lng0 + dlng(i);
    
    %UPDATE MQ9 KML
    rpy = [-30 0 mod(az(i)+90,360)]; %yaw only valid from -360 to 360
    MQ9kmlname = [cd '\GEfiles\MQ9.kml'];

    %h.SetCameraParams(lat-1,lng,alt,altMode,range(i),tilt(i),heading(i),speed);
    fcnwriteMQ9([lat lng alt], rpy, heading(i), tilt(i), range(i)); pause(2); winopen(MQ9kmlname); pause(2);
   
    %UPDATE GE SCREEN
    %h.SetCameraParams(lat-1,lng,alt,altMode,range(i),tilt(i),heading(i),speed);
    fcnFinishStreamProgress(h)
    pause(.001)
    
    %GET SCREENSHOT!
    F = getscreen(ssge); %image(x.cdata); axis equal
    aviobj = addframe(aviobj,F);
end

%CLOSE STUFF!
aviobj = close(aviobj);
delete(h)




