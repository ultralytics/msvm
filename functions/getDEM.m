function DEM = getDEM(cam,lla,geoidstr)
if nargin==2; geoidstr = 'ellipsoid'; end %elevations off ellipsoid instead of geoid

ni=20; %DEM points along x and y
mlla = mean(lla,1);  r = max(mlla(3)/1000*1.5, 1); %km extension about lla edges
kmperdeg = fcnmperLLAdeg(mlla)/1000;  dext = r./kmperdeg; % degrees of extension
latv = linspace(min(lla(:,1)) - dext(1), max(lla(:,1)) + dext(1), ni);
lngv = linspace(min(lla(:,2)) - dext(2), max(lla(:,2)) + dext(2), ni);

if isfield(cam,'DEM') && isfield(cam.DEM,'Fned')
    fprintf('Discovered Existing DEM in getDEM()\n')
    DEM = cam.DEM;
else %get DEM
    DEM = fcnGoogleElevationAPIDEM(latv,lngv,geoidstr);
end