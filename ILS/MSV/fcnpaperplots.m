function [] = fcnpaperplots(cam)

try
    vfr=VideoReader([cam.filename '.mp4']);
catch
    vfr=VideoReader([cam.filename '.avi']);
end

centerpixel = [(cam.width+1)/2 (cam.height+1)/2];
sc = pixel2camsc(cam,centerpixel); 
rpy_focus2centerpixel = [0 sc(2:3)*d2r];
C_focus2centerpixel = fcnRPY2DCM_W2B(rpy_focus2centerpixel);

try
    r = mean(fcnrange(mean(cam.tpnedhat,1),cam.apriori.ned))*.4; %focal length
catch
    r = mean(fcnrange(cam.apriori.ned(1,:),cam.apriori.ned(end,:)))/3;
end
ov = ones(1,10);
vv = cam.fovv/2; %half fov vertical
vh = cam.fovh/2; %half fov horizontal
w = cam.width/20;
h = cam.height/20;

try 
    frames = cam.msv.vf0;
catch
    frames = linspace(1,cam.frames,3);
    cam.frameID = 1:cam.frames;
end

for i = frames
    Ci = C_focus2centerpixel*fcnRPY2DCM_W2B(cam.aposteriori.rpy(i,:)*d2r);
    
    el = [ov*-vv, linspace(-vv, vv, 10), ov*vv, linspace(vv, -vv, 10)]*d2r;
    az = [linspace(-vh, vh, 10), ov*vh, linspace(vh, -vh, 10), ov*-vh]*d2r;
    x = [fcnSC2CC([ones(40,1)*r el' az']); [0 0 0]; fcnSC2CC([r vv*d2r -vh*d2r]); [0 0 0]; fcnSC2CC([r vv*d2r vh*d2r]); [0 0 0]; fcnSC2CC([r -vv*d2r vh*d2r])];
    x = bsxfun(@plus,cam.aposteriori.ned(i,:),x*Ci);
    plot3(x(:,1),x(:,2),x(:,3),'-','color',[0 0 0]);
    
    [X,Y] = ndgrid(linspace(-vh, vh, w), linspace(vv, -vv, h));  x = fcnSC2CCd([Y(:)*0+r, Y(:), X(:)]);
    x = bsxfun(@plus,cam.aposteriori.ned(i,:),x*Ci);
    X=reshape(x(:,1),[w h]);  Y=reshape(x(:,2),[w h]);  Z=reshape(x(:,3),[w h]);
    C = read(vfr,cam.frameID(i));
    surf(X',Y',Z',C,'edgecolor','none','facecolor','texture');
end

fcnplot3(cam.true.ned,'b.','markersize',10); 
%trimesh(cam.DEM.delaunayfaces,cam.DEM.ned(:,1),cam.DEM.ned(:,2),cam.DEM.ned(:,3),cam.DEM.nedz(:),'edgecolor',[.7 .7 .7]); hold on


fcntight('xyz')