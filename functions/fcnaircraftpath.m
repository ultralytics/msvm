% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [cam, ekf] = fcnaircraftpath(DEM,cam)
% %pick points
% ni = 1000;
% [h1, hf]=fig();
% [xmin,xmax]=fcnminmax(DEM.lng);  [ymin,ymax]=fcnminmax(DEM.lat); 
% button=1; i=0; x=[]; y=[]; h2=[];   
% plot3(h1,cam.true.focus.lla(2),cam.true.focus.lla(1),10,'.w','Markersize',20); hold on;  image([xmin xmax],[ymin ymax],flipdim(imrotate(DEM.cdata,00),1),'parent',h1); 
% alpha(.5); view(0,90); axis(h1,[xmin xmax ymin ymax]); daspect(h1,[1 1 1]); drawnow
% while button==1
%     i=i+1;
%     %cla(h1);
%     [x(i),y(i),button] = ginput(1);
%     plot(h1,x,y,'b.-');
%     if i==1;  continue;  end
%     xi = interp1(linspace(0,1,i),x,linspace(0,1,ni),'spline')';
%     yi = interp1(linspace(0,1,i),y,linspace(0,1,ni),'spline')';
%     delete(h2);  h2=plot(h1,xi,yi,'m-');
%     axis(h1,[xmin xmax ymin ymax]); daspect(h1,[1 1 1]);  view(0,90);
% end
% close(hf)
% lla = [yi, xi, linspace(4000,3000,ni)'];  lla(:,3)=lla(:,3)+DEM.F(lla(:,2),lla(:,1));  save airtrack.mat lla
% ned = ecef2ned(DEM,lla2ecef(lla));


%load airtrackline.mat
load airtrack2circles.mat
lla = lla(180:end-500,:); %#ok<NODEF>
ni1 = size(lla,1);

lla(:,1:2) = lla(:,1:2)+ones(ni1,1)*(cam.true.LLAstart(1:2) - mean(lla(:,1:2))); %center track on cam.true.LLAstart

lla = [lla(:,1), lla(:,2), linspace(2500,2000,ni1)'];  lla(:,3)=lla(:,3)+DEM.F(lla(:,2),lla(:,1));
vel = .05; %km/s
llasmoothed = interp1(linspace(0,1,ni1),lla,linspace(0,1,1E4),'spline');
cam.traveldistance = sum(fcnrange(diff(lla2ecef(llasmoothed)))); %km
cam.traveltime = cam.traveldistance/vel; %s

%NEW WAY
cam.true.dt = 1/cam.fps; %time between frames
cam.frames = floor(cam.traveltime/cam.true.dt);
[ekf, cam] = fcnbuilderrors(cam);

%FOCUS POINT
LLA = cam.true.LLAstart; %central focus point
mperdeg = fcnmperLLAdeg(LLA);
fpw = bsxfun(@rdivide, cam.true.focus.stochasticwander(:,1:2)*1E3, mperdeg);
cam.true.focus.lla=[LLA(1)+fpw(:,1), LLA(2)+fpw(:,2), zeros(cam.frames,1)];  cam.true.focus.lla(:,3)=cam.true.focus.lla(:,3)+DEM.F(cam.true.focus.lla(:,2),cam.true.focus.lla(:,1));  %sol, madrid
cam.true.focus.ecef = lla2ecef(cam.true.focus.lla);
cam.true.focus.ned = ecef2ned(DEM,cam.true.focus.ecef);

cam.true.lla = interp1(linspace(0,1,ni1),lla,linspace(0,1,cam.frames),'spline'); %upsample or downsample
cam.true.ecef = lla2ecef(cam.true.lla);
cam.true.ned = ecef2ned(DEM,cam.true.ecef);
cam.apriori.ned = cam.true.ned - ekf.x([1 3 5],:)';
cam.apriori.ecef = ned2ecef(DEM,cam.apriori.ned);

[~, dx] = fcnrange(cam.true.ned,cam.true.focus.ned);
sc = fcnCC2SC(dx);  %sc(:,2)=sc(:,2)+90;
cam.true.focus.sc = sc;
cam.true.focus.gsd = sc(:,1).*tand(cam.fovh/cam.width)*1000;

cam.true.rpy = [randn(cam.frames,1)*cam.roll1s sc(:,2:3)];
cam.apriori.rpy = cam.true.rpy - ekf.x([7 9 11],:)';

cam.aposteriori.ned  = cam.apriori.ned;
cam.aposteriori.ecef = cam.apriori.ecef;
cam.aposteriori.rpy  = cam.apriori.rpy;
end


function [ekf, cam] = fcnbuilderrors(cam)
dt = cam.true.dt;
ns = 12; %number of states
ni = cam.frames;
r = randn(ns,ni);
m2km = 1/1000;
cam.true.t = 0:dt:(dt*cam.frames-dt);

%position and velocity ----------------------------------------------------
B1a = 1/108; %seconds
B2a = 1/162; %seconds
[s2a, P0a] = KF2_simple_baseline_april_10_2012_init_cov(1/B1a, 1/B2a, cam.xyz1s*m2km); %meters
%s2a = .18*m2km; %km/s^2
s1a = s2a/2; %km/s
[Phi61, Q61] = fcngetPhiQ(s1a,s2a,B1a,B2a,dt);

%rotational position and velocity -----------------------------------------
B1b = 1/36; %seconds
B2b = 1/54; %seconds
[s2b, P0b] = KF2_simple_baseline_april_10_2012_init_cov(1/B1b, 1/B2b, cam.rpy1s); %deg
%s2b = 1/100; %deg/s^2
s1b = s2b/2; %deg/s
[Phi62, Q62] = fcngetPhiQ(s1b,s2b,B1b,B2b,dt);

%focus point position wander ----------------------------------------------
B1c = 1/8; %seconds
B2c = 1/12; %seconds
[s2c, P0c] = KF2_simple_baseline_april_10_2012_init_cov(1/B1c, 1/B2c, cam.focuspointwander*m2km); %meters
s1c = 0*s2c/2; %deg/s
[Phi63, Q63] = fcngetPhiQ(s1c,s2c,B1c,B2c,dt);
wfpw = sqrtm(Q63)*randn(6,ni);
Xfpw = zeros(6,ni);

%combine ------------------------------------------------------------------
zm2=zeros(2);
zm6=zeros(6);
Phi = [ Phi61     zm6
        zm6       Phi62  ];
Q =   [ Q61       zm6
        zm6       Q62    ];
    
Qsr = sqrtm(Q);
w = Qsr*r;
%initialize X and P -------------------------------------------------------
%P = zeros(ns,ns,ni);
P(:,:,1) = [ P0a     zm2     zm2     zm2     zm2     zm2
             zm2     P0a     zm2     zm2     zm2     zm2
             zm2     zm2     P0a     zm2     zm2     zm2
             zm2     zm2     zm2     P0b     zm2     zm2
             zm2     zm2     zm2     zm2     P0b     zm2
             zm2     zm2     zm2     zm2     zm2     P0b];

X = zeros(ns,ni);    
X(1:12,1) = sqrtm(P(:,:,1))*randn(ns,1) + [0 0 0 0 0 0 0 0 1 0 0 0]'*0; %initial errors
    
%propagate X and P --------------------------------------------------------
for i=2:ni
    X(:,i) = Phi*X(:,i-1) + w(:,i);
    %P(:,:,i) = Phi*P(:,:,i-1)*Phi' + Q;
    Xfpw(:,i) = Phi63*Xfpw(:,i-1) + wfpw(:,i);
end

ekf.x=X;
ekf.Phi12=Phi;
ekf.Q12=Q;
ekf.P0_12 = P(:,:,1);
cam.true.focus.stochasticwander = Xfpw([1 3 5],:)';

% h=fig(2,2);
% t = cam.true.t;
% P11 = reshape(P(1,1,ti),[frames 1])';
% P22 = reshape(P(2,2,ti),[frames 1])';
% P33 = reshape(P(7,7,ti),[frames 1])';
% P44 = reshape(P(8,8,ti),[frames 1])';
% plot(h(1),t,X(1,:),'-r',t,X(3,:),'-g',t,X(5,:),'-b',t,sqrt(P11),'-k',t,-sqrt(P11),'-k','linewidth',2); title(h(1),'xyz errors')
% plot(h(2),t,X(2,:),'-r',t,X(4,:),'-g',t,X(6,:),'-b',t,sqrt(P22),'-k',t,-sqrt(P22),'-k','linewidth',2); title(h(2),'dxyz errors')
% plot(h(3),t,X(7,:),'-r',t,X(9,:),'-g',t,X(11,:),'-b',t,sqrt(P33),'-k',t,-sqrt(P33),'-k','linewidth',2); title(h(3),'rpy errors')
% plot(h(4),t,X(8,:),'-r',t,X(10,:),'-g',t,X(12,:),'-b',t,sqrt(P44),'-k',t,-sqrt(P44),'-k','linewidth',2); title(h(4),'drpy errors')
end

function [Phi61, Q61] = fcngetPhiQ(s1,s2,B1,B2,dt)
zm = zeros(2);
Phi2 = [ exp(-B1*dt)                      (B2-B1)^-1*(exp(-B1*dt)-exp(-B2*dt))
        0                                 exp(-B2*dt)                                 ];
Phi61 = [ Phi2    zm      zm
          zm      Phi2    zm
          zm      zm      Phi2 ];

Q2 = zm;
Q2(1,1) = s2^2*(B2-B1)^-2*(.5*B1^-1*(1-exp(-2*B1*dt))+.5*B2^-1*(1-exp(-2*B2*dt))-2*(B1+B2)^-1*(1-exp(-(B1+B2)*dt)))  +  s1^2*.5*B1^-1*(1-exp(-2*B1*dt));
Q2(1,2) = s2^2*(B2-B1)^-1*((B1+B2)^-1*(1-exp(-(B1+B2)*dt)) - .5*B2^-1*(1-exp(-2*B2*dt)));
Q2(2,1) = Q2(1,2);
Q2(2,2) = s2^2*.5*B2^-1*(1-exp(-2*B2*dt));
Q61 = [ Q2    zm    zm
        zm    Q2    zm 
        zm    zm    Q2 ];  
end