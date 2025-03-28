% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [h] = fcnplotcam(ha,xhat,C,fovh,fovv,c,R,t)
r2d = 180/pi;
d2r = pi/180;
fv=fovv/2;
fh=fovh/2;

r = .5;

no = 4; %number of ones
ov = ones(1,no);

el = [ov*-fv, linspace(-fv, fv, no), ov*fv, linspace(fv, -fv, no)]*d2r;
az = [linspace(-fh, fh, no), ov*fh, linspace(fh, -fh, no), ov*-fh]*d2r;

x = [fcnSC2CC([ones(no*4,1)*r el' az']); [0 0 0]; fcnSC2CC([r fv*d2r -fh*d2r]); [0 0 0]; fcnSC2CC([r fv*d2r fh*d2r]); [0 0 0]; fcnSC2CC([r -fv*d2r fh*d2r])];
x=fcnxform(x,xhat);
%h=plot3(ha,x(:,1),x(:,2),x(:,3),'-','color',[.85 .85 1]); hold(ha,'on')
h=plot3(ha,t(1),t(2),t(3),'.','color',c); hold(ha,'on')

C = fcnRPY2DCM_W2B(xhat(4:6));

%R = fcnVEC2DCM_W2B(-t)';
%R=R';

% x = [0 0 0; .3 0 0]';
% x = R*x + [t; t]';
% h=plot3(ha,x(1,:),x(2,:),x(3,:),'-','color','r');
% 
% x = [0 0 0; 0 .3 0]';
% x = R*x + [t; t]';
% h=plot3(ha,x(1,:),x(2,:),x(3,:),'-','color','g');

x = [0 0 0; 0 0 5]';
x = R*x + [t; t]';
h=plot3(ha,x(1,:),x(2,:),x(3,:),'-','color','b');


% %plot image
% [X,Y] = ndgrid(linspace(-fh, fh, 640), linspace(fv, -fv, 480));
% 
% x = fcnSC2CCd([ones(640*480,1)*r, reshape(Y,[640*480 1]), reshape(X,[640*480 1])]);
% x=fcnxform(x,xhat);
% X=reshape(x(:,1),[640 480]);  Y=reshape(x(:,2),[640 480]);  Z=reshape(x(:,3),[640 480]);
% 
% %h = [h, surf(ha,X',Y',Z',C,'edgecolor','none')];
% %h = [h, surf(ha,X',Y',Z',C,'edgecolor','none','facecolor','texture','cdata',C)];

end


