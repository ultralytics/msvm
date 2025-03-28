% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [hf] = fcncompareILS(cam,ekf,ils,a)
m2km = 1/1000;
km2m = 1000;
dt = cam.true.dt;
nframes = size(ils.xhat12,1);

%EKF ----------------------------------------------------------------------
m2km = 1/1000;
km2m = 1000;
dt = cam.true.dt;

%prepare
x = 1:nframes;
t = x*dt;
sigmagain = erfinv(.90)*sqrt(2); %1sigma to 90% confidence
mpd = sum3(abs(diff(a.upx,1,2)))/numel(a.upx); %mean pixel difference between frames
%[h, hf] = fig(2,3);
% 
% %plot aposteriori pos error with 1sigmas
% ha=h(1);
% v1=[1 3 5];  true=ekf.x(v1,x)'*km2m;  aposteriori=ekf.xhat(v1,x)'*km2m;  
% y1=true-aposteriori;   y2=ekf.xhatsigma(v1,x)'*km2m*sigmagain;   fcnploterrors3(ha,t,x,y1,y2)
% title(ha,'EKF {\itaposteriori} {\bfposition} error with 90% confidence bounds'); ylabel(ha,'error (m)'); xlabel(ha,sprintf('time (s), %.1ffps, %.2f pixels mean TP shift between frames',cam.fps,mpd)); 
% legend(ha,'x','y','z'); set(ha,'ylim',[-1 1]*12); alims{1} = axis;
% 
% %plot aposteriori vel error with 1sigmas
% ha=h(4); 
% v1=[2 4 6];  true=ekf.x(v1,x)'*km2m;  aposteriori=ekf.xhat(v1,x)'*km2m; 
% y1=true-aposteriori;  y2=ekf.xhatsigma(v1,x)'*km2m*sigmagain;  fcnploterrors3(ha,t,x,y1,y2)
% title(ha,'EKF {\itaposteriori} {\bfvelocity} error with 90% confidence bounds'); ylabel(ha,'error (m/s)'); xlabel(ha,sprintf('time (s), %.1ffps',cam.fps))
% legend(ha,'x','y','z'); set(ha,'ylim',[-1 1]*1.2); alims{2} = axis;
% 
% %plot aposteriori rpy error with 1sigmas
% ha=h(2); 
% v1=[7 9 11];  true=ekf.x(v1,x)';  aposteriori=ekf.xhat(v1,x)';
% y1=true-aposteriori;  y2=ekf.xhatsigma(v1,x)'*sigmagain;  fcnploterrors3(ha,t,x,y1,y2)
% title(ha,'EKF {\itaposteriori} {\bfrpy} error with 90% confidence bounds'); ylabel(ha,'error (deg)'); xlabel(ha,sprintf('time (s), %.1ffps',cam.fps))
% legend(ha,'roll','pitch','yaw'); set(ha,'ylim',[-1 1]*.5); alims{3} = axis;
% 
% ha=h(5);
% v1=[8 10 12];  true=ekf.x(v1,x)';  aposteriori=ekf.xhat(v1,x)'; 
% y1=true-aposteriori;  y2=ekf.xhatsigma(v1,x)'*sigmagain;  fcnploterrors3(ha,t,x,y1,y2)
% title(ha,'EKF {\itaposteriori} {\bfrpy rate} error with 90% confidence bounds'); ylabel(ha,'error (deg/s)'); xlabel(ha,sprintf('time (s), %.1ffps',cam.fps))
% legend(ha,'roll','pitch','yaw'); set(ha,'ylim',[-1 1]*.05); alims{4} = axis;
% 
% %tp1
% ha=h(3);
% v1=[13 14 15];  true=ones(numel(x),1)*a.ipned(1,:)*km2m;  aposteriori=ekf.xhat(v1,x)'*km2m; 
% y1=true-aposteriori;  y2=ekf.xhatsigma(v1,x)'*sigmagain*km2m;  fcnploterrors3(ha,t,x,y1,y2)
% title(ha,'EKF {\itaposteriori} TP1 {\bfposition} error with 90% confidence bounds'); ylabel(ha,'error (m)'); xlabel(ha,sprintf('time (s), %.1ffps',cam.fps)); 
% legend(ha,'x','y','z'); set(ha,'ylim',[-1 1]*40); alims{5} = axis;
% 
% ha=h(6);  cla(ha);
% y=true;         plot(ha,t,y(x,1),'r',t,y(x,2),'g',t,y(x,3),'b','linesmoothing','off','linewidth',2); hold(ha,'on')
% y=aposteriori;  plot(ha,t,y(:,1),'color',[1 .8 .8],'linesmoothing','on'); plot(ha,t,y(:,2),'color',[.8 1 .8],'linesmoothing','on'); plot(ha,t,y(:,3),'color',[.8 .8 1],'linesmoothing','on');
% title(ha,'EKF true and {\itaposteriori} tp1 position'); ylabel(ha,'m'); xlabel(ha,sprintf('time (s), %.1ffps',cam.fps)); 
% axis(ha,'tight')
% box(ha,'off')
% drawnow
%
%EKF ----------------------------------------------------------------------

%prepare
x = 1:nframes;
t = x*dt;
sigmagain = erfinv(.90)*sqrt(2); %1sigma to 90% confidence
mpd = sum3(abs(diff(a.upx,1,2)))/numel(a.upx); %mean pixel difference between frames
[h, hf] = fig(2,3);

%plot aposteriori pos error with 1sigmas
ha=h(1);
v1=[1 3 5];  true=ekf.x(v1,x)'*km2m;  aposteriori=ils.xhat12(:,v1)*km2m;  
y1=true-aposteriori;   y2=ils.xhatsigma12(x,v1)*km2m*sigmagain;   fcnploterrors3(ha,t,x,y1,y2)
title(ha,'ILS {\itaposteriori} {\bfposition} error with 90% confidence bounds'); ylabel(ha,'error (m)'); xlabel(ha,sprintf('time (s), %.1ffps, %.2f pixels mean TP shift between frames',cam.fps,mpd)); 
legend(ha,'x','y','z'); set(ha,'ylim',[-1 1]*12);

%plot aposteriori vel error with 1sigmas
ha=h(4); 
v1=[2 4 6];  true=ekf.x(v1,x)'*km2m;  aposteriori=ils.xhat12(:,v1)*km2m;  
y1=true-aposteriori;  y2=ils.xhatsigma12(x,v1)*km2m*sigmagain;  fcnploterrors3(ha,t,x,y1,y2)
title(ha,'ILS {\itaposteriori} {\bfvelocity} error with 90% confidence bounds'); ylabel(ha,'error (m/s)'); xlabel(ha,sprintf('time (s), %.1ffps',cam.fps))
legend(ha,'x','y','z'); set(ha,'ylim',[-1 1]*1.2);

%plot aposteriori rpy error with 1sigmas
ha=h(2); 
v1=[7 9 11];  true=ekf.x(v1,x)';  aposteriori=ils.xhat12(:,v1);
y1=true-aposteriori;  y2=ils.xhatsigma12(x,v1)*sigmagain;  fcnploterrors3(ha,t,x,y1,y2)
title(ha,'ILS {\itaposteriori} {\bfrpy} error with 90% confidence bounds'); ylabel(ha,'error (deg)'); xlabel(ha,sprintf('time (s), %.1ffps',cam.fps))
legend(ha,'roll','pitch','yaw'); set(ha,'ylim',[-1 1]*.5);

ha=h(5);
v1=[8 10 12];  true=ekf.x(v1,x)';  aposteriori=ils.xhat12(:,v1); 
y1=true-aposteriori;  y2=ils.xhatsigma12(x,v1)*sigmagain;  fcnploterrors3(ha,t,x,y1,y2)
title(ha,'ILS {\itaposteriori} {\bfrpy rate} error with 90% confidence bounds'); ylabel(ha,'error (deg/s)'); xlabel(ha,sprintf('time (s), %.1ffps',cam.fps))
legend(ha,'roll','pitch','yaw'); set(ha,'ylim',[-1 1]*.05);

%tp1
ha=h(3);
v1=[1];  true=ones(numel(x),1)*a.ipned(1,:)*km2m;  aposteriori=ones(numel(x),1)*ils.xhattp(v1,:)*km2m; 
y1=true-aposteriori;  y2=ones(numel(x),1)*ils.xhattpsigma(v1,:)*sigmagain*km2m;  fcnploterrors3(ha,t,x,y1,y2)
title(ha,'ILS {\itaposteriori} TP1 {\bfposition} error with 90% confidence bounds'); ylabel(ha,'error (m)'); xlabel(ha,sprintf('time (s), %.1ffps',cam.fps)); 
legend(ha,'x','y','z'); set(ha,'ylim',[-1 1]*40);

% ha=h(6);  cla(ha);
% y=true;         plot(ha,t,y(x,1),'r',t,y(x,2),'g',t,y(x,3),'b','linesmoothing','off','linewidth',2); hold(ha,'on')
% y=aposteriori;  plot(ha,t,y(:,1),'color',[1 .8 .8],'linesmoothing','on'); plot(ha,t,y(:,2),'color',[.8 1 .8],'linesmoothing','on'); plot(ha,t,y(:,3),'color',[.8 .8 1],'linesmoothing','on');
% title(ha,'true and {\itaposteriori} tp1 position'); ylabel(ha,'m'); xlabel(ha,sprintf('time (s), %.1ffps',cam.fps)); 
% box(ha,'off')
% drawnow

fcntight(h,'x')


drawnow
end

function fcnploterrors3(ha,t,x,y1,y2)
cla(ha)

plot(ha,t,y1(x,1),'r',t,y1(x,2),'g',t,y1(x,3),'b','linesmoothing','on');  hold(ha,'on')
plot(ha,t,y2(x,1),'r',t,-y2(x,1),'r','linewidth',2,'linesmoothing','on')
plot(ha,t,y2(x,2),'g',t,-y2(x,2),'g','linewidth',2,'linesmoothing','on')
plot(ha,t,y2(x,3),'b',t,-y2(x,3),'b','linewidth',2,'linesmoothing','on')

axis(ha,'tight')
box(ha,'off')
end

