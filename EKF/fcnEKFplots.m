% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function []=fcnEKFplots(cam,ekf,a,index)
km2m = 1000;
nf = cam.frames;
x = 1:nf;
t = cam.true.t;
sg = erfinv(.683)*sqrt(2); %1sigma gain to 90% confidence
[h, hf] = fig(2,3);
ekf.x = zeros(size(ekf.xhat));

%POS
ha=h(1);  i=index.ned;  true=ekf.x(i,x)'*km2m;  aposteriori=ekf.xhat(i,x)'*km2m;  
e=true-aposteriori;   s=ekf.S(i,x)*km2m*sg;   fcnploterrors3(ha,t,x,e,s)
title(ha,'{\itaposteriori} {\bfposition} error with 90% confidence bounds'); ylabel(ha,'error (m)'); xlabel(ha,sprintf('time (s), %.1ffps',cam.fps)); 
legend(ha,'x','y','z'); %set(ha,'ylim',[-1 1]*12)

%VEL
ha=h(4);  i=index.dned;  true=ekf.x(i,x)'*km2m;  aposteriori=ekf.xhat(i,x)'*km2m; e=true-aposteriori;  s=ekf.S(i,x)*km2m*sg;  fcnploterrors3(ha,t,x,e,s)
title(ha,'{\itaposteriori} {\bfvelocity} error with 90% confidence bounds'); ylabel(ha,'error (m/s)'); xlabel(ha,sprintf('time (s), %.1ffps',cam.fps))
legend(ha,'x','y','z'); %set(ha,'ylim',[-1 1]*1.2)

%RPY
ha=h(2);  i=index.rpy;  true=ekf.x(i,x)';  aposteriori=ekf.xhat(i,x)';  e=true-aposteriori;  s=ekf.S(i,x)*sg;  fcnploterrors3(ha,t,x,e,s)
title(ha,'{\itaposteriori} {\bfrpy} error with 90% confidence bounds'); ylabel(ha,'error (deg)'); xlabel(ha,sprintf('time (s), %.1ffps',cam.fps))
legend(ha,'roll','pitch','yaw'); %set(ha,'ylim',[-1 1]*5)

%DRPY
ha=h(5);  i=index.drpy;  true=ekf.x(i,x)';  aposteriori=ekf.xhat(i,x)'; e=true-aposteriori;  s=ekf.S(i,x)*sg;  fcnploterrors3(ha,t,x,e,s)
title(ha,'{\itaposteriori} {\bfrpy rate} error with 90% confidence bounds'); ylabel(ha,'error (deg/s)'); xlabel(ha,sprintf('time (s), %.1ffps',cam.fps))
legend(ha,'roll','pitch','yaw'); %set(ha,'ylim',[-1 1]*30)

%K
ha=h(6);  i=index.K;  true=ekf.x(i,x)';  aposteriori=ekf.xhat(i,x)'; e=true-aposteriori;  s=ekf.S(i,x)*sg;  fcnploterrors3(ha,t,x,e,s)
title(ha,'{\itaposteriori} {\bfK} with 90% confidence bounds'); ylabel(ha,'value (pixels)'); xlabel(ha,sprintf('time (s), %.1ffps',cam.fps))
legend(ha,'xppo','yppo','focalLength');

%TP
ha=h(3);  s=zeros(3,nf);
for f=1:nf
    j = a.state(:,f)==1;  nj = sum(j);  xi = fcnai2xi(a,j);
    s(:,f) = mean( reshape(ekf.S(xi,f),[nj 3])); %sigmas
end
tx = ekf.xhat(index.tp0(:,1),:)*km2m;  tx(tx==0) = nan;
ty = ekf.xhat(index.tp0(:,2),:)*km2m;  ty(ty==0) = nan;
tz = ekf.xhat(index.tp0(:,3),:)*km2m;  tz(tz==0) = nan;
plot(ha,t,tx,'r',t,ty,'g',t,tz,'b');  axis(ha,'tight'); hold(ha,'on');  e=[]; fcnploterrors3(ha,t,x,e,s*km2m*sg)
title(ha,'{\itaposteriori} {\bfTP position} with 90% confidence bounds'); ylabel(ha,'NED Position (m)'); xlabel(ha,sprintf('time (s), %.1ffps',cam.fps))
%legend(ha,'x','y','z');

drawnow
end

function fcnploterrors3(ha,t,x,e,s)
if ~isempty(e);  plot(ha,t,e(x,1),'r',t,e(x,2),'g',t,e(x,3),'b');  hold(ha,'on'); end
plot(ha,t,s(1,x),'r',t,-s(1,x),'r','linewidth',2);  hold(ha,'on');
plot(ha,t,s(2,x),'g',t,-s(2,x),'g','linewidth',2)
plot(ha,t,s(3,x),'b',t,-s(3,x),'b','linewidth',2)

axis(ha,'tight')
box(ha,'off')
end