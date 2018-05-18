function [] = fcnPlotCorrelations()
tic
close all
load ECEFxyz40hr.mat
X = ecef2ned(x);

%i = 1:100;  x=x(i,:);  X=X(i,:);
vt = 1:1:(140000); 
nt = numel(vt);
nx = size(X,1);

h=fig(2,2,1.5);;

%ECEF
% c = zeros(nt,3);
% for i = 1:nt
%     j = vt(i);
%     a = zeros(nx-j,6);
%     a(:,1:3) = x(1:nx-j,:);
%     a(:,4:6) = x(1+j:nx,:);
%     C = corrcoef(a);
%     c(i,:) = [C(4,1) C(5,2) C(6,3)];
% end
c = zeros(nx,3);
c(:,1) = autocorr(x(:,1),nx-1);
c(:,2) = autocorr(x(:,2),nx-1);
c(:,3) = autocorr(x(:,3),nx-1);
c = c(vt,:);

axes(h(2))
plot(vt,c(:,1),'r-'); hold on; plot(vt,c(:,2),'g-'); plot(vt,c(:,3),'b-');
plot(vt,vt*0,'-','color',[1 1 1]*.8)
legend('X','Y','Z')
ylabel('Correlation Coefficient')
xlabel('\Deltat (s)')
title('uBlox LEA-5 {\bfECEF} Position Correlations'); axis tight

%NED
% c = zeros(nt,3);
% for i = 1:nt
%     j = vt(i);
%     a = zeros(nx-j,6);
%     a(:,1:3) = X(1:nx-j,:);
%     a(:,4:6) = X(1+j:nx,:);
%     C = corrcoef(a);
%     c(i,:) = [C(4,1) C(5,2) C(6,3)];
% end
c = zeros(nx,3);
c(:,1) = autocorr(X(:,1),nx-1);
c(:,2) = autocorr(X(:,2),nx-1);
c(:,3) = autocorr(X(:,3),nx-1);
c = c(vt,:);

axes(h(4))
plot(vt,c(:,1),'r-'); hold on; plot(vt,c(:,2),'g-'); plot(vt,c(:,3),'b-');
plot(vt,vt*0,'-','color',[1 1 1]*.8)
legend('X','Y','Z')
ylabel('Correlation Coefficient')
xlabel('\Deltat (s)')
title('uBlox LEA-5 {\bfNED} Position Correlations'); axis tight

axes(h(1))
plot3(X(:,1),X(:,2),X(:,3),'-','color',[.7 .7 .7],'markersize',6); 
set(gca,'zdir','reverse','ydir','reverse');
box on;  axis tight vis3d;
axis tight; title(sprintf('NED Errors\n 1\\sigma = [%.1f  %.1f  %.1f]m NED',std(X)))
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');

axes(h(3))
plot(1:nx,X(:,1),'r','linewidth',.5); hold on
plot(1:nx,X(:,2),'g','linewidth',.5); 
plot(1:nx,X(:,3),'b','linewidth',.5);
plot(1:nx,(1:nx)*0,'-','color',[1 1 1]*.8)
legend('X','Y','Z')
ylabel('error (m)')
xlabel('t (s)')
title('uBlox LEA-5 {\bfNED} Position Errors'); axis tight

set(findobj(gcf,'type','line'),'linesmoothing','on')
fcnfontsize
toc

figure
autocorr(X(:,3),nx-1);           % Inspect the ACF with 95% confidence


% % % clear all
% % % f1=8;
% % % f2=64;
% % % NFFT=1024;
% % % dt=1/max([f1 f2])/8;
% % % Freq=1/dt;
% % % df=Freq/NFFT;
% % % f=(0:NFFT/2)*df;
% % % t=(0:NFFT-1)*dt;
% % % y=sin(2*pi*f1*t)+2*sin(2*pi*f2*t)+1.1;
% 
% h = fig(2,2,1.5);;
% NFFT = nx;
% t = (1:nx)'; f=t;
% y = x(:,1);
% 
% fy1=fft(y)*2/NFFT;
% fy2=fft(y-mean(y))*2/NFFT;
% 
% axes(h(1))
% plot(t,y-mean(y))
% xlabel('time(s)')
% ylabel('response')
% title ('time history w/ mean removed')
% 
% axes(h(2))
% plot(f(1:NFFT/4),abs(fy1(1:NFFT/4)))
% xlabel('Freq(Hz)')
% ylabel('response')
% title ('Frequency Response');grid
% 
% axes(h(3))
% plot(f(1:NFFT/4),abs(fy2(1:NFFT/4)))
% xlabel('Freq(Hz)')
% ylabel('response')
% title ('Frequency Response w/DC removed')
% 
% loglog(f(1:NFFT/4),abs(fy2(1:NFFT/4)))

end