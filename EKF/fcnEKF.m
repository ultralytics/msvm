function [cam] = fcnEKF(cam,ekf,a)
%[x y z, dx dy dz, r p y, dr dp dy, xppo yppo focalLength,  tp1x tp1y tp1z, ...  tpnz   tpny   tpnz]', (15+3n x 1)
% 1 2 3  4   5  6  7 8 9  10 11 12    13   14          15     16   17   18, ... 10+3n  11+3n  12+3n]
startclock = clock;
ekf = fcnT2PPhiQ(cam,ekf,cam.true.dt(2));
[ntp, nf] = size(a.state);  nx = ntp*3+15;  ntp0 = sum(a.state(:,1)==1);
m2km = .001;

index.ned = [1 3 5];  index.dned = [2 4 6];
index.rpy = [7 9 11];  index.drpy = [8 10 12];  index.K = [13 14 15];
index.tp0 = reshape(1:ntp*3, [3 ntp])'+15;

xhat = [zeros(12,1); cam.xppo; cam.yppo; cam.focalLength; zeros(ntp0*3,1)];
P = eye(nx)*(500*m2km).^2;  P(index.K,index.K) = diag(cam.sigmas.K.^2);  P(1:12,1:12) = ekf.P0_12; % km, km/s, deg, deg/s

a.zi = zeros(ntp,nf);
X = zeros(nx,nf);
S = zeros(nx,nf);
fprintf('\nStarting EKF: %d tie points over %d images (%d measurements of %d states)...\n',ntp,nf,sum3(a.state==1),nx)
fprintf('image    measurements (tie points)    states (tie points)    residual RMS   residual mean\n')
for f=1:nf
    dt = cam.true.dt(f);
    as = a.state(:,f);
    new = as==1 & a.state(:,max(f-1,1))==0;
    
    index.za    = as==1 & ~new;
    index.xhata = find(as~=0 & ~new);   nxtp = numel(index.xhata);
    %index.xhata = find(as==1 & ~new);   nxtp = numel(index.xhata);
    m = reshape(1:nxtp*3,[3 nxtp])+15;   index.tp = m(:,as(index.xhata)==1)'; 
    index.X = [(1:15)'; fcnai2xi(a,index.xhata)];
    
    if f>1
        [xhat, P, index] = fcnnewTP(cam,X(index.X,f-1),P,a,index,f);
        
        z = [a.upx(index.za,f) a.upy(index.za,f)]';  z = z(:);  nz = numel(z);
        R = speye(nz)*(ekf.zsigma.^2);
        [cam, ekf, xhat, P(index.X,index.X), res] = fcnEKF1(cam, ekf, index, xhat, P(index.X,index.X), dt, z, R, f);
        fprintf('I%-4d%16d  %-4d%17d  %-4d%23.2f%16.2f\n',f,nz,nz/2,nxtp*3+15,nxtp,fcnrms(res),mean(res(:)))
        if fcnrms(res)>300; break; end
    end

    a.zi(:,f) = index.za;
    X(index.X,f) = xhat;
    S(index.X,f) = sqrt(diag(P(index.X,index.X)));
end
fprintf('Done (%.0fs).\n\n',etime(clock,startclock))

ekf.xhat = X(:,1:f);
ekf.S = S(:,1:f);
fcnEKFplots(cam,ekf,a,index)
ekf.index = index;
cam.ekf = ekf;

fcnplotmatrix(real(cam.ekf.S(cam.ekf.index.tp0(:,3),:))*1000); set(gca,'clim',[0 10]); axis normal
xlabel('image'); ylabel('Tie Point ID'); title('Tie Point Z 1\sigma'); axis xy; fcncolorbar(1,'m')
end


function [xhat, P, index] = fcnnewTP(cam,xhat,P,a,index,f) %discover new tie points and give them an initial guess on the ground
new = all(xhat(index.tp)==0,2); %x, y and z = 0
index.xhatnewtp = index.tp(new,:);

if any(new)
    ai = find(index.za);
    b = fcncropTP(cam,a,ai(new));
    j = index.tp(new,:);

    %APRIORI
    %xhat(j) = cam.tpnedhat(ai(new),:);
    
    %DEM INTERCEPT 
    %xhat(j)  = pixel2ned(cam.DEM,cam,f,[b.upx(:,f) b.upy(:,f)],'apriori');
    xhat(j)  = pixel2ned(cam.DEM,cam,f-1,[b.upx(:,f-1) b.upy(:,f-1)],'aposteriori');
    
    %MIGMSV
    %[cam, b] = fcncropimages(cam,b,[f-1 f]);
    %xhat(j) = fcnMIGMSV(cam,b);
    
    if f>10
        j = j(:);
        P(index.X(j),index.X(j)) = eye(numel(j))*(.5).^2; %km
    end
end
end


function []=fcnEKFplots(cam,ekf,a,index)
km2m = 1000;
nf = numel(ekf.xhat(1,:));
x = 1:nf;
t = x; %cam.true.t;
xstr = sprintf('image, %.1ffps',cam.fps);

sg = erfinv(.683)*sqrt(2); %1sigma gain to 90% confidence
[h, hf] = fig(2,3);
ekf.x = zeros(size(ekf.xhat));

%POS
ha=h(1);  i=index.ned;  true=ekf.x(i,x)'*km2m;  aposteriori=ekf.xhat(i,x)'*km2m;  
e=true-aposteriori;   s=ekf.S(i,x)*km2m*sg;   fcnploterrors3(ha,t,x,e,s)
title(ha,'{\itaposteriori} {\bfposition} error with 90% confidence bounds'); ylabel(ha,'error (m)'); xlabel(ha,xstr); 
legend(ha,'x','y','z'); %set(ha,'ylim',[-1 1]*12)

%VEL
ha=h(4);  i=index.dned;  true=ekf.x(i,x)'*km2m;  aposteriori=ekf.xhat(i,x)'*km2m; e=true-aposteriori;  s=ekf.S(i,x)*km2m*sg;  fcnploterrors3(ha,t,x,e,s)
title(ha,'{\itaposteriori} {\bfvelocity} error with 90% confidence bounds'); ylabel(ha,'error (m/s)'); xlabel(ha,xstr)
legend(ha,'x','y','z'); %set(ha,'ylim',[-1 1]*1.2)

%RPY
ha=h(2);  i=index.rpy;  true=ekf.x(i,x)';  aposteriori=ekf.xhat(i,x)';  e=true-aposteriori;  s=ekf.S(i,x)*sg;  fcnploterrors3(ha,t,x,e,s)
title(ha,'{\itaposteriori} {\bfrpy} error with 90% confidence bounds'); ylabel(ha,'error (deg)'); xlabel(ha,xstr)
legend(ha,'roll','pitch','yaw'); %set(ha,'ylim',[-1 1]*5)

%DRPY
ha=h(5);  i=index.drpy;  true=ekf.x(i,x)';  aposteriori=ekf.xhat(i,x)'; e=true-aposteriori;  s=ekf.S(i,x)*sg;  fcnploterrors3(ha,t,x,e,s)
title(ha,'{\itaposteriori} {\bfrpy rate} error with 90% confidence bounds'); ylabel(ha,'error (deg/s)'); xlabel(ha,xstr)
legend(ha,'roll','pitch','yaw'); %set(ha,'ylim',[-1 1]*30)

%K
ha=h(6);  i=index.K;  true=ekf.x(i,x)';  aposteriori=ekf.xhat(i,x)'; e=aposteriori;  s=ekf.S(i,x)*sg;  fcnploterrors3(ha,t,x,e,s)
title(ha,'{\itaposteriori} {\bfK} with 90% confidence bounds'); ylabel(ha,'value (pixels)'); xlabel(ha,xstr)
legend(ha,'xppo','yppo','focalLength');

%TP
ha=h(3);  s=zeros(3,nf);
for f=1:nf
    j = a.zi(:,f)==1;  nj = sum(j);  xi = fcnai2xi(a,j);
    s(:,f) = mean( reshape(ekf.S(xi,f),[3 nj]), 2); %sigmas
end
tx = ekf.xhat(index.tp0(:,1),:)*km2m;  tx(tx==0) = nan;
ty = ekf.xhat(index.tp0(:,2),:)*km2m;  ty(ty==0) = nan;
tz = ekf.xhat(index.tp0(:,3),:)*km2m;  tz(tz==0) = nan;
plot(ha,t,tx,'r',t,ty,'g',t,tz,'b');  axis(ha,'tight'); hold(ha,'on');  e=[]; fcnploterrors3(ha,t,x,e,s*km2m*sg)
title(ha,'{\itaposteriori} {\bfTP position} with 90% confidence bounds'); ylabel(ha,'NED Position (m)'); xlabel(ha,xstr)
%legend(ha,'x','y','z');

drawnow; 
end


function fcnploterrors3(ha,t,x,e,s)
if ~isempty(e);  plot(ha,t,e(x,1),'r',t,e(x,2),'g',t,e(x,3),'b');  hold(ha,'on'); end
plot(ha,t,s(1,x),'r',t,-s(1,x),'r','linewidth',2);  hold(ha,'on');
plot(ha,t,s(2,x),'g',t,-s(2,x),'g','linewidth',2)
plot(ha,t,s(3,x),'b',t,-s(3,x),'b','linewidth',2)

axis(ha,'tight')
box(ha,'off')
end


function xi = fcnai2xi(a,ai)
ntp = size(a.upx,1);
m = reshape(1:ntp*3, [3 ntp])+15;    xi= m(:,ai); 
xi = xi(:);
end