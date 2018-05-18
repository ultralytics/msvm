function [cam, rmse] = fcnMSV(cam,a)
startclock = clock;
rmse = [];
r2d = 180/pi;  d2r = pi/180;
ntp = size(a.upx,1); %number of tie points

flagsMIG = false;
fli = find(a.state(:,1)==1 & a.state(:,end)==1); %tie points exising in first and last frames


%BEST 3 ONLY
% [~, j] = max(a.score(fli));  fli = [fli(j); fli];
% [j(1), j(2), j(3)] = fcnget3distantframes([a.upx(fli,1) a.upy(fli,1) zeros(numel(fli),1)]);
% fli = fli(j);

if numel(fli)~=ntp
    flagsMIG = true;
    a0 = a;  a=fcncropTP(cam,a,fli);  ntp = size(a.upx,1);
end

[n1, n2, n3] = fcnget3distantframes(cam.apriori.ned);
% n1 = 1;
% n2 = 53;
% n3 = 105;
cam.msv.vf0 = [n1 n2 n3];

A=cam.apriori.ned(n1,:);  B=cam.apriori.ned(n2,:);  C=cam.apriori.ned(n3,:);
fprintf('Starting MSV: %3g points over frames [%i, %i, %i]...',ntp,n1,n2,n3) %figure; fcnplot3(cam.apriori.ned,'g.'); hold on; fcnplot3(cam.apriori.ned([n1 n2 n3],:),'bo')
st = [ntp 1];

[~, tpi] = max(a.score);  cam.msv.tp1i = tpi;
cam.u.xyz=cell(cam.frames,1); cam.uc.xyz=cam.u.xyz; zv=zeros(cam.frames,ntp); cam.u.x=zv; cam.u.y=zv; cam.u.z=zv; cam.uc.x=zv; cam.uc.y=zv; cam.uc.z=zv; a.upxi=zv'; a.upyi=zv';
for i=1:cam.frames
    [sc, uic] = pixel2camsc(cam, [a.upx(:,i) a.upy(:,i)]);  sc1 = sc(tpi,:)*d2r;  ui = fcnrotateB2W(0,sc1(2),sc1(3),uic);
    xy = camcc2pixel(cam,ui);  a.upxi(:,i) = xy(:,1);  a.upyi(:,i) = xy(:,2);
    cam.u.x(i,:) = ui(:,1);    cam.u.y(i,:) = ui(:,2);    cam.u.z(i,:) = ui(:,3);    cam.u.xyz{i} = ui;
    cam.uc.x(i,:) = uic(:,1);  cam.uc.y(i,:) = uic(:,2);  cam.uc.z(i,:) = uic(:,3);  cam.uc.xyz{i} = uic;
end

% %RECOVER ROLL ANGLE BETWEEN N1 and N2 -------------------------------------
% clc
% cam.true.rpy(n2,:)-cam.true.rpy(n1,:);
% f1=[a.upx(:,n1) a.upy(:,n1)];  f2=[a.upx(:,n2) a.upy(:,n2)];
% 
% Cw1 = fcnRPY2DCM_W2B(cam.true.rpy(n1,:)*d2r); C1w=Cw1';
% Cw3 = fcnRPY2DCM_W2B(cam.true.rpy(n2,:)*d2r);
% C13 = Cw3*C1w;  rpytrue = fcnW2BDCM2RPY(C13)*r2d
% 
% ov = ones(ntp,1);
% f1s = [f1 ov];
% f2s = [f2 ov];
% 
% % T = cp2tform(f2,f1,'projective'); % to overhead pixel
% % f1hat = tformfwd(T, f2);
% % f2hat = tforminv(T, f1);
% % 
% % T = cp2tform(f2,f1,'nonreflective similarity');
% % f1hat = tformfwd(T, f2);
% % f2hat = tforminv(T, f1);
% % 
% % u = [0 1];
% % v = [0 0];
% % [x, y] = tformfwd(T, u, v);
% % dx = x(2) - x(1);
% % dy = y(2) - y(1);
% % angle = (180/pi) * atan2(dy, dx);
% % scale = 1 / sqrt(dx^2 + dy^2);
% 
% %SFM
% Ccn = [0 0 1; 1 0 0; 0 -1 0];% DCM cam (z x -y) to ned (x y z)
% K = fcnK(cam.width,cam.height,cam.fovh,cam.fovv); %camera calibration matrix
% F = fcn8ptF(f1s', f2s');
% E = fcnF2E(F,K);
% P = fcnE2P(E);
% P = getCorrectCameraMatrix(P, K, K,f1s',f2s');
% R = Ccn*P(:,1:3)*Ccn'; %C1to3
% rpyhat = fcnW2BDCM2RPY(R)*r2d
% %Cw3hat = fcnW2BDCM2RPY(R*Cw1)*r2d
% t = P(:,4);

%RECOVER ROLL ANGLE BETWEEN N1 and N2 -------------------------------------
nx = cam.msvgridpoints;
nz = round(nx/2);
[xv,yv,zv,xva,yva,zva,z] = fcndefinegrid(cam,nx,nz);

np = numel(xva);
zm1 = zeros(size(xva));
urm = repmat(cam.u.xyz{n1},[np 1]);  cc1=[xva-A(1), yva-A(2), zm1]; cc1xys=(cc1(:,1).^2+cc1(:,2).^2)';
vrm = repmat(cam.u.xyz{n2},[np 1]);  cc2=[xva-B(1), yva-B(2), zm1]; cc2xys=(cc2(:,1).^2+cc2(:,2).^2)';
wrm = repmat(cam.u.xyz{n3},[np 1]);  cc3=[xva-C(1), yva-C(2), zm1]; cc3xys=(cc3(:,1).^2+cc3(:,2).^2)';

yaw1=repmat(fcnaz(cc1)',st);  yaw1=yaw1(:); cy1=cos(yaw1); sy1=sin(yaw1);
yaw2=repmat(fcnaz(cc2)',st);  yaw2=yaw2(:); cy2=cos(yaw2); sy2=sin(yaw2);
yaw3=repmat(fcnaz(cc3)',st);  yaw3=yaw3(:); cy3=cos(yaw3); sy3=sin(yaw3);

rollv = [0]*d2r;
nr = numel(rollv); %number of roll angles;
Uc = cell(nr,nz);
Vc = cell(nr,nz);
Wc = cell(nr,nz);
fxc = zeros(nr^3,1);
xc = cell(nr^3,1);
for i = 1:nr
    roll1=rollv(i);  sr1=sin(roll1);  cr1=cos(roll1);  %fprintf('%.0f\n',i)
    for zi = 1:nz;
        z1 = zv(zi) + zva;
        sp1 = (A(3)-z1)./sqrt(cc1xys + (z1-A(3)).^2);  cp1 = sqrt(1-sp1.*sp1); %sin pitch and cos pitch
        sp2 = (B(3)-z1)./sqrt(cc2xys + (z1-B(3)).^2);  cp2 = sqrt(1-sp2.*sp2);
        sp3 = (C(3)-z1)./sqrt(cc3xys + (z1-C(3)).^2);  cp3 = sqrt(1-sp3.*sp3);

        sp1=repmat(sp1,st);  cp1=repmat(cp1,st);  Uc{i,zi} = fcnrotateW2Bwsc(sr1,sp1(:),sy1,cr1,cp1(:),cy1,urm);
        sp2=repmat(sp2,st);  cp2=repmat(cp2,st);  Vc{i,zi} = fcnrotateW2Bwsc(sr1,sp2(:),sy2,cr1,cp2(:),cy2,vrm);
        sp3=repmat(sp3,st);  cp3=repmat(cp3,st);  Wc{i,zi} = fcnrotateW2Bwsc(sr1,sp3(:),sy3,cr1,cp3(:),cy3,wrm);
    end
end

n = 0;
for i = 1:nr
    for j = 1:nr
        for k = 1:nr
            n=n+1;  %fprintf('%3.0f',[i j k]); fprintf('--------------------------------------------------------------------------------------------------------------------------------------\n'); drawnow
            P = zeros(nx,nx,nz);
            for zi = 1:nz
                fx1 = fcnc03(A, Uc{i,zi}, B, Vc{j,zi}, C, Wc{k,zi});
                %fx1 = fcnc02(A, Uc{i,zi}, C, Wc{j,zi});
                fx = sum(reshape(fx1,[ntp np]));  P(:,:,zi) = reshape(fx, [nx nx])';
            end
            [fxc(n),mi]= min3(P);  if nz==1; mi(3)=1; end
            xc{n} = [xv(mi(2)) yv(mi(1)) zv(mi(3))+z(mi(2),mi(1)) rollv([i j k])]; %ML focus point (x-y swVERIFIED)
        end
    end
end
[~, i]=min(fxc);  mx=xc{i};  c0=mx(1:3);


%GE PLOT ------------------------------------------------------------------
zi = max(min(fcnindex1(zv,cam.DEM.Fned(0,0)/1000),nz),1);  lla = ned2lla(cam.DEM,[xv' yv' yv'*0]);
cam.google.kml.MSV = fcnGenerateKMLoverlayPNG(pwd, sprintf('%s MSV.png',cam.filename), lla(:,1), lla(:,2), exp(-P(:,:,zi)'));  %figure; pcolor(lla(:,2),lla(:,1),log(P(:,:,zi)')); shading flat

%ROLL ANGLES --------------------------------------------------------------
% fcnplot3D(xv,yv,zv,exp(-P),.02); 
% rolls = fcngetroll(cam,c0,A,C,n1,n3,st);
% rolls = rolls*d2r;

%DEFINE OPTIMIZER CONSTANTS
vf = unique(round( linspace(n1,n3,min(10,n3)) ));
nf=numel(vf);  ova=ones(1,nf);  ovb=ones(1,ntp);  cam.msv.vf=vf;  cam.msv.nf=nf;  cam.msv.ova=ova;  cam.msv.ovb=ovb;  cam.msv.ntp=numel(cam.u.x(1,:)); 
cam.msv.ux=cam.u.x(vf,:);  cam.msv.uy=cam.u.y(vf,:);  cam.msv.uz=cam.u.z(vf,:);  cam.msv.A=cam.aposteriori.ned(vf,:);
[~, cam.tpnedhat] = fcnminxyz(cam,a,c0);
if cam.syntheticVideoFlag; [rmse,error] = fcngeterrors(cam,a,vf); fprintf('%.1fm GRID RMSE (%.1fm std), ',rmse.xyz,mean(std(error.xyz))); end

%FCNMINXYZ ----------------------------------------------------------------
xhat1= lsqnonlin(@(x) fcnminxyz(cam,a,x),c0,[],[],cam.optim.options3);    [~, cam.tpnedhat] = fcnminxyz(cam,a,xhat1);
if cam.syntheticVideoFlag; [rmse,error] = fcngeterrors(cam,a,vf); fprintf('%.2fm NLS3 RMSE (%.2fm std), ',rmse.xyz,mean(std(error.xyz))); end

%DEFINE OPTIMIZER CONSTANTS
vf = unique(round( linspace(n1,n3,min(20,n3)) ));
nf=numel(vf);  ova=ones(1,nf);  ovb=ones(1,ntp);  cam.msv.vf=vf;  cam.msv.nf=nf;  cam.msv.ova=ova;  cam.msv.ovb=ovb;  cam.msv.ntp=numel(cam.u.x(1,:)); 
cam.msv.ux=cam.u.x(vf,:);  cam.msv.uy=cam.u.y(vf,:);  cam.msv.uz=cam.u.z(vf,:);  cam.msv.A=cam.aposteriori.ned(vf,:);

% %FCNMINXYZRRR -------------------------------------------------------------
% x0 = [xhat1 zeros(1,nf)];  [~, cam.tpnedhat] = fcnminxyzrrr(cam,a,x0);
% xhat2= lsqnonlin(@(x) fcnminxyzrrr(cam,a,x),x0,[],[],cam.optim.options3);    [~, cam.tpnedhat] = fcnminxyzrrr(cam,a,xhat2);
% if cam.syntheticVideoFlag; [rmse,error] = fcngeterrors(cam,a,vf); fprintf('%.2fm NLS5+ RMSE (%.2fm std), ',rmse.xyz,mean(std(error.xyz))); end
% 
% %RESECTION ----------------------------------------------------------------
% prxy = zeros(ntp,2,cam.frames); %pixel residuals xy
% for i=1:cam.frames
%     z = fcnvec2uvec(cam.tpnedhat - cam.aposteriori.ned(i*ovb,:));  %z = fcnvec2uvec(a.ipned - cam.true.ned(i*ovb,:))';
%     H=cam.uc.xyz{i};  %u = fcnvec2uvec(ned2camcc(cam,i,a.ipned,'true'));
%     DCM = (H'*H)\H'*z;  %DCM = DCM./sqrt(fcnrange(DCM)*fcnrange(DCM')');
%     cam.aposteriori.rpy(i,:) = fcnW2BDCM2RPY(DCM)*r2d;
%     cam.aposteriori.rpy(i,:) = fcnnlsresection(cam,a,i,cam.aposteriori.rpy(i,:))*r2d; %nonlinear resection
%     prxy(:,:,i) = [a.upx(:,i), a.upy(:,i)] - ned2pixel(cam,i,cam.tpnedhat,'aposteriori'); %z-zhat
% end
% cam.msv.residuals.pixelxy = prxy; %pixel residuals x-y
% cam.msv.residuals.pixelr = reshape( sqrt(prxy(:,1,:).^2 + prxy(:,2,:).^2), [ntp cam.frames]);  %pixel residuals range
% if cam.syntheticVideoFlag; rmse = fcngeterrors(cam,a,1:cam.frames); fprintf('%.2fdeg resection RMSE, ',rmse.rpy); end
% 
% %FCNMINRPYRPYRPY ----------------------------------------------------------
% x0 = cam.aposteriori.rpy(vf,:)'*d2r;
% xhat1 = lsqnonlin(@(x) fcnminrpyrpyrpy(cam,a,x),x0,[],[],cam.optim.options3);    [~, cam.tpnedhat] = fcnminrpyrpyrpy(cam,a,xhat1);
% if cam.syntheticVideoFlag; [rmse,error] = fcngeterrors(cam,a,vf); fprintf(' %.2fm RMSE (%.2fm std), ',rmse.xyz,mean(std(error.xyz))); end

%NLS RESECTION ------------------------------------------------------------
for i=1:cam.frames
    cam.aposteriori.rpy(i,:) = fcnnlsresection(cam,a,i,cam.aposteriori.rpy(i,:))*r2d;
    prxy(:,:,i) = [a.upx(:,i) a.upy(:,i)] - ned2pixel(cam,i,cam.tpnedhat,'aposteriori'); %z-zhat
end
if cam.syntheticVideoFlag; rmse = fcngeterrors(cam,a,1:cam.frames); fprintf('%.3fdeg in %.0fs',rmse.rpy,etime(clock,startclock)); end
fprintf(' Done (%.0fs). Residuals RMS = %.4f pixels\n',etime(clock,startclock),fcnrms(prxy))

%RUN MIG IF MISSING SOME GROUND POINTS ------------------------------------
if flagsMIG; cam.tpnedhat = fcnMIGMSV(cam,a0); end
if cam.syntheticVideoFlag; rmse = fcngeterrors(cam,a,1:cam.frames); fprintf('%.2fdeg MIG RMSE, ',rmse.rpy); end
fprintf(' Done (%.0fs). Residuals RMS = %.4f pixels\n',etime(clock,startclock),fcnrms(prxy))


%PLOT ---------------------------------------------------------------------
delete(findobj(0,'Tag','MSV'))
[ha, hf] = fig(1,2,1.98); set(hf,'Tag','MSV'); h=[];

axes(ha(1))
aned = cam.apriori.ned(1:cam.frames,:);
h(1)=fcnplot3(aned,'.b'); hold on;  h(2)=fcnplot3(cam.tpnedhat,'.r');  h(3)=plot3(aned([n1 n2 n3],1),aned([n1 n2 n3],2),aned([n1 n2 n3],3),'.b','markersize',15);
text(aned(n1,1),aned(n1,2),aned(n1,3),sprintf('  C_{%.0f}',n1),'fontsize',8);  text(aned(n2,1),aned(n2,2),aned(n2,3),sprintf(' C_{%.0f}',n2),'fontsize',8);  text(aned(n3,1),aned(n3,2),aned(n3,3),sprintf(' C_{%.0f}',n3),'fontsize',8)
axis tight;
title(sprintf('%s\nf(X_a) = %.2f, X_a = [%4.2f %4.2f %4.2f]km NED',cam.filename,min3(P)*1E3,cam.tpnedhat(tpi,:)))
n=1; ov=ones(n,1);
c = ov*cam.apriori.ned(n1,:); d=ov*cam.tpnedhat(tpi,:);  plot3([c(:,1) d(:,1)]',[c(:,2) d(:,2)]',[c(:,3) d(:,3)]','r-')
c = ov*cam.apriori.ned(n2,:); d=ov*cam.tpnedhat(tpi,:);  plot3([c(:,1) d(:,1)]',[c(:,2) d(:,2)]',[c(:,3) d(:,3)]','r-')
c = ov*cam.apriori.ned(n3,:); d=ov*cam.tpnedhat(tpi,:);  h(4) = plot3([c(:,1) d(:,1)]',[c(:,2) d(:,2)]',[c(:,3) d(:,3)]','r-');
sigma = 1000/cam.width*cam.fovh*d2r;
v = exp(-P/(2*sigma^2));
h(5:6)=fcnplot3D(xv,yv,zv,v,0.95,cam.tpnedhat(tpi,3));

zlim = fcnminmax(zv);  xlim = fcnminmax(xv);  ylim = fcnminmax(yv);
[~, i] = max3(v);  x=xv(i(2));  y=yv(i(1));  z=zv(i(3));
h(7) = plot3(xlim, [y y], [z z], '-k','LineWidth',1,'color',[.7 .7 .7]); hold on
plot3([x x], ylim, [z z], '-k','LineWidth',1,'color',[.7 .7 .7]);
plot3([x x], [y y], zlim, '-k','LineWidth',1,'color',[.7 .7 .7]); plot3(x, y, z, '.','MarkerSize',20,'Color',[.7 .7 .7]);

zlim = fcnminmax(zv);  xlim = fcnminmax(xv);  ylim = fcnminmax(yv);
x=cam.tpnedhat(tpi,1);  y=cam.tpnedhat(tpi,2);  z=cam.tpnedhat(tpi,3);
plot3(xlim, [y y], [z z], '-g','LineWidth',1);
plot3([x x], ylim, [z z], '-g','LineWidth',1);
h(7) = plot3([x x], [y y], zlim, '-g','LineWidth',1); plot3(x, y, z, '.g','MarkerSize',20);
%h(5) = surf(xv,yv,cam.tpnedhat(tpi,3)+(xv'*xv)*0,v(:,:,zi)); shading flat;

hl = legend(h,'Camera Positions C','Tie Points X','Cameras C_k','LOS Vectors L = C_i - X_a','X_a Likelihood',...
    'ML X_a (Zero Roll)','ML X_a (Resection+Triangulation)','Location','Best');
set(hl,'fontsize',7)
xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)'); view(-70,40); axis tight equal vis3d; set(gca,'cameraviewangle',10)

axes(ha(2)); h=[];
aned = cam.apriori.ned(1:cam.frames,:);
h(1)=fcnplot3(aned,'.b'); hold on;  h(2)=fcnplot3(cam.tpnedhat,'.r');  h(3)=plot3(aned([n1 n2 n3],1),aned([n1 n2 n3],2),aned([n1 n2 n3],3),'.b','markersize',15);
text(aned(n1,1),aned(n1,2),aned(n1,3),sprintf(' C_{%.0f}',n1),'fontsize',8);  text(aned(n2,1),aned(n2,2),aned(n2,3),sprintf(' C_{%.0f}',n2),'fontsize',8);  text(aned(n3,1),aned(n3,2),aned(n3,3),sprintf(' C_{%.0f}',n3),'fontsize',8)

x=0; y=0; z=0; onesMat = ones(2);
plot3(xlim, [y y], [z z], '-k','LineWidth',1,'color',[.7 .7 .7]); hold on
plot3([x x], ylim, [z z], '-k','LineWidth',1,'color',[.7 .7 .7]);
plot3([x x], [y y], zlim, '-k','LineWidth',1,'color',[.7 .7 .7]); plot3(x, y, z, '.','MarkerSize',20,'Color',[.7 .7 .7]);
surf(xlim,ylim,onesMat*z,'FaceColor',[.4 .4 .4],'EdgeColor',[.7 .7 .7]','CData', onesMat ,'facealpha',.1);

%PLOT 3D grid
%[xm,ym,zm]=ndgrid(xv,yv,zv);  h(8) = plot3(xm(:),ym(:),zm(:),'k.','MarkerSize',0.8);  hl=legend(h([1 2 3 8]),'Camera Positions C','Tie Points X','Cameras C_k','100x100x50 Grid Points','Location','Best');
hl=legend(h([1 2 3]),'Camera Positions C','Tie Points X','Cameras C_k','Location','Best'); %no grid
set(hl,'fontsize',7)

% %FIGURE 3
colors = rand(cam.msv.ntp,3);
n=min(100,cam.msv.ntp);
r = mean(fcnrange(mean(cam.tpnedhat,1),cam.apriori.ned))*2; %focal length
for j=[n1 n2 n3]
    c = cam.apriori.ned(j,:); 
    for i=1:n;
        rpy=cam.aposteriori.rpy(j,:)*d2r;
        d=c+fcnrotateW2B(rpy(1),rpy(2),rpy(3),r*[cam.uc.x(j,i)' cam.uc.y(j,i)' cam.uc.z(j,i)']);  plot3([c(:,1) d(:,1)]',[c(:,2) d(:,2)]',[c(:,3) d(:,3)]','-','Color',colors(i,:),'linewidth',.2)
        x=cam.tpnedhat(i,:); plot3(x(1),x(2),x(3),'.','Color',colors(i,:),'markersize',10);
    end
    d=c+fcnrotateW2B(rpy(1),rpy(2),rpy(3),r*[cam.uc.x(j,tpi)' cam.uc.y(j,tpi)' cam.uc.z(j,tpi)']);  plot3([c(:,1) d(:,1)]',[c(:,2) d(:,2)]',[c(:,3) d(:,3)]','r-','linewidth',2)
end
x=cam.tpnedhat(tpi,:); plot3(x(1),x(2),x(3),'r.','markersize',20);
surf(xv,yv,cam.tpnedhat(tpi,3)+(xv'*xv)*0,v(:,:,zi),'facealpha',.7); shading flat;
xyzlabel('x_{NED} (km)','y_{NED} (km)','z_{NED} (km)'); axis tight equal vis3d; set(gca,'cameraviewangle',10)
title('Initial Guess Region'); view(-70,40);
fcnpaperplots(cam)


ha = fig(1,1,2); 
popoutsubplot(gca,ha(1))

%FIGURE 4
i=1;
ntp = cam.msv.ntp;  zv = zeros(ntp,1)';
h=fig(2,1,2); axes(h(1));
hn=plot3([zv; cam.uc.x(i,:)],[zv; cam.uc.y(i,:)],[zv; cam.uc.z(i,:)],'-r'); hold on
ha=plot3([0; cam.uc.x(i,tpi)],[0; cam.uc.y(i,tpi)],[0; cam.uc.z(i,tpi)],'-b','linewidth',3); hold on
axis tight equal; box on
xlabel('x'); ylabel('y'); zlabel('z'); view(-165,20); 
set(h(1),'xlim',[0 1],'ylim',[-.05 .05],'zdir','reverse')
title('Image i=1, Normalized Image Coordinates X_i {\bfexpressed in Camera Reference Frame}')
legend([hn(1) ha(1)],'X_{i1...n} unit vectors','X_{ia} unit vector','location','northeast')

axes(h(2));
hn=plot3([zv; cam.u.x(i,:)],[zv; cam.u.y(i,:)],[zv; cam.u.z(i,:)],'-r'); hold on
ha=plot3([0; cam.u.x(i,tpi)],[0; cam.u.y(i,tpi)],[0; cam.u.z(i,tpi)],'-b','linewidth',3); hold on
axis tight equal; box on
xlabel('x'); ylabel('y'); zlabel('z'); view(-165,20);
set(h(2),'xlim',get(h(1),'xlim'),'ylim',get(h(1),'ylim'),'zlim',get(h(1),'zlim'),'zdir','reverse')
title('Image i=1, Normalized Image Coordinates X_i {\bfexpressed in X_a Reference Frame}')
legend([hn(1) ha(1)],'X_{i1...n} unit vectors','X_{ia} unit vector','location','northeast')


% %FIGURE 5, 3 images
% h=fig(3,1);
% try vfr=VideoReader([cam.filename '.avi']); catch vfr=VideoReader([cam.filename '.mp4']); end
% j = [n1 n2 n3];
% ki = [1 2 3 4];
% for i=1:3
%     ji = j(i);
%     rgb = single(read(vfr,ji))/255;  axes(h(i));
%     imshow(rgb,'parent',h(i));  title(sprintf('I_{%.0f}',ji)); hold on
%     plot(a.upx(ki(1),ji),a.upy(ki(1),ji),'r.'); plot(a.upx(ki(2),ji),a.upy(ki(2),ji),'g.'); plot(a.upx(ki(3),ji),a.upy(ki(3),ji),'b.'); plot(a.upx(ki(4),ji),a.upy(ki(4),ji),'m.')
%     axis on tight
% end
% fcnmarkersize(20)

%END
cam.apriori.rpy = cam.aposteriori.rpy;
cam.aposteriori.rpyMSV = cam.aposteriori.rpy;
drawnow
fprintf('\n')
end

function [xv,yv,zv,xva,yva,zva,z] = fcndefinegrid(cam,nx,nz)
demflag = true;
ned = cam.apriori.ned(1:cam.frames,:);
if demflag
    a = cam.DEM.Fned(ned(:,1),ned(:,2));
    a = -mean(a(~isnan(a)))/1000; %mean ned z value of DEM
else
    a = 0;
end

dx = a - mean(ned(:,3)) + 1;
xv = linspace(min(ned(:,1))-dx,        max(ned(:,1))+dx,     nx);
yv = linspace(min(ned(:,2))-dx,        max(ned(:,2))+dx,     nx);
zv = linspace(-1,                                     1,     nz)*.5;  if nz==1; zv=0; end
[x,y]=ndgrid(xv,yv);  xva=x(:); yva=y(:);

% if demflag
%     %zva = cam.DEM.Fned(xva,yva)'/1000;
%     zva = ones([1 nx^2])*a;
% else
     zva = zeros([1 nx^2]);
% end
z = reshape(zva,[nx nx]);

zv = zv + a;
end

function mx = fcngetroll(cam,c0,A,C,n1,n3,st)

py = fcnelaz(c0-A);  cpy = cos(py);  cp1= cpy(1);  cy1 = cpy(2);  spy = sin(py);  sp1= spy(1);  sy1 = spy(2);
py = fcnelaz(c0-C);  cpy = cos(py);  cp2= cpy(1);  cy2 = cpy(2);  spy = sin(py);  sp2= spy(1);  sy2 = spy(2);

nx = 500;
rollv = linspace(-20,20,nx)*d2r;
[rm1, rm2] = ndgrid(rollv,rollv);  

np = numel(rm1);
urm = repmat(cam.u.xyz{n1},[np 1]); 
wrm = repmat(cam.u.xyz{n3},[np 1]);

roll1 = rm1(:)';  cr1 = repmat(cos(roll1),st);  sr1 = repmat(sin(roll1),st);
roll2 = rm2(:)';  cr2 = repmat(cos(roll2),st);  sr2 = repmat(sin(roll2),st);

Uc = fcnrotateW2Bwsc(sr1(:),sp1,sy1,cr1(:),cp1,cy1,urm);
Wc = fcnrotateW2Bwsc(sr2(:),sp2,sy2,cr2(:),cp2,cy2,wrm);

fx1 = fcnc02(A, Uc, C, Wc);
fx = sum(reshape(fx1,[st(1) np]));  Pr = reshape(fx, [nx nx])';
figure; pcolor(rm1,rm2,(Pr)); shading flat

[mx, mi] = min3(Pr);  mx = [rollv(mi(1)) rollv(mi(2))]*r2d;
end

%LOCAL FUNCTIONS ----------------------------------------------------------
function x = fcnnlsresection(cam,a,i,x)
fz = pixel2camsc(cam, [a.upx(:,i) a.upy(:,i)]); fz=fz(:,2:3)*d2r;
x=x(:)*d2r; %x = cam.apriori.rpy(i,:)';
dx_ned = cam.tpnedhat-cam.aposteriori.ned(i*cam.msv.ovb,:);  TAx=dx_ned(:,1); TAy=dx_ned(:,2); TAz=dx_ned(:,3);  deldrpy=zeros(cam.msv.ntp,3);  dazdrpy=deldrpy;
for j=1:4
    r = x(1);  sr=sin(r); cr=cos(r);
    p = x(2);  sp=sin(p); cp=cos(p);
    y = x(3);  sy=sin(y); cy=cos(y);
    k1 = (TAx*(sr*sy + cr*cy*sp) - TAy*(cy*sr - cr*sp*sy) + TAz*cp*cr);
    k2 = (TAy*(cr*cy + sp*sr*sy) - TAx*(cr*sy - cy*sp*sr) + TAz*cp*sr);
    k3 = (TAx*cp*cy - TAz*sp + TAy*cp*sy);
    k4 = (k3.^2 + k1.^2 + k2.^2);
    fhat1 = fcnelaz(fcnrotateB2W(r,p,y,dx_ned));
%     dx = 1E-8;  f0=fcnelaz(dx_ned*fcnRPY2DCM_B2W(x+[0 0 0]'));  f1=fcnelaz(dx_ned*fcnRPY2DCM_B2W(x+[dx 0 0]'));  f2=fcnelaz(dx_ned*fcnRPY2DCM_B2W(x+[0 dx 0]'));  f3=fcnelaz(dx_ned*fcnRPY2DCM_B2W(x+[0 0 dx]'));
%     Bnumerical = [f1(:)-f0(:) f2(:)-f0(:) f3(:)-f0(:)]/dx;
    deldrpy(:,1) = k2./sqrt(k4 - k1.^2);
    deldrpy(:,2) = -((TAx*cp*cr*cy - TAz*cr*sp + TAy*cp*cr*sy)./sqrt(k4) - (k1.*(2*k2.*(TAx*cp*cy*sr - TAz*sp*sr + TAy*cp*sr*sy) - 2*k3.*(TAz*cp + TAx*cy*sp + TAy*sp*sy) + 2*k1.*(TAx*cp*cr*cy - TAz*cr*sp + TAy*cp*cr*sy)))./(2*k4.^(3./2)))./sqrt(1 - k1.^2./k4);
    deldrpy(:,3) = -((TAx*(cy*sr - cr*sp*sy) + TAy*(sr*sy + cr*cy*sp))./sqrt(k4) - (k1.*(2*(TAx*(cy*sr - cr*sp*sy) + TAy*(sr*sy + cr*cy*sp)).*k1 + 2*(TAy*cp*cy - TAx*cp*sy).*k3 - 2*(TAx*(cr*cy + sp*sr*sy) + TAy*(cr*sy - cy*sp*sr)).*k2))./(2*k4.^(3./2)))./sqrt(1 - k1.^2./k4);
    dazdrpy(:,1) = k1./((k2.^2./k3.^2 + 1).*k3);
    dazdrpy(:,2) = ((TAx*cp*cy*sr - TAz*sp*sr + TAy*cp*sr*sy)./k3 + (k2.*(TAz*cp + TAx*cy*sp + TAy*sp*sy))./k3.^2)./(k2.^2./k3.^2 + 1);
    dazdrpy(:,3) = -((TAx*(cr*cy + sp*sr*sy) + TAy*(cr*sy - cy*sp*sr))./k3 + ((TAy*cp*cy - TAx*cp*sy).*k2)./k3.^2)./(k2.^2./k3.^2 + 1);
    B = [deldrpy; dazdrpy];

    f = fz - fhat1;
    Bt = B';
    x = x+(Bt*B)\(Bt*f(:));
end
end

function [rmse, error] = fcngeterrors(cam,a,i)
if isfield(a,'iis') %include only specific features in error statistics
    error.xyz = (cam.tpnedhat(a.iis,:)-a.ipned(a.iis,:))*1000;
else
    error.xyz = (cam.tpnedhat-a.ipned)*1000;
end
error.rpy = fcndrpyd(cam.true.rpy(i,:),cam.aposteriori.rpy(i,:));
rmse.xyz = fcnrms(error.xyz);
rmse.rpy = fcnrms(error.rpy); %tie point RMSE, rpy RMSE
rmse.each = [fcnrms(error.xyz(:,1)) fcnrms(error.xyz(:,2)) fcnrms(error.xyz(:,3)) fcnrms(error.rpy(:,1)) fcnrms(error.rpy(:,2)) fcnrms(error.rpy(:,3))];
end

function x2 = fcnrotateW2B(r,p,y,x)
sr=sin(r);  sp=sin(p);  sy=sin(y);
cr=cos(r);  cp=cos(p);  cy=cos(y);
x2 = zeros(size(x));
x2(:,1)=x(:,1).*(cp.*cy)+x(:,2).*(sr.*sp.*cy-cr.*sy)+x(:,3).*(cr.*sp.*cy+sr.*sy);
x2(:,2)=x(:,1).*(cp.*sy)+x(:,2).*(sr.*sp.*sy+cr.*cy)+x(:,3).*(cr.*sp.*sy-sr.*cy);
x2(:,3)=x(:,1).*(-sp)+x(:,2).*(sr.*cp)+x(:,3).*(cr.*cp);
end

function [x2, y2, z2] = fcnrotateW2Bxyz(r,p,y,x)
sr=sin(r);  sp=sin(p);  sy=sin(y);
cr=cos(r);  cp=cos(p);  cy=cos(y);
x2=x(:,1).*(cp.*cy)+x(:,2).*(sr.*sp.*cy-cr.*sy)+x(:,3).*(cr.*sp.*cy+sr.*sy);
y2=x(:,1).*(cp.*sy)+x(:,2).*(sr.*sp.*sy+cr.*cy)+x(:,3).*(cr.*sp.*sy-sr.*cy);
z2=x(:,1).*(-sp)+x(:,2).*(sr.*cp)+x(:,3).*(cr.*cp);
end

function [x2, y2, z2] = fcnrotateW2Brpyxyz(sr,sp,sy,cr,cp,cy,x,y,z)
x2=x.*(cp.*cy) +y.*(sr.*sp.*cy-cr.*sy)  +z.*(cr.*sp.*cy+sr.*sy);
y2=x.*(cp.*sy) +y.*(sr.*sp.*sy+cr.*cy)  +z.*(cr.*sp.*sy-sr.*cy);
z2=x.*(-sp)    +y.*(sr.*cp)             +z.*(cr.*cp);
end

function x2 = fcnrotateW2Bwsc(sr,sp,sy,cr,cp,cy,xin) %with sin and cosine supplied
x2 = zeros(size(xin));
x2(:,1)=xin(:,1).*(cp.*cy)+xin(:,2).*(sr.*sp.*cy-cr.*sy)+xin(:,3).*(cr.*sp.*cy+sr.*sy);
x2(:,2)=xin(:,1).*(cp.*sy)+xin(:,2).*(sr.*sp.*sy+cr.*cy)+xin(:,3).*(cr.*sp.*sy-sr.*cy);
x2(:,3)=xin(:,1).*(-sp)   +xin(:,2).*(sr.*cp)           +xin(:,3).*(cr.*cp);
end

% function x2 = fcnrotateB2W(r,p,y,x)
% sr=sin(r); sp=sin(p); sy=sin(y);
% cr=cos(r); cp=cos(p); cy=cos(y);
% x2 = zeros(size(x));
% x2(:,1)=x(:,1).*(cp.*cy)+x(:,2).*(cp.*sy )+x(:,3).*(-sp);
% x2(:,2)=x(:,1).*(sr.*sp.*cy-cr.*sy)+x(:,2).*(sr.*sp.*sy+cr.*cy)+x(:,3).*(sr.*cp);
% x2(:,3)=x(:,1).*(cr.*sp.*cy+sr.*sy)+x(:,2).*(cr.*sp.*sy-sr.*cy)+x(:,3).*(cr.*cp);
% end

function [x2, y2, z2] = fcnrotateB2Wxyz(r,p,y,x)
sr=sin(r); sp=sin(p); sy=sin(y);
cr=cos(r); cp=cos(p); cy=cos(y);
x2=x(:,1).*(cp.*cy)+x(:,2).*(cp.*sy )+x(:,3).*(-sp);
y2=x(:,1).*(sr.*sp.*cy-cr.*sy)+x(:,2).*(sr.*sp.*sy+cr.*cy)+x(:,3).*(sr.*cp);
z2=x(:,1).*(cr.*sp.*cy+sr.*sy)+x(:,2).*(cr.*sp.*sy-sr.*cy)+x(:,3).*(cr.*cp);
end

function [res, C0] = fcnminrpyrpyrpy(cam,a,x) %[rpyrpyrpy...] CENTERED
%LOAD CONSTANTS
vf = cam.msv.vf;  nf = cam.msv.nf;  ntp = cam.msv.ntp;  ovb = cam.msv.ovb;  A = cam.msv.A;

%RPY
rpy = reshape(x,[3 nf])';

%VECTOR ROTATIONS
r = rpy(:,1); cr=cos(r); sr=sin(r); cr=cr(:,ovb); sr=sr(:,ovb);
p = rpy(:,2); cp=cos(p); sp=sin(p); cp=cp(:,ovb); sp=sp(:,ovb);
y = rpy(:,3); cy=cos(y); sy=sin(y); cy=cy(:,ovb); sy=sy(:,ovb);
x0=cam.uc.x(vf,:); y0=cam.uc.y(vf,:); z0=cam.uc.z(vf,:);%x0=cam.msv.ux; y0=cam.msv.uy; z0=cam.msv.uz;
ux1=x0.*(cp.*cy) +y0.*(sr.*sp.*cy-cr.*sy)  +z0.*(cr.*sp.*cy+sr.*sy);  %[ux1,uy1,uz1] = fcnrotateW2Brpyxyz(sr,sp,sy,cr,cp,cy,x,y,z);
uy1=x0.*(cp.*sy) +y0.*(sr.*sp.*sy+cr.*cy)  +z0.*(cr.*sp.*sy-sr.*cy);
uz1=x0.*(-sp)    +y0.*(sr.*cp)             +z0.*(cr.*cp);

C0 = fcnNvintercept(A,ux1,uy1,uz1);
C0x=C0(:,1)';  C0y=C0(:,2)';  C0z=C0(:,3)';

%IMAGE ANGLE RESIDUALS
%ova = cam.msv.ova;  CAx = C0x(ova,:)-A(:,ovb);  CAy = C0y(ova,:)-A(:,ovb*2);  CAz = C0z(ova,:)-A(:,ovb*3);
%ts1 = cam.msv.nf*ntp - sum(sum(  (ux1.*CAx + uy1.*CAy + uz1.*CAz).^2./(CAx.*CAx+CAy.*CAy+CAz.*CAz)  ));  %sum of sqrt of angle sines  %sin = sqrt(1-ct^2)
%ts1 = asin(sqrt(abs(  1 - (ux1.*CAx + uy1.*CAy + uz1.*CAz).^2./(CAx.*CAx+CAy.*CAy+CAz.*CAz)  )));  %sum of sqrt of angle sines  %sin = sqrt(1-ct^2)
%ts1 = 1 - (ux1.*CAx + uy1.*CAy + uz1.*CAz).^2./(CAx.*CAx+CAy.*CAy+CAz.*CAz);  %sum of sqrt of angle sines  %sin = sqrt(1-ct^2)
%res = ts1;

%IMAGE PIXEL RESIDUALS
zx = [a.upx(:,vf); a.upy(:,vf)];
res = zeros(ntp*2,nf); %residuals
for i=1:nf
    ui = [C0x-A(i,1); C0y-A(i,2); C0z-A(i,3)]';
    uic = fcnrotateB2W(r(i),p(i),y(i),ui);
    z1 = camcc2pixel(cam, uic);
    res(:,i) = zx(:,i)-z1(:); 
end
end

function [res, C0] = fcnminxyzrrr(cam,a,x) %[xyzrrr...]
%LOAD CONSTANTS
vf = cam.msv.vf;  nf = cam.msv.nf;  ntp = cam.msv.ntp;  ovb = cam.msv.ovb;  A = cam.msv.A;

%RPY
r = x(4:3+nf)'; %roll
xt = x(1)-A(:,1);  yt = x(2)-A(:,2);  zt = x(3)-A(:,3);  y = 2*atan((sqrt(xt.^2+yt.^2)-xt)./yt);
sp = -zt./sqrt(xt.^2 + yt.^2 + zt.^2);  cp = sqrt(1 - sp.^2);  p = asin(-zt./sqrt(xt.^2 + yt.^2 + zt.^2));

% %HANK STUFF ---------------------------------------------------------------
% vf=cam.msv.vf; x=cam.uc.x(vf,:); y=cam.uc.y(vf,:); z=cam.uc.z(vf,:);
% rpy=zeros(cam.msv.nf,3);
% cam.tpnedhat(1,:)=x(1:3);
% cam.aposteriori.rpy(vf,1)=r*r2d;
% for i=1:cam.msv.nf
%     [~,~,Mned] = solveAzPt_givenRoll(cam,a,cam.msv.vf(i));
%     rpy(i,:) = fcnB2WDCM2RPY(Mned); %compare to cam.true.rpy(vf,:) 
% end
% r=rpy(:,1);     cr=cos(r); sr=sin(r);       cr=cr(:,ovb); sr=sr(:,ovb);
% p=rpy(:,2);     cp=cos(p); sp=sin(p);       cp=cp(:,ovb); sp=sp(:,ovb);
% yaw=rpy(:,3);   cy=cos(yaw); sy=sin(yaw);   cy=cy(:,ovb); sy=sy(:,ovb);
% %HANK STUFF ---------------------------------------------------------------

%VECTOR ROTATIONS
cr=cos(r); sr=sin(r); cr=cr(:,ovb); sr=sr(:,ovb);
                      cp=cp(:,ovb); sp=sp(:,ovb);
cy=cos(y); sy=sin(y); cy=cy(:,ovb); sy=sy(:,ovb);
x0=cam.msv.ux; y0=cam.msv.uy; z0=cam.msv.uz;
ux1=x0.*(cp.*cy) +y0.*(sr.*sp.*cy-cr.*sy)  +z0.*(cr.*sp.*cy+sr.*sy);  %[ux1,uy1,uz1] = fcnrotateW2Brpyxyz(sr,sp,sy,cr,cp,cy,x,y,z);
uy1=x0.*(cp.*sy) +y0.*(sr.*sp.*sy+cr.*cy)  +z0.*(cr.*sp.*sy-sr.*cy);
uz1=x0.*(-sp)    +y0.*(sr.*cp)             +z0.*(cr.*cp);

C0 = fcnNvintercept(A,ux1,uy1,uz1);  %C0 = fcn2vintercept(A,ux1,uy1,uz1);
C0x=C0(:,1)';  C0y=C0(:,2)';  C0z=C0(:,3)';

%IMAGE ANGLE RESIDUALS
% ova = cam.msv.ova;  CAx = C0x(ova,:)-A(:,ovb);  CAy = C0y(ova,:)-A(:,ovb*2);  CAz = C0z(ova,:)-A(:,ovb*3);
% ts1 = cam.msv.nf*ntp - sum(sum(  (ux1.*CAx + uy1.*CAy + uz1.*CAz).^2./(CAx.*CAx+CAy.*CAy+CAz.*CAz)  ));  %sum of sqrt of angle sines  %sin = sqrt(1-ct^2)
% ts1 = asin(sqrt(abs(  1 - (ux1.*CAx + uy1.*CAy + uz1.*CAz).^2./(CAx.*CAx+CAy.*CAy+CAz.*CAz)  )));  %sum of sqrt of angle sines  %sin = sqrt(1-ct^2)
% ts1 = 1 - (ux1.*CAx + uy1.*CAy + uz1.*CAz).^2./(CAx.*CAx+CAy.*CAy+CAz.*CAz);  %sum of sqrt of angle sines  %sin = sqrt(1-ct^2)
% res = ts1;

% %IMAGE PIXEL RESIDUALS
zx = [a.upxi(:,vf); a.upyi(:,vf)];
res = zeros(ntp*2,nf); %residuals
for i=1:nf
    ui = [C0x-A(i,1); C0y-A(i,2); C0z-A(i,3)]'; 
    uic = fcnrotateB2W(r(i),p(i),y(i),ui);
    z1 = camcc2pixel(cam, uic);
    res(:,i) = zx(:,i)-z1(:); 
end
end


function [res, C0] = fcnminxyz(cam,a,x) %[xyz]
%LOAD CONSTANTS
vf = cam.msv.vf;  nf = cam.msv.nf;  ntp = cam.msv.ntp;  ovb = cam.msv.ovb;  A = cam.msv.A;

%RPY
r = zeros(nf,1); %roll
xt = x(1)-A(:,1);  yt = x(2)-A(:,2);  zt = x(3)-A(:,3);  y = 2*atan((sqrt(xt.^2+yt.^2)-xt)./yt);
sp = -zt./sqrt(xt.^2 + yt.^2 + zt.^2);  cp = sqrt(1 - sp.^2);  p = asin(-zt./sqrt(xt.^2 + yt.^2 + zt.^2));

%VECTOR ROTATIONS
cr=cos(r); sr=sin(r); cr=cr(:,ovb); sr=sr(:,ovb);
                      cp=cp(:,ovb); sp=sp(:,ovb);
cy=cos(y); sy=sin(y); cy=cy(:,ovb); sy=sy(:,ovb);
x0=cam.msv.ux; y0=cam.msv.uy; z0=cam.msv.uz;
ux1=x0.*(cp.*cy) +y0.*(sr.*sp.*cy-cr.*sy)  +z0.*(cr.*sp.*cy+sr.*sy);  %[ux1,uy1,uz1] = fcnrotateW2Brpyxyz(sr,sp,sy,cr,cp,cy,x,y,z);
uy1=x0.*(cp.*sy) +y0.*(sr.*sp.*sy+cr.*cy)  +z0.*(cr.*sp.*sy-sr.*cy);
uz1=x0.*(-sp)    +y0.*(sr.*cp)             +z0.*(cr.*cp);

C0 = fcnNvintercept(A,ux1,uy1,uz1);  %C0 = fcn2vintercept(A,ux1,uy1,uz1);
C0x=C0(:,1)';  C0y=C0(:,2)';  C0z=C0(:,3)';

%IMAGE ANGLE RESIDUALS
%ova = cam.msv.ova;  CAx = C0x(ova,:)-A(:,ovb);  CAy = C0y(ova,:)-A(:,ovb*2);  CAz = C0z(ova,:)-A(:,ovb*3);
%ts1 = cam.msv.nf*ntp - sum(sum(  (ux1.*CAx + uy1.*CAy + uz1.*CAz).^2./(CAx.*CAx+CAy.*CAy+CAz.*CAz)  ));  %sum of sqrt of angle sines  %sin = sqrt(1-ct^2)
%ts1 = asin(sqrt(abs(  1 - (ux1.*CAx + uy1.*CAy + uz1.*CAz).^2./(CAx.*CAx+CAy.*CAy+CAz.*CAz)  )));  %sum of sqrt of angle sines  %sin = sqrt(1-ct^2)
%ts1 = 1 - (ux1.*CAx + uy1.*CAy + uz1.*CAz).^2./(CAx.*CAx+CAy.*CAy+CAz.*CAz);  %sum of sqrt of angle sines  %sin = sqrt(1-ct^2)
%res = ts1;

% %IMAGE PIXEL RESIDUALS
zx = [a.upxi(:,vf); a.upyi(:,vf)];
res = zeros(ntp*2,nf); %residuals
for i=1:nf
    ui = [C0x-A(i,1); C0y-A(i,2); C0z-A(i,3)]'; 
    uic = fcnrotateB2W(r(i),p(i),y(i),ui);
    z1 = camcc2pixel(cam, uic);
    res(:,i) = zx(:,i)-z1(:); 
end
end

function rs = fcnc03(A,u,B,v,C,w)
x=A-B;  y=A-C;  z=B-C;

%A-B
d = u(:,1).*v(:,1) + u(:,2).*v(:,2) + u(:,3).*v(:,3);
e = u(:,1).*x(1) + u(:,2).*x(2) + u(:,3).*x(3);
f = v(:,1).*x(1) + v(:,2).*x(2) + v(:,3).*x(3);
g = 1 - d.*d;
s1 = (d.*f - e)./g;
t1 = (f - d.*e)./g;
%rs = (t1.*v(:,1)-x(1)-s1.*u(:,1)).^2 + (t1.*v(:,2)-x(2)-s1.*u(:,2)).^2 + (t1.*v(:,3)-x(3)-s1.*u(:,3)).^2;

%A-C
d = u(:,1).*w(:,1) + u(:,2).*w(:,2) + u(:,3).*w(:,3);
e = u(:,1).*y(1) + u(:,2).*y(2) + u(:,3).*y(3);
f = w(:,1).*y(1) + w(:,2).*y(2) + w(:,3).*y(3);
g = 1 - d.*d;
s2 = (d.*f - e)./g;
t2 = (f - d.*e)./g;
%rs = rs + (t2.*w(:,1)-y(1)-s2.*u(:,1)).^2 + (t2.*w(:,2)-y(2)-s2.*u(:,2)).^2 + (t2.*w(:,3)-y(3)-s2.*u(:,3)).^2;

%B-C
d = v(:,1).*w(:,1) + v(:,2).*w(:,2) + v(:,3).*w(:,3);
e = v(:,1).*z(1) + v(:,2).*z(2) + v(:,3).*z(3);
f = w(:,1).*z(1) + w(:,2).*z(2) + w(:,3).*z(3);
g = 1 - d.*d;
s3 = (d.*f - e)./g;
t3 = (f - d.*e)./g;
%rs = rs + (t3.*w(:,1)-z(1)-s3.*v(:,1)).^2 + (t3.*w(:,2)-z(2)-s3.*v(:,2)).^2 + (t3.*w(:,3)-z(3)-s3.*v(:,3)).^2;

D = (A+B+C)*2;  t1s3=t1+s3;  s1s2=s1+s2;  t2t3=t2+t3;
C0x = (D(1) + t1s3.*v(:,1)+s1s2.*u(:,1) + t2t3.*w(:,1))*(1/6);
C0y = (D(2) + t1s3.*v(:,2)+s1s2.*u(:,2) + t2t3.*w(:,2))*(1/6);
C0z = (D(3) + t1s3.*v(:,3)+s1s2.*u(:,3) + t2t3.*w(:,3))*(1/6);

%ALTERNATE
%C0 = fcn2vintercept([A; B],[u(:,1) v(:,1)]',[u(:,2) v(:,2)]',[u(:,3) v(:,3)]');
%C0x=C0(:,1);  C0y=C0(:,2);  C0z=C0(:,3);

%ANGLE RESIDUALS
ct1sq = (u(:,1).*(C0x-A(1)) + u(:,2).*(C0y-A(2)) + u(:,3).*(C0z-A(3))).^2./((C0x-A(1)).^2+(C0y-A(2)).^2+(C0z-A(3)).^2);  %sin = sqrt(1-ct^2)
ct2sq = (v(:,1).*(C0x-B(1)) + v(:,2).*(C0y-B(2)) + v(:,3).*(C0z-B(3))).^2./((C0x-B(1)).^2+(C0y-B(2)).^2+(C0z-B(3)).^2);
ct3sq = (w(:,1).*(C0x-C(1)) + w(:,2).*(C0y-C(2)) + w(:,3).*(C0z-C(3))).^2./((C0x-C(1)).^2+(C0y-C(2)).^2+(C0z-C(3)).^2);  

%rs = 3-(ct1sq.^2+ct2sq.^2+ct3sq.^2); %sum of sqrt of angle sines  %ts4 = sqrt(ts1)+sqrt(ts2)+sqrt(ts3); %sum of angle sines
rs = acos(sqrt(min(ct1sq,1))).^2 + acos(sqrt(min(ct2sq,1))).^2 + acos(sqrt(min(ct3sq,1))).^2; %sum of sqrt of angle sines  %ts4 = sqrt(ts1)+sqrt(ts2)+sqrt(ts3); %sum of angle sines
end

function rs = fcnc02(A,u,B,v)
x=A-B;

%A-B
d = u(:,1).*v(:,1) + u(:,2).*v(:,2) + u(:,3).*v(:,3);
e = u(:,1).*x(1) + u(:,2).*x(2) + u(:,3).*x(3);
f = v(:,1).*x(1) + v(:,2).*x(2) + v(:,3).*x(3);
g = 1 - d.*d;
s1 = (d.*f - e)./g;
t1 = (f - d.*e)./g;
%rs = (t1.*v(:,1)-x(1)-s1.*u(:,1)).^2 + (t1.*v(:,2)-x(2)-s1.*u(:,2)).^2 + (t1.*v(:,3)-x(3)-s1.*u(:,3)).^2;

D = (A+B);
C0x = (D(1) + t1.*v(:,1)+s1.*u(:,1))*(1/2);
C0y = (D(2) + t1.*v(:,2)+s1.*u(:,2))*(1/2);
C0z = (D(3) + t1.*v(:,3)+s1.*u(:,3))*(1/2);

%ALTERNATE
%C0 = fcn2vintercept([A; B],[u(:,1) v(:,1)]',[u(:,2) v(:,2)]',[u(:,3) v(:,3)]');
%C0x=C0(:,1);  C0y=C0(:,2);  C0z=C0(:,3);

%ANGLE RESIDUALS
dot1 = (u(:,1).*(C0x-A(1)) + u(:,2).*(C0y-A(2)) + u(:,3).*(C0z-A(3))).^2./((C0x-A(1)).^2+(C0y-A(2)).^2+(C0z-A(3)).^2);  
dot2 = (v(:,1).*(C0x-B(1)) + v(:,2).*(C0y-B(2)) + v(:,3).*(C0z-B(3))).^2./((C0x-B(1)).^2+(C0y-B(2)).^2+(C0z-B(3)).^2);  

%angle1 = acos(dot1); %sin = sqrt(1-ct^2)
%angle2 = acos(dot2);  
%rs = real(angle1+angle2);
rs = acos(dot1).^2 + acos(dot2).^2;
%rs = (1-dot1).^2 + (1-dot2).^2;
end

function h=fcnplot3D(xv,yv,zv,v,percentile,greenv,redv)
isovalue          = fcndsearch(v,percentile);
[ml.val, i]       = max3(v);  ml.row=i(1);  ml.col=i(2);  ml.layer=i(3);
onesMat = ones(2);  zlim = fcnminmax(zv);  xlim = fcnminmax(xv);  ylim = fcnminmax(yv);

%plot mlpoint
x=xv(ml.col);  y=yv(ml.row);  z=zv(ml.layer);
plot3(xlim, [y y], [z z], '-k','LineWidth',1,'color',[.7 .7 .7]); hold on
plot3([x x], ylim, [z z], '-k','LineWidth',1,'color',[.7 .7 .7]);
plot3([x x], [y y], zlim, '-k','LineWidth',1,'color',[.7 .7 .7]); plot3(x, y, z, '.','MarkerSize',20,'Color',[.7 .7 .7]);
h(1) = surf(xlim,ylim,onesMat*z,'FaceColor',[.4 .4 .4],'EdgeColor',[.7 .7 .7]','CData', onesMat ); alpha(h(1),0.1);

[xm,ym,zm]=meshgrid(xv,yv,zv);
n = isosurface(xm,ym,zm,v,isovalue);
p1 = patch(n,'AmbientStrength',1);
set(p1,'FaceColor',[.5 .5 .5],'EdgeColor','none','FaceAlpha',.85)
%reducepatch(p1, .8)
isonormals(xv,yv,zv,v,p1) %WTF does this actually do?
lighting phong;  material shiny
mrp = mean(zv);%mean reactor power
mrp = mrp + randn*mrp/3;
light('Position',[0 0 mrp-100]);  light('Position',[0 0 mrp+100]);  light('Position',[-500 0 mrp]);  light('Position',[0 -500 mrp]);  light('Position',[500 0 mrp]);  light('Position',[0 500 mrp]);  camlight
p2 = patch(isocaps(xm,ym,zm,v,isovalue));
set(p2,'EdgeColor','none','FaceColor','interp','SpecularColorReflectance',0,'SpecularStrength',.3);

%PLOT MAX LIKELIHOOD CONTOUR
hcont = contourslice(xm,ym,zm,v,[],[],zv(ml.layer), [isovalue isovalue]);
set(hcont,'EdgeColor',[.6 .6 .6],'LineWidth',1)

if exist('redv','var')
    hcont = contourslice(xm,ym,zm,v,[],[],redv, [isovalue isovalue]);
    set(hcont,'EdgeColor',[1 0 0],'LineWidth',1)
end
if exist('greenv','var')
    hcont = contourslice(xm,ym,zm,v,[],[],greenv, [isovalue isovalue]);
    set(hcont,'EdgeColor',[0 1 0],'LineWidth',1)
end
fcnfontsize; box on; axis equal vis3d tight; set(gca,'zdir','reverse')
h(2) = p1;
end


