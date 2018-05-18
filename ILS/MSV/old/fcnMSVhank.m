function [cam, rmse] = fcnMSVhank(cam,a)
activetp = a.upx~=0 | a.upy~=0;
ntp = find(any(a.upy~=0 & a.upx~=0,2),1,'last'); %number of tie points
aned = cam.true.ned(1:cam.frames,:);
r2d = 180/pi;
d2r = pi/180;
t = [ntp 1];

n1 = 1;%ceil(rand*(cam.frames-50));
n2 = round(cam.frames/2);
n3 = cam.frames;%cam.frames;
A=cam.apriori.ned(n1,:);  B=cam.apriori.ned(n2,:);  C=cam.apriori.ned(n3,:);

%[~, tpi] = min(sum([a.upx(:,[n1 n2 n3])-cam.width/2, a.upy(:,[n1 n2 n3])-cam.height/2].^2,2)); %best tp
tpi = 1; 
cam.u.xyz=cell(cam.frames,1); cam.uc.xyz=cam.u.xyz; zv=zeros(cam.frames,ntp); cam.u.x=zv; cam.u.y=zv; cam.u.z=zv; cam.uc.x=zv; cam.uc.y=zv; cam.uc.z=zv;
for i=1:cam.frames
    f1=[a.upx(:,i)-a.upx(tpi,i)+cam.width/2 a.upy(:,i)-a.upy(tpi,i)+cam.height/2]+1/2;
    ui=fcnSC2CCd(pixel2camsc(cam, f1));
    cam.u.x(i,:)=ui(:,1);  cam.u.y(i,:)=ui(:,2);  cam.u.z(i,:)=ui(:,3);  cam.u.xyz{i}=ui;
    
    ui=fcnSC2CCd(pixel2camsc(cam, [a.upx(:,i) a.upy(:,i)]));
    cam.uc.x(i,:)=ui(:,1);  cam.uc.y(i,:)=ui(:,2);  cam.uc.z(i,:)=ui(:,3);  cam.uc.xyz{i}=ui;
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
% % T = cp2tform(f2,f1,'projective'); %pixel to overhead pixel
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

%t = P(:,4);
%RECOVER ROLL ANGLE BETWEEN N1 and N2 -------------------------------------

P0 = mean(cam.apriori.ned(1:cam.frames,:));  %P0 = [0 0 0];
dx = 3;
nx = cam.msvgridpoints;
nz = cam.msvgridpoints/2;
xv = linspace(P0(1)-dx,                        P0(1)+dx,     nx);
yv = linspace(P0(2)-dx,                        P0(2)+dx,     nx);
zv = linspace(min(cam.apriori.ned(:,3))-.1,    1,            nz);
[x,y]=ndgrid(xv,yv);  xva=x(:); yva=y(:);

np = numel(x);
zm1 = zeros(size(xva));
urm = repmat(cam.u.xyz{n1},[np 1]);  cc1=[xva-A(1), yva-A(2), zm1]; cc1xys=(cc1(:,1).^2+cc1(:,2).^2)';
vrm = repmat(cam.u.xyz{n2},[np 1]);  cc2=[xva-B(1), yva-B(2), zm1]; cc2xys=(cc2(:,1).^2+cc2(:,2).^2)';
wrm = repmat(cam.u.xyz{n3},[np 1]);  cc3=[xva-C(1), yva-C(2), zm1]; cc3xys=(cc3(:,1).^2+cc3(:,2).^2)';

yaw1=repmat(fcnaz(cc1)',t);  yaw1=yaw1(:); cy1=cos(yaw1); sy1=sin(yaw1);
yaw2=repmat(fcnaz(cc2)',t);  yaw2=yaw2(:); cy2=cos(yaw2); sy2=sin(yaw2);
yaw3=repmat(fcnaz(cc3)',t);  yaw3=yaw3(:); cy3=cos(yaw3); sy3=sin(yaw3);

rollv = [0]*5*d2r;
nr = numel(rollv); %number of roll angles;
Uc = cell(nr,nz);
Vc = cell(nr,nz);
Wc = cell(nr,nz);
fxc = zeros(nr^3,1);
xc = cell(nr^3,1);
for i = 1:nr
    %fprintf('%.0f\n',i)
    roll1=rollv(i);  sr1=sin(roll1);  cr1=cos(roll1);
    for zi = 1:nz;
        z=zv(zi);
        sp1 = (A(3)-z)./sqrt(cc1xys + (z-A(3))^2);  cp1 = sqrt(1-sp1.*sp1); %sin pitch and cos pitch
        sp2 = (B(3)-z)./sqrt(cc2xys + (z-B(3))^2);  cp2 = sqrt(1-sp2.*sp2);
        sp3 = (C(3)-z)./sqrt(cc3xys + (z-C(3))^2);  cp3 = sqrt(1-sp3.*sp3);
        
        sp1=repmat(sp1,t);  cp1=repmat(cp1,t);  Uc{i,zi} = fcnrotateW2Bwsc(sr1,sp1(:),sy1,cr1,cp1(:),cy1,urm);
        sp2=repmat(sp2,t);  cp2=repmat(cp2,t);  Vc{i,zi} = fcnrotateW2Bwsc(sr1,sp2(:),sy2,cr1,cp2(:),cy2,vrm);
        sp3=repmat(sp3,t);  cp3=repmat(cp3,t);  Wc{i,zi} = fcnrotateW2Bwsc(sr1,sp3(:),sy3,cr1,cp3(:),cy3,wrm);
    end
end


n = 0;
for i = 1:nr
    cam.aposteriori.rpy(n1,1)=0;
    for j = 1:nr
        cam.aposteriori.rpy(n2,1)=0;
        for k = 1:nr
            cam.aposteriori.rpy(n3,1)=0;
            n = n+1;
            %fprintf('%3.0f',[i j k]); fprintf('--------------------------------------------------------------------------------------------------------------------------------------\n'); drawnow
            P = zeros(nx,nx,nz);
            for zi = 1:nz
                for yi = 1:nx
                    for xi = 1:nx
                        cam.tpnedhat(1,:) = [xv(xi) yv(yi) zv(zi)];
                        
                        rpy1 = solveAzPt_givenRoll(cam,a,n1);  %rpy1 = fcnB2WDCM2RPY(C); %compare to cam.true.rpy(vf,:)
                        u = fcnrotateW2B(rpy1(1),rpy1(2),rpy1(3),cam.uc.xyz{n1});
                        
                        %[~, D] = solveAzPt_givenRoll_fast(cam,a,n1);  %rpy1 = fcnB2WDCM2RPY(C); %compare to cam.true.rpy(vf,:)
                        %u2 = cam.uc.xyz{n1}*D';
                        
                        rpy1 = solveAzPt_givenRoll(cam,a,n2);  %rpy1 = fcnB2WDCM2RPY(C); %compare to cam.true.rpy(vf,:)
                        v = fcnrotateW2B(rpy1(1),rpy1(2),rpy1(3),cam.uc.xyz{n2});
                        
                        %[~, D] = solveAzPt_givenRoll_fast(cam,a,n2);  %rpy1 = fcnB2WDCM2RPY(C); %compare to cam.true.rpy(vf,:)
                        %v2 = cam.uc.xyz{n2}*D';
                        
                        rpy1 = solveAzPt_givenRoll(cam,a,n3);  %rpy1 = fcnB2WDCM2RPY(C); %compare to cam.true.rpy(vf,:)
                        w = fcnrotateW2B(rpy1(1),rpy1(2),rpy1(3),cam.uc.xyz{n3});
                        
                        %[~, D] = solveAzPt_givenRoll_fast(cam,a,n3);  %rpy1 = fcnB2WDCM2RPY(C); %compare to cam.true.rpy(vf,:)
                        %w2 = cam.uc.xyz{n3}*D';
                        
                        fx1 = fcnc03(A, u, B, v, C, w);
                        P(yi,xi,zi) = sum(fx1);
                    end
                end
                fprintf('%.0f,',zi)
            end
            [fxc(n),mi]= min3(P);  xc{n} = [xv(mi(2)) yv(mi(1)) zv(mi(3)) rollv([i j k])]; %ML focus point (x-y swVERIFIED)
        end
    end
end
[~, i]=min(fxc);  mx=xc{i};  c0=mx(1:3);  fprintf('\n')


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


vf = unique(round( linspace(n1,n3,min(30,n3)) ));
nf=numel(vf);  ova=ones(1,nf); ovb=ones(1,ntp); cam.msv.vf=vf;
cam.msv.ux=cam.u.x(vf,:); cam.msv.uy=cam.u.y(vf,:); cam.msv.uz=cam.u.z(vf,:);
D = cam.aposteriori.ned(vf,:);  cam.msv.A=D;  cam.msv.Ax=D(:,1);  cam.msv.Ay=D(:,2);  cam.msv.Az=D(:,3);
%DEFINE OPTIMIZER CONSTANTS
cam.msv.nf=nf; nx=cam.msv.nf*3;%number of frames
cam.msv.ntp=numel(cam.u.x(1,:)); %number of tie points
cam.msv.ova = ova;
cam.msv.ovb = ovb; 
%DEFINE PERMUTATIONS
x = 1:nf; 
j = tril(x(ova,:), -1); j=j(j~=0); cam.msv.j=uint16(j);
k = tril(x(ova,:)',-1); k=k(k~=0); cam.msv.k=uint16(k); cam.msv.njk=numel(j);
%DEFINE VECTOR ORIGINS
BAx=D(j,1)-D(k,1);  BAy=D(j,2)-D(k,2);  BAz=D(j,3)-D(k,3);  cam.msv.BAx=BAx(:,ovb);  cam.msv.BAy=BAy(:,ovb);  cam.msv.BAz=BAz(:,ovb);
cam.msv.BAxj=repmat(cam.msv.BAx,[1 1 nx+1]);  cam.msv.BAyj=repmat(cam.msv.BAy,[1 1 nx+1]);  cam.msv.BAzj=repmat(cam.msv.BAz,[1 1 nx+1]);

% x0 = [zeros(nf,1) fcnelaz(c0(ones(nf,1),:)-D)];
% xhat1 = lsqnonlin(@(x) fcnminseprjacobian(cam,vf,x),x0,[],[],cam.optim.options4);
% cam.tpnedhat = fcnminsepnq(cam,vf,xhat1);  [rmse,error] = fcngeterrors(cam,a,vf);
% fprintf('%.1fm RMSE tie points (%.1fm std), ',rmse.xyz,mean(std(error.xyz)))

%x0 = [c0 zeros(1,nf)];
%xhat4 = lsqnonlin(@(x) fcnminsepr2(cam,a,x),x0,[],[],cam.optim.options3);

x0 = [c0 zeros(1,nf)];
xhat4 = lsqnonlin(@(x) fcnminsepangle(cam,a,x),x0,[],[],cam.optim.options3);

%xhat1 = [xhat4(4:end)' fcnelaz(xhat4(ova,1:3)-cam.msv.A)]
rpy=zeros(cam.msv.nf,3);
cam.tpnedhat(1,:)=xhat4(1:3);
cam.aposteriori.rpy(cam.msv.vf,1)=xhat4(4:end)'*(180/pi);
for i=1:cam.msv.nf
    rpy(i,:) = solveAzPt_givenRoll(cam,a,cam.msv.vf(i));
end
xhat1 = rpy;
cam.tpnedhat = fcnminsepnq(cam,vf,xhat1);  [rmse,error] = fcngeterrors(cam,a,vf);
fprintf('%.1fm RMSE tie points (%.1fm std), ',rmse.xyz,mean(std(error.xyz)))



%RESECTION ----------------------------------------------------------------
zm = zeros(3*ntp,9);
for i=1:cam.frames
    z = fcnvec2uvec(cam.tpnedhat - cam.aposteriori.ned(i*ovb,:))';  %z = fcnvec2uvec(a.ipned - cam.true.ned(i*ovb,:))';
    u = cam.uc.xyz{i};  H=zm;  H(1:3:ntp*3,1:3)=u;  H(2:3:ntp*3,4:6)=u;  H(3:3:ntp*3,7:9)=u;    %u = fcnvec2uvec(fcnSC2CC(ned2camsc(cam,i,a.ipned,'true')*d2r));  H=zm;  H(1:3:ntp*3,1:3)=u;  H(2:3:ntp*3,4:6)=u;  H(3:3:ntp*3,7:9)=u;
    xhat = (H'*H)\H'*z(:); %LLS fast simple for equal measurement noise zh = H*xhat;
    DCM = reshape(xhat,[3 3]);  %C./sqrt(fcnrange(C)*fcnrange(C')')
    cam.aposteriori.rpy(i,:) = fcnW2BDCM2RPY(DCM)*r2d; %compare to cam.true.rpy(j,:)
end
rmse = fcngeterrors(cam,a,1:cam.frames);
fprintf('%.3fdeg after MSV resection\n',rmse.rpy)

% %NLS RESECTION -------------------------------------------------------------
% for i=1:cam.frames
%     cam.aposteriori.rpy(i,:) = fcnnlsresection(cam,a,i,cam.aposteriori.rpy(i,:))*r2d;
% end
% rmse = fcngeterrors(cam,a,1:cam.frames);
% fprintf(' (%.3fdeg after MSV NLS resection)\n',rmse.rpy)

% %PLOT ---------------------------------------------------------------------
% h=fig; %if j==1; h=fig;  axes(h(1)); end; 
% aned = cam.apriori.ned(1:cam.frames,:);
% fcnplot3(aned,'.b'); hold on;  fcnplot3(cam.tpnedhat(1:ntp,:),'.r'); fcnplot3(a.ipned(1:ntp,:),'.g');  plot3(aned([n1 n2 n3],1),aned([n1 n2 n3],2),aned([n1 n2 n3],3),'.r','markersize',15);
% text(aned(n1,1),aned(n1,2),aned(n1,3),sprintf(' %.0f',n1),'fontsize',8);  text(aned(n2,1),aned(n2,2),aned(n2,3),sprintf(' %.0f',n2),'fontsize',8);  text(aned(n3,1),aned(n3,2),aned(n3,3),sprintf(' %.0f',n3),'fontsize',8)
% axis tight; fcn3label; view(-60,30);  title(sprintf('HANK\n MinVal=%.3g, MinLocation=[%.2f %.2f %.2f]km',min3(P),c0(1),c0(2),c0(3)));
% n=1; ov=ones(n,1);
% c = ov*cam.apriori.ned(n1,:); d=ov*cam.tpnedhat(tpi,:);  plot3([c(:,1) d(:,1)]',[c(:,2) d(:,2)]',[c(:,3) d(:,3)]','r-')
% c = ov*cam.apriori.ned(n2,:); d=ov*cam.tpnedhat(tpi,:);  plot3([c(:,1) d(:,1)]',[c(:,2) d(:,2)]',[c(:,3) d(:,3)]','r-')
% c = ov*cam.apriori.ned(n3,:); d=ov*cam.tpnedhat(tpi,:);  plot3([c(:,1) d(:,1)]',[c(:,2) d(:,2)]',[c(:,3) d(:,3)]','r-')
% %plot3(cam.true.focus.ned(:,1), cam.true.focus.ned(:,2), cam.true.focus.ned(:,3), 'k*')
% fcnplot3D(xv,yv,zv,1/P,1-.8,cam.tpnedhat(tpi,3),a.ipned(tpi,3));
% 
% onesMat = ones(2);  zlim = fcnminmax(zv);  xlim = fcnminmax(xv);  ylim = fcnminmax(yv);
% x=cam.tpnedhat(tpi,1);  y=cam.tpnedhat(tpi,2);  z=cam.tpnedhat(tpi,3);
% plot3(xlim, [y y], [z z], '-r','LineWidth',1);
% plot3([x x], ylim, [z z], '-r','LineWidth',1);
% plot3([x x], [y y], zlim, '-r','LineWidth',1); plot3(x, y, z, '.r','MarkerSize',20);
% h = surf(xlim,ylim,onesMat*z,'FaceColor','r','EdgeColor',[1 0 0]','CData', onesMat ); alpha(h,0.1); %best estimate
% 
% x=a.ipned(tpi,1);  y=a.ipned(tpi,2);  z=a.ipned(tpi,3);
% plot3(xlim, [y y], [z z], '-g','LineWidth',1);
% plot3([x x], ylim, [z z], '-g','LineWidth',1);
% plot3([x x], [y y], zlim, '-g','LineWidth',1); plot3(x, y, z, '.g','MarkerSize',20);
% h = surf(xlim,ylim,onesMat*z,'FaceColor','g','EdgeColor',[0 1 0]','CData', onesMat ); alpha(h,0.1); %truth
% 
% legend('Apriori Aircraft Path','Estimates after Optimization','True Values','Location','Best')

rmse = rmse.each;
fprintf('\n')
end


function x = fcnnlsresection(cam,a,i,x)
fz = pixel2camsc(cam, [a.upx(:,i) a.upy(:,i)]); fz=fz(:,2:3)*d2r;
x=x(:)*d2r;
dx_ned = cam.tpnedhat-cam.apriori.ned(i*cam.msv.ovb,:);  TAx=dx_ned(:,1); TAy=dx_ned(:,2); TAz=dx_ned(:,3); 

%x = cam.apriori.rpy(i,:)';
for j=1:10
    r = x(1);  sr=sin(r); cr=cos(r);
    p = x(2);  sp=sin(p); cp=cos(p);
    y = x(3);  sy=sin(y); cy=cos(y);
    k1 = (TAx*(sr*sy + cr*cy*sp) - TAy*(cy*sr - cr*sp*sy) + TAz*cp*cr);
    k2 = (TAy*(cr*cy + sp*sr*sy) - TAx*(cr*sy - cy*sp*sr) + TAz*cp*sr);
    k3 = (TAx*cp*cy - TAz*sp + TAy*cp*sy);
    k4 = (k3.^2 + k1.^2 + k2.^2);
        
    C_NED2CAM = fcnRPY2DCM_B2W(x); %W2B but need transpose!
    sc = fcnCC2SCr(dx_ned*C_NED2CAM);  fhat1 = sc(:,2:3);
    
    dazdrpy = [ k1./((k2.^2./k3.^2 + 1).*k3), ((TAx*cp*cy*sr - TAz*sp*sr + TAy*cp*sr*sy)./k3 + (k2.*(TAz*cp + TAx*cy*sp + TAy*sp*sy))./k3.^2)./(k2.^2./k3.^2 + 1), -((TAx*(cr*cy + sp*sr*sy) + TAy*(cr*sy - cy*sp*sr))./k3 + ((TAy*cp*cy - TAx*cp*sy).*k2)./k3.^2)./(k2.^2./k3.^2 + 1)];
    deldrpy = [ k2./((1 - k1.^2./k4).^.5.*k4.^.5), -((TAx*cp*cr*cy - TAz*cr*sp + TAy*cp*cr*sy)./k4.^.5 - (k1.*(2*k2.*(TAx*cp*cy*sr - TAz*sp*sr + TAy*cp*sr*sy) - 2*k3.*(TAz*cp + TAx*cy*sp + TAy*sp*sy) + 2*k1.*(TAx*cp*cr*cy - TAz*cr*sp + TAy*cp*cr*sy)))./(2*k4.^(3./2)))./(1 - k1.^2./k4).^.5, -((TAx*(cy*sr - cr*sp*sy) + TAy*(sr*sy + cr*cy*sp))./k4.^.5 - (k1.*(2*(TAx*(cy*sr - cr*sp*sy) + TAy*(sr*sy + cr*cy*sp)).*k1 + 2*(TAy*cp*cy - TAx*cp*sy).*k3 - 2*(TAx*(cr*cy + sp*sr*sy) + TAy*(cr*sy - cy*sp*sr)).*k2))./(2*k4.^(3./2)))./(1 - k1.^2./k4).^.5];
    B = [deldrpy
        dazdrpy];
    
%     dx = 1E-8;
%     sc = fcnCC2SCr(dx_ned*fcnRPY2DCM_B2W(x+[0 0 0]'));   fhat = sc(:,2:3); f0=fhat(:);
%     sc = fcnCC2SCr(dx_ned*fcnRPY2DCM_B2W(x+[dx 0 0]'));  fhat = sc(:,2:3); f1=fhat(:);
%     sc = fcnCC2SCr(dx_ned*fcnRPY2DCM_B2W(x+[0 dx 0]'));  fhat = sc(:,2:3); f2=fhat(:);
%     sc = fcnCC2SCr(dx_ned*fcnRPY2DCM_B2W(x+[0 0 dx]'));  fhat = sc(:,2:3); f3=fhat(:);
%     B = [f1-f0 f2-f0 f3-f0]/dx;
    
    f = fz - fhat1;
    Bt = B';
    N = Bt*B;
    x = x+N\(Bt*f(:));
    %fcnrms(f)
end
end

function [rmse, error] = fcngeterrors(cam,a,i)
drp = cam.aposteriori.rpy(i,1:2) - cam.true.rpy(i,1:2);
daz = fcndaz(cam.aposteriori.rpy(i,3)*d2r,cam.true.rpy(i,3)*d2r)*r2d;
error.xyz = (cam.tpnedhat-a.ipned)*1000;
error.rpy = [drp daz];
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

function x2 = fcnrotateB2W(r,p,y,x)
sr=sin(r); sp=sin(p); sy=sin(y);
cr=cos(r); cp=cos(p); cy=cos(y);
x2 = [x(:,1).*(cp.*cy)+x(:,2).*(cp.*sy )+x(:,3).*(-sp),   x(:,1).*(sr.*sp.*cy-cr.*sy)+x(:,2).*(sr.*sp.*sy+cr.*cy)+x(:,3).*(sr.*cp),   x(:,1).*(cr.*sp.*cy+sr.*sy)+x(:,2).*(cr.*sp.*sy-sr.*cy)+x(:,3).*(cr.*cp)];
end

function rs = fcnminsepr(cam,rpy) %[rpyrpyrpy...]
%LOAD CONSTANTS
ovb = cam.msv.ovb;  vf = cam.msv.vf;
j = cam.msv.j;  k = cam.msv.k;
BAx=cam.msv.BAx; BAy=cam.msv.BAy; BAz=cam.msv.BAz;
x=cam.u.x(vf,:); y=cam.u.y(vf,:); z=cam.u.z(vf,:);

%VECTOR ROTATIONS
r = rpy(:,1); cr=cos(r); sr=sin(r); cr=cr(:,ovb); sr=sr(:,ovb);
p = rpy(:,2); cp=cos(p); sp=sin(p); cp=cp(:,ovb); sp=sp(:,ovb);
a = rpy(:,3); cy=cos(a); sy=sin(a); cy=cy(:,ovb); sy=sy(:,ovb); %az
ux1=x.*(cp.*cy) +y.*(sr.*sp.*cy-cr.*sy)  +z.*(cr.*sp.*cy+sr.*sy);  %[ux1,uy1,uz1] = fcnrotateW2Brpyxyz(sr,sp,sy,cr,cp,cy,x,y,z);
uy1=x.*(cp.*sy) +y.*(sr.*sp.*sy+cr.*cy)  +z.*(cr.*sp.*sy-sr.*cy);
uz1=x.*(-sp)    +y.*(sr.*cp)             +z.*(cr.*cp);
ux=ux1(j,:);  uy=uy1(j,:);  uz=uz1(j,:);
vx=ux1(k,:);  vy=uy1(k,:);  vz=uz1(k,:);

%VECTOR INTERCEPTS
d = ux.*vx  + uy.*vy  + uz.*vz;
e = ux.*BAx + uy.*BAy + uz.*BAz;
f = vx.*BAx + vy.*BAy + vz.*BAz;
g = 1 - d.*d;
s1 = (d.*f - e)./g; %multiply times u
t1 = (f - d.*e)./g; %multiply times v

%MISCLOSURE VECTOR RANGE RESIDUALS
rs =  sqrt( (t1.*vx-BAx-s1.*ux).^2 + (t1.*vy-BAy-s1.*uy).^2 + (t1.*vz-BAz-s1.*uz).^2 ); %sum of the squared ranges of the misclosure vectors
end

function [rs, B] = fcnminseprjacobian(cam,vf,rpy) %[rpyrpyrpy...]
%LOAD CONSTANTS
ovb = cam.msv.ovb; ntp=numel(ovb); nx=numel(rpy);
j = cam.msv.j;  k = cam.msv.k;  njk=numel(j);
BAx=cam.msv.BAxj; BAy=cam.msv.BAyj; BAz=cam.msv.BAzj;
x=cam.u.x(vf,:); y=cam.u.y(vf,:); z=cam.u.z(vf,:);
ux=zeros(njk,ntp,nx+1);  uy=ux;  uz=ux;  vx=ux;  vy=ux;  vz=ux;

%VECTOR ROTATIONS
dx = 1E-8;
vroll = 1:cam.msv.nf;  vpitch = cam.msv.nf+1:cam.msv.nf*2;  vyaw = cam.msv.nf*2+1:cam.msv.nf*3;
for i=1:nx+1
    rpy1=rpy(:);  
    if i>1
        rpy1(i-1)=rpy1(i-1)+dx;
    end

    r = rpy1(vroll);  cr=cos(r); sr=sin(r); cr=cr(:,ovb); sr=sr(:,ovb);
    p = rpy1(vpitch); cp=cos(p); sp=sin(p); cp=cp(:,ovb); sp=sp(:,ovb);
    a = rpy1(vyaw);   cy=cos(a); sy=sin(a); cy=cy(:,ovb); sy=sy(:,ovb); %az
    ux1=x.*(cp.*cy) +y.*(sr.*sp.*cy-cr.*sy)  +z.*(cr.*sp.*cy+sr.*sy);  %[ux1,uy1,uz1] = fcnrotateW2Brpyxyz(sr,sp,sy,cr,cp,cy,x,y,z);
    uy1=x.*(cp.*sy) +y.*(sr.*sp.*sy+cr.*cy)  +z.*(cr.*sp.*sy-sr.*cy);
    uz1=x.*(-sp)    +y.*(sr.*cp)             +z.*(cr.*cp);
    ux(:,:,i)=ux1(j,:);  uy(:,:,i)=uy1(j,:);  uz(:,:,i)=uz1(j,:);
    vx(:,:,i)=ux1(k,:);  vy(:,:,i)=uy1(k,:);  vz(:,:,i)=uz1(k,:);
end

%VECTOR INTERCEPTS
d = ux.*vx  + uy.*vy  + uz.*vz;
e = ux.*BAx + uy.*BAy + uz.*BAz;
f = vx.*BAx + vy.*BAy + vz.*BAz;
g = 1 - d.*d;
s1 = (d.*f - e)./g; %multiply times u
t1 = (f - d.*e)./g; %multiply times v

%MISCLOSURE VECTOR RANGE RESIDUALS
rs =  sqrt( (t1.*vx-BAx-s1.*ux).^2 + (t1.*vy-BAy-s1.*uy).^2 + (t1.*vz-BAz-s1.*uz).^2 ); %sum of the squared ranges of the misclosure vectors

drs = rs(:,:,2:end)-rs(:,:,ones(nx,1));  B = reshape(drs,[njk*ntp nx])./dx;
rs = rs(:,:,1);
end


function rs = fcnminsepr2(cam,a,xhat) %[xyzrrr...]
ovb = cam.msv.ovb;
j = cam.msv.j;  k = cam.msv.k;
BAx=cam.msv.BAx; BAy=cam.msv.BAy; BAz=cam.msv.BAz;
r = xhat(4:end)';

% %ORIGINAL VECTOR ROTATIONS
%xt = x(1)-cam.msv.Ax;  yt = x(2)-cam.msv.Ay;  zt = x(3)-cam.msv.Az;  yaw = 2*atan((sqrt(xt.^2+yt.^2)-xt)./yt);
%sp = -zt./sqrt(xt.^2 + yt.^2 + zt.^2);  cp = sqrt(1 - sp.^2);  %p = asin(-zt./sqrt(xt.^2 + yt.^2 + zt.^2));
% cr=cos(r); sr=sin(r);       cr=cr(:,ovb); sr=sr(:,ovb);
%                             cp=cp(:,ovb); sp=sp(:,ovb);
% cy=cos(yaw); sy=sin(yaw);   cy=cy(:,ovb); sy=sy(:,ovb);
% x=cam.msv.ux; y=cam.msv.uy; z=cam.msv.uz;

%HANK STUFF ---------------------------------------------------------------
vf=cam.msv.vf; x=cam.uc.x(vf,:); y=cam.uc.y(vf,:); z=cam.uc.z(vf,:);
rpy=zeros(cam.msv.nf,3);
cam.tpnedhat(1,:)=xhat(1:3);
cam.aposteriori.rpy(vf,1)=r*(180/pi);
for i=1:cam.msv.nf
    %[~,~,Mned] = solveAzPt_givenRoll(cam,a,cam.msv.vf(i));
    %rpy(i,:) = fcnB2WDCM2RPY(Mned); %compare to cam.true.rpy(vf,:) 
    rpy(i,:) = solveAzPt_givenRoll(cam,a,cam.msv.vf(i));
end
r=rpy(:,1);     cr=cos(r); sr=sin(r);       cr=cr(:,ovb); sr=sr(:,ovb);
p=rpy(:,2);     cp=cos(p); sp=sin(p);       cp=cp(:,ovb); sp=sp(:,ovb);
yaw=rpy(:,3);   cy=cos(yaw); sy=sin(yaw);   cy=cy(:,ovb); sy=sy(:,ovb);
%HANK STUFF ---------------------------------------------------------------

%VECTOR ROTATIONS
ux1=x.*(cp.*cy) +y.*(sr.*sp.*cy-cr.*sy)  +z.*(cr.*sp.*cy+sr.*sy);  %[ux1,uy1,uz1] = fcnrotateW2Brpyxyz(sr,sp,sy,cr,cp,cy,x,y,z);
uy1=x.*(cp.*sy) +y.*(sr.*sp.*sy+cr.*cy)  +z.*(cr.*sp.*sy-sr.*cy);
uz1=x.*(-sp)    +y.*(sr.*cp)             +z.*(cr.*cp);
vx=ux1(k,:);  vy=uy1(k,:);  vz=uz1(k,:);
ux=ux1(j,:);  uy=uy1(j,:);  uz=uz1(j,:);

%VECTOR INTERCEPTS
d = ux.*vx + uy.*vy + uz.*vz;
e = ux.*BAx + uy.*BAy + uz.*BAz;
f = vx.*BAx + vy.*BAy + vz.*BAz;
g = 1 - d.*d;
s1 = (d.*f - e)./g; %multiply times u
t1 = (f - d.*e)./g; %multiply times v

%MISCLOSURE VECTOR RANGE RESIDUALS
rs =  sqrt( (t1.*vx-BAx-s1.*ux).^2 + (t1.*vy-BAy-s1.*uy).^2 + (t1.*vz-BAz-s1.*uz).^2 ); %sum of the squared ranges of the misclosure vectors
end

function ts1 = fcnminsepangle(cam,a,xhat)
%LOAD CONSTANTS
ova = cam.msv.ova;  ovb = cam.msv.ovb;
j = cam.msv.j;  k = cam.msv.k;
A = cam.msv.A;  BAx=cam.msv.BAx; BAy=cam.msv.BAy; BAz=cam.msv.BAz;
r = xhat(4:end)';

% %ORIGINAL VECTOR ROTATIONS
% xt = x(1)-cam.msv.Ax;  yt = x(2)-cam.msv.Ay;  zt = x(3)-cam.msv.Az;  y = 2*atan((sqrt(xt.^2+yt.^2)-xt)./yt);
% sp = -zt./sqrt(xt.^2 + yt.^2 + zt.^2);  cp = sqrt(1 - sp.^2);  %p = asin(-zt./sqrt(xt.^2 + yt.^2 + zt.^2));
% cr=cos(r); sr=sin(r); cr=cr(:,ovb); sr=sr(:,ovb);
%                       cp=cp(:,ovb); sp=sp(:,ovb);
% cy=cos(y); sy=sin(y); cy=cy(:,ovb); sy=sy(:,ovb);
% x=cam.msv.ux; y=cam.msv.uy; z=cam.msv.uz;

%HANK STUFF ---------------------------------------------------------------
vf=cam.msv.vf; x=cam.uc.x(vf,:); y=cam.uc.y(vf,:); z=cam.uc.z(vf,:);
rpy=zeros(cam.msv.nf,3);
cam.tpnedhat(1,:)=xhat(1:3);
cam.aposteriori.rpy(vf,1)=r*(180/pi);
for i=1:cam.msv.nf
    %[~,~,Mned] = solveAzPt_givenRoll(cam,a,cam.msv.vf(i));
    %rpy(i,:) = fcnB2WDCM2RPY(Mned); %compare to cam.true.rpy(vf,:) 
    rpy(i,:) = solveAzPt_givenRoll(cam,a,cam.msv.vf(i));
end
r=rpy(:,1);     cr=cos(r); sr=sin(r);       cr=cr(:,ovb); sr=sr(:,ovb);
p=rpy(:,2);     cp=cos(p); sp=sin(p);       cp=cp(:,ovb); sp=sp(:,ovb);
yaw=rpy(:,3);   cy=cos(yaw); sy=sin(yaw);   cy=cy(:,ovb); sy=sy(:,ovb);
%HANK STUFF ---------------------------------------------------------------

ux1=x.*(cp.*cy) +y.*(sr.*sp.*cy-cr.*sy)  +z.*(cr.*sp.*cy+sr.*sy);  %[ux1,uy1,uz1] = fcnrotateW2Brpyxyz(sr,sp,sy,cr,cp,cy,x,y,z);
uy1=x.*(cp.*sy) +y.*(sr.*sp.*sy+cr.*cy)  +z.*(cr.*sp.*sy-sr.*cy);
uz1=x.*(-sp)    +y.*(sr.*cp)             +z.*(cr.*cp);
vx=ux1(k,:);  vy=uy1(k,:);  vz=uz1(k,:);
ux=ux1(j,:);  uy=uy1(j,:);  uz=uz1(j,:);

%VECTOR INTERCEPTS
d = ux.*vx + uy.*vy + uz.*vz;
e = ux.*BAx + uy.*BAy + uz.*BAz;
f = vx.*BAx + vy.*BAy + vz.*BAz;
g = 1 - d.*d;
s1 = (d.*f - e)./g; %multiply times u
t1 = (f - d.*e)./g; %multiply times v

%MISCLOSURE VECTOR RANGE RESIDUALS
%rs =  (t1.*vx-BAx-s1.*ux).^2 + (t1.*vy-BAy-s1.*uy).^2 + (t1.*vz-BAz-s1.*uz).^2 ; %sum of the squared ranges of the misclosure vectors

%TIE POINT CENTERS
den = cam.msv.njk*2; %denominator = number of permutations times 2
D = sum(A)*(cam.msv.nf-1);
C0x = (sum(t1.*vx+s1.*ux)+D(1)) / den;
C0y = (sum(t1.*vy+s1.*uy)+D(2)) / den;
C0z = (sum(t1.*vz+s1.*uz)+D(3)) / den;
%C0 = zeros(ntp,3); C0(:,1)=C0x; C0(:,2)=C0y; C0(:,3)=C0z;

%IMAGE ANGLE RESIDUALS
CAx = C0x(ova,:)-A(:,ovb);  CAy = C0y(ova,:)-A(:,ovb*2);  CAz = C0z(ova,:)-A(:,ovb*3);
%ts1 = nf*ntp - sum(sum(  (ux1.*CAx + uy1.*CAy + uz1.*CAz).^2./(CAx.*CAx+CAy.*CAy+CAz.*CAz)  ));  %sum of sqrt of angle sines  %sin = sqrt(1-ct^2)
ts1 = 1 - (ux1.*CAx + uy1.*CAy + uz1.*CAz).^2./(CAx.*CAx+CAy.*CAy+CAz.*CAz) ;  %sum of sqrt of angle sines  %sin = sqrt(1-ct^2)
end

function [C0, rs, ts1] = fcnminsepnq(cam,vf,rpy)
%LOAD CONSTANTS
nf = cam.msv.nf;  ntp = cam.msv.ntp;
ova = cam.msv.ova;  ovb = cam.msv.ovb;
j = cam.msv.j;  k = cam.msv.k;
A = cam.msv.A;  BAx=cam.msv.BAx; BAy=cam.msv.BAy; BAz=cam.msv.BAz;

%VECTOR ROTATIONS
r = rpy(:,1); cr=cos(r); sr=sin(r); cr=cr(:,ovb); sr=sr(:,ovb);
p = rpy(:,2); cp=cos(p); sp=sin(p); cp=cp(:,ovb); sp=sp(:,ovb);
y = rpy(:,3); cy=cos(y); sy=sin(y); cy=cy(:,ovb); sy=sy(:,ovb);
x=cam.uc.x(vf,:); y=cam.uc.y(vf,:); z=cam.uc.z(vf,:);
ux1=x.*(cp.*cy) +y.*(sr.*sp.*cy-cr.*sy)  +z.*(cr.*sp.*cy+sr.*sy);  %[ux1,uy1,uz1] = fcnrotateW2Brpyxyz(sr,sp,sy,cr,cp,cy,x,y,z);
uy1=x.*(cp.*sy) +y.*(sr.*sp.*sy+cr.*cy)  +z.*(cr.*sp.*sy-sr.*cy);
uz1=x.*(-sp)    +y.*(sr.*cp)             +z.*(cr.*cp);
ux=ux1(j,:);  uy=uy1(j,:);  uz=uz1(j,:);
vx=ux1(k,:);  vy=uy1(k,:);  vz=uz1(k,:);

%VECTOR INTERCEPTS
d = ux.*vx + uy.*vy + uz.*vz;
e = ux.*BAx + uy.*BAy + uz.*BAz;
f = vx.*BAx + vy.*BAy + vz.*BAz;
g = 1 - d.*d;
s1 = (d.*f - e)./g; %multiply times u
t1 = (f - d.*e)./g; %multiply times v

%MISCLOSURE VECTOR RANGE RESIDUALS
rs =  (t1.*vx-BAx-s1.*ux).^2 + (t1.*vy-BAy-s1.*uy).^2 + (t1.*vz-BAz-s1.*uz).^2 ; %sum of the squared ranges of the misclosure vectors

%TIE POINT CENTERS
den = numel(j)*2; %denominator = number of permutations times 2
D = sum(A)*(nf-1);
C0x = (sum(t1.*vx+s1.*ux)+D(1)) / den;
C0y = (sum(t1.*vy+s1.*uy)+D(2)) / den;
C0z = (sum(t1.*vz+s1.*uz)+D(3)) / den;
C0 = zeros(ntp,3); C0(:,1)=C0x; C0(:,2)=C0y; C0(:,3)=C0z;

%IMAGE ANGLE RESIDUALS
CAx = C0x(ova,:)-A(:,ovb);  CAy = C0y(ova,:)-A(:,ovb*2);  CAz = C0z(ova,:)-A(:,ovb*3);
%ts1 = nf*ntp - sum(sum(  (ux1.*CAx + uy1.*CAy + uz1.*CAz).^2./(CAx.*CAx+CAy.*CAy+CAz.*CAz)  ));  %sum of sqrt of angle sines  %sin = sqrt(1-ct^2)
ts1 = 1 - (ux1.*CAx + uy1.*CAy + uz1.*CAz).^2./(CAx.*CAx+CAy.*CAy+CAz.*CAz) ;  %sum of sqrt of angle sines  %sin = sqrt(1-ct^2)
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

ts1 = (u(:,1).*(C0x-A(1)) + u(:,2).*(C0y-A(2)) + u(:,3).*(C0z-A(3))).^2./((C0x-A(1)).^2+(C0y-A(2)).^2+(C0z-A(3)).^2);  %sin = sqrt(1-ct^2)
ts2 = (v(:,1).*(C0x-B(1)) + v(:,2).*(C0y-B(2)) + v(:,3).*(C0z-B(3))).^2./((C0x-B(1)).^2+(C0y-B(2)).^2+(C0z-B(3)).^2);
ts3 = (w(:,1).*(C0x-C(1)) + w(:,2).*(C0y-C(2)) + w(:,3).*(C0z-C(3))).^2./((C0x-C(1)).^2+(C0y-C(2)).^2+(C0z-C(3)).^2);  rs = 3-(ts1+ts2+ts3); %sum of sqrt of angle sines  %ts4 = sqrt(ts1)+sqrt(ts2)+sqrt(ts3); %sum of angle sines
end

function fcnplot3D(xv,yv,zv,v,percentile,redv,greenv)
isovalue          = fcndsearch(v,percentile);
[ml.val, i]       = max3(v);  ml.row=i(1);  ml.col=i(2);  ml.layer=i(3);
onesMat = ones(2);  zlim = fcnminmax(zv);  xlim = fcnminmax(xv);  ylim = fcnminmax(yv);

%plot mlpoint
x=xv(ml.col);  y=yv(ml.row);  z=zv(ml.layer);
plot3(xlim, [y y], [z z], '-k','LineWidth',1,'color',[.7 .7 .7]); hold on
plot3([x x], ylim, [z z], '-k','LineWidth',1,'color',[.7 .7 .7]);
plot3([x x], [y y], zlim, '-k','LineWidth',1,'color',[.7 .7 .7]); plot3(x, y, z, '.','MarkerSize',20,'Color',[.7 .7 .7]);
h = surf(xlim,ylim,onesMat*z,'FaceColor',[.4 .4 .4],'EdgeColor',[.7 .7 .7]','CData', onesMat ); alpha(h,0.1);

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

hcont = contourslice(xm,ym,zm,v,[],[],redv, [isovalue isovalue]);
set(hcont,'EdgeColor',[1 0 0],'LineWidth',1)

hcont = contourslice(xm,ym,zm,v,[],[],greenv, [isovalue isovalue]);
set(hcont,'EdgeColor',[0 1 0],'LineWidth',1)


fcn3label; fcnfontsize; box on; axis equal vis3d tight
end


