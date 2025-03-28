% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [cam, h] = fcnSFM(cam,a,varargin)
%cam.optim.options3 = optimset('Algorithm','levenberg-marquardt','display','iter','TolX',1E-16,'Tolfun',1E-16,'MaxFunEvals',9000,'MaxIter',3000);
sequentialflag = false;  sequentialstr = 'non-sequential'; %default is non-sequential SFM
method='Horn'; %default is Horn method in fcnE2P
if nargin==3
    str = varargin{1};
    if ~isempty(regexpi(str,'zisserman'));  method='Zisserman';  end
    if isempty(regexpi(str,'non')) && ~isempty(regexpi(str,'sequential')); sequentialflag=true; sequentialstr='sequential';  end
end
startclock = clock; fprintf('SFM... ')

Ccn = [0 0 1; 1 0 0; 0 -1 0];% DCM cam (z x -y) to ned (x y z)
ntp = size(a.upx,1); %number of tie points
K = fcnK(cam.width,cam.height,cam.xppo,cam.yppo,cam.focalLength); %camera calibration matrix
nf = cam.frames;  vf0 = cam.msv.vf0;
af = 1; %special keyframe
if sequentialflag; af=1; end%sequential flag. 1=sequential SFM, 0=keyframe SFM

afi = af;
Ci2af = cell(nf,1);
zv=zeros(nf,ntp); cam.u.x=zv; cam.u.y=zv; cam.u.z=zv; cam.uc.x=zv; cam.uc.y=zv; cam.uc.z=zv;
for i=1:nf
    fprintf('%.0f/%.0f\n',i,nf)
    f2 = [a.upx(:,i)   a.upy(:,i)]; %frame 2
    
    if i>1 || sequentialflag==0
        if sequentialflag;  afi = i-1;  end %#ok<*UNRCH>
        f1 = [a.upx(:,afi)  a.upy(:,afi)]; %frame 1 (af)
        j = fcnfindvalidtp(a,[afi i]);  ov = ones(sum(j),1);
        x1s = [f1(j,:) ov]'; %x1squiggle
        x2s = [f2(j,:) ov]';
        
        F = fcn8ptF(x1s, x2s);
        E = fcnF2E(F,K);
        P = fcnE2P(E,method);
        P = getCorrectCameraMatrix(P, K, K,x1s,x2s);
        R = P(:,1:3); %Ri2af
        %t = P(:,4); %t = cam.true.ned(f,:)';
        Ci2h = Ccn*R'*Ccn'; %C frame i to frame af 
        if sequentialflag
           Ci2af{i} = Ci2af{i-1} * Ci2h;    
        else
           Ci2af{i} = Ci2h; 
        end
    end
    if i==af
        Ci2af{i} = eye(3);
    end

    %GET CAMERA UNIT VECTORS ----------------------------------------------
    ui = pixel2camcc(cam, f2);  cam.uc.x(i,:)=ui(:,1);  cam.uc.y(i,:)=ui(:,2);  cam.uc.z(i,:)=ui(:,3);
    ui = ui * Ci2af{i}';        cam.u.x(i,:) =ui(:,1);  cam.u.y(i,:) =ui(:,2);  cam.u.z(i,:) =ui(:,3);
end


%OPTIMIZER CONSTANTS ------------------------------------------------------
vf = unique(round( linspace(1,nf,min(30,nf)) )); %at most use 100 frames
vtp = fcnfindvalidtp(a,vf);  b = fcncropTP(cam,a, vtp);  ntp = size(b.upx,1); %number of tie points
nf=numel(vf);  ova=ones(1,nf);  ovb=ones(1,ntp);  cam.msv.vf=vf;  cam.msv.nf=nf;  cam.msv.ova=ova;  cam.msv.ovb=ovb;  cam.msv.ntp=ntp; 
cam.msv.ux=cam.uc.x(vf,vtp);  cam.msv.uy=cam.uc.y(vf,vtp);  cam.msv.uz=cam.uc.z(vf,vtp);  cam.msv.A=cam.aposteriori.ned(vf,:);

% %OPTIMIZE -----------------------------------------------------------------
h=fig(1,2,.9); axes(h(1));  %h=fig(2,1,1.98); axes(h(1)); 

% ni=50*8; nj=25*8; f1=zeros(ni,nj);
% azvec = linspace(-pi,pi,ni); %azvec = linspace(-pi,pi-2*pi/ni,ni);
% elvec = linspace(-pi/2,pi/2,nj);
% 
% % rpy = zeros(3,nf);
% % Caf2w = fcnRPY2DCM_B2W(cam.true.rpy(af,:)*d2r);
% % for i = 1:numel(vf);
% %     rpy(:,i) = fcnB2WDCM2RPY( (Caf2w * Ci2af{vf(i)}) );
% % end
% % fcnminrpyrpyrpy(cam,b, rpy(:)')
% 
% zm = zeros(3,nf);
% parfor i=1:ni
%     for j=1:nj
%         rpy = zm;
%         Caf2w = fcnRPY2DCM_B2W([0 elvec(j) azvec(i)]);
%         for k = 1:numel(vf);
%             rpy(:,k) = fcnB2WDCM2RPY( Caf2w * Ci2af{vf(k)} );
%         end
%         
%         f1(i,j) = sum3( fcnminrpyrpyrpy(cam,b,rpy(:)').^2 );
%     end; fprintf('%.0f/%.0f\n',i,ni)
% end
% %[~,j] = min3(f1); x0=[0 elvec(j(2)) azvec(j(1))];  %x0 = lsqnonlin(@(x) fcnminrpyrpyrpy(cam,a,x),x0,[],[],cam.optim.options3); 
% pcolor(azvec*r2d,elvec*r2d,-log(f1')); shading flat; axis equal tight; xlabel('keyframe yaw (deg)'); ylabel('keyframe pitch (deg)'); 
% title(sprintf('{\\bf%s %s SFM} keyframe search',method,sequentialstr))
% x0 = cam.true.rpy(af,:);  plot3(x0(:,3),x0(:,2),0,'g+','linewidth',1.5,'markersize',30,'color',[.8 .8 .8]); colorbar
% set(gca,'ytick',-90:45:90); set(gca,'xtick',-180:45:180)

if cam.syntheticVideoFlag
    x0 = cam.true.rpy(af,:);
else
    x0 = cam.aposteriori.rpy(af,:);  
end
Caf2w = fcnRPY2DCM_B2W(x0*d2r);
for i=1:cam.frames
    Ci2w = Caf2w * Ci2af{i};
    cam.aposteriori.rpySFM(i,:) = fcnB2WDCM2RPY(Ci2w)*r2d;
end
cam.tpnedhat = fcnMIGMSV(cam,a);  %fcnrms(a.ipned - cam.tpnedhat)*1000


%PLOT ---------------------------------------------------------------------
f = 1:cam.frames;
ned = cam.aposteriori.ned(f,:);
ov = ones(cam.frames,1)*1;
vt = fcnSC2CCd([ov cam.true.rpy(f,2:3)]);
g = zeros(5,1);
if cam.syntheticVideoFlag
    nedFocus = cam.true.focus.ned;
    nedTP = a.ipned;
else
    cam.aposteriori.focus.lla = zeros(cam.frames,3);
    for i=1:cam.frames
        cam.aposteriori.focus.lla(i,:) = pixel2lla(cam.DEM,cam,i,cam.true.focus.pixel,'aposteriori');
    end
    nedFocus = lla2ned(cam.DEM,cam.aposteriori.focus.lla);
    nedTP = cam.tpnedhat;
end
r = mean(fcnrange(nedFocus,cam.apriori.ned))*2;


axes(h(2))
v = fcnSC2CCd([ov cam.aposteriori.rpy(f,2:3)])*r;
x1 = ned(:,1); x2 = x1 + v(:,1)*10;
y1 = ned(:,2); y2 = y1 + v(:,2)*10;
z1 = ned(:,3); z2 = z1 + v(:,3)*10;
g(1)=quiver3(x1,y1,z1,v(:,1),v(:,2),v(:,3),0,'b');
g(2)=quiver3(x1(vf0),y1(vf0),z1(vf0),v(vf0,1),v(vf0,2),v(vf0,3),0,'linewidth',2,'color',[1 .6 0]); %keyframe
box on; axis equal vis3d; fcn3label(' (km)'); set(gca,'zdir','Reverse'); view(-90,90)
g(3)=fcnplot3(ned,'b.');  g(4)=fcnplot3(nedTP,'g.'); g(5)=fcnplot3(nedFocus(f,:),'r.'); fcnmarkersize(gca,6)
legend(g,'LOS','keyframe LOS','cameras C','tie points X','focus points F','location','SouthWest');
if cam.syntheticVideoFlag
    title(sprintf('{\\bfMSVM} solution (%.1f^o RMSE)',fcnrms(fcnangle(vt,v)*r2d)))
else
    title(sprintf('{\\bfMSVM} solution (unknown RMSE)'))
end


axes(h(1)) %axes(h(3))
v = fcnSC2CCd([ov cam.aposteriori.rpySFM(f,2:3)])*r;
x1 = ned(:,1); x2 = x1 + v(:,1)*10;
y1 = ned(:,2); y2 = y1 + v(:,2)*10;
z1 = ned(:,3); z2 = z1 + v(:,3)*10;
g(1)=quiver3(x1,y1,z1,v(:,1),v(:,2),v(:,3),0,'b');
g(2)=quiver3(x1(af),y1(af),z1(af),v(af,1),v(af,2),v(af,3),0,'linewidth',2,'color',[1 .6 0]); %keyframe
box on; axis equal vis3d; fcn3label(' (km)'); set(gca,'zdir','Reverse'); view(-90,90)
g(3)=fcnplot3(ned,'b.');  g(4)=fcnplot3(nedTP,'g.'); g(5)=fcnplot3(nedFocus(f,:),'r.'); fcnmarkersize(gca,6)
legend(g,'LOS','keyframe LOS','cameras C','tie points X','focus points F','location','SouthWest');
if cam.syntheticVideoFlag
    title(sprintf('{\\bf%s %s SFM} (%.1f^o RMSE)',method,sequentialstr,fcnrms(fcnangle(vt,v)*r2d)))
else
    title(sprintf('{\\bf%s %s SFM} (unknown RMSE)',method,sequentialstr))
end
fcnmatchlims(h(2),h(1))%fcnmatchlims(h(2),h(3))


% str = [];
% lla1 = ecef2lla(ned2ecef(cam.DEM,[x1 y1 z1])); 
% lla2 = ecef2lla(ned2ecef(cam.DEM,[x2 y2 z2]));
% if ~cam.syntheticVideoFlag; lla1=lla2llag(lla1); lla2=lla2llag(lla2); end
% for i=f
%     if i==af
%         color = [1 .7 0];
%         linewidth = 4;
%     else
%         color = [1 0 0];
%         linewidth = 2;
%     end
%     str1 = ge_plot3([lla1(i,2) lla2(i,2)]',[lla1(i,1) lla2(i,1)]',[lla1(i,3) lla2(i,3)]','extrude',0,'lineWidth',linewidth,'altitudeMode','absolute','lineColor',['FF' fcnhexcolor(color)]);
%     str = [str str1]; %#ok<*AGROW>
% end
% cam.google.kml.SFMfocuspoint = ge_folder(sprintf('%s %s SFM LOS',method,sequentialstr),str);

%[cam] = fcngetcamcornersLLA(cam);
fprintf(' Done. (%.0fs)\n',etime(clock,startclock))
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
%x0=cam.uc.x(vf,:); y0=cam.uc.y(vf,:); z0=cam.uc.z(vf,:);  
x0=cam.msv.ux; y0=cam.msv.uy; z0=cam.msv.uz;
ux1=x0.*(cp.*cy) +y0.*(sr.*sp.*cy-cr.*sy)  +z0.*(cr.*sp.*cy+sr.*sy);  %[ux1,uy1,uz1] = fcnrotateW2Brpyxyz(sr,sp,sy,cr,cp,cy,x,y,z);
uy1=x0.*(cp.*sy) +y0.*(sr.*sp.*sy+cr.*cy)  +z0.*(cr.*sp.*sy-sr.*cy);
uz1=x0.*(-sp)    +y0.*(sr.*cp)             +z0.*(cr.*cp);


% fig;
% Ax = repmat(A(:,1),[1,ntp]);  Ay = repmat(A(:,2),[1,ntp]);  Az = repmat(A(:,3),[1,ntp]);
% quiver3(Ax(:),Ay(:),Az(:),ux1(:),uy1(:),uz1(:),100,'b');  box on; axis equal vis3d; fcn3label; set(gca,'zdir','Reverse'); 

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
