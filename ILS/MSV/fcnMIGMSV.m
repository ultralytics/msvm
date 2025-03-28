% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [C0, msvr] = fcnMIGMSV(cam,a)
ntp = size(a.upx,1);

cam.u.xyz=cell(cam.frames,1); cam.uc.xyz=cam.u.xyz; zv=zeros(cam.frames,ntp); cam.uc.x=zv; cam.uc.y=zv; cam.uc.z=zv;for i=1:cam.frames
    uic = pixel2camcc(cam, [a.upx(:,i) a.upy(:,i)]);    %uic=fcnSC2CCd(pixel2camsc(cam, [a.upx(:,i) a.upy(:,i)]));
    cam.uc.x(i,:) = uic(:,1);  cam.uc.y(i,:) = uic(:,2);  cam.uc.z(i,:) = uic(:,3);  cam.uc.xyz{i} = uic;
end

ntp1 = 1;
C0 = zeros(ntp,3);
msvr = zeros(ntp,1);
for i=1:ntp
    vf = find(a.state(i,:)==1);  nf=numel(vf);
    %[n1, n2, n3] = fcnget3distantframes(cam.aposteriori.ned(vf,:));  vf = vf([n1 n2 n3]);
    if nf>5;  vf = vf(unique(round(linspace(1,nf,5))));  nf = numel(vf);  end
        
    %DEFINE OPTIMIZER CONSTANTS
    ova=ones(1,nf); ovb=ones(1,ntp1); cam.msv.vf=vf;
    cam.msv.ucx=cam.uc.x(:,i);  cam.msv.ucy=cam.uc.y(:,i);  cam.msv.ucz=cam.uc.z(:,i); 
    D = cam.aposteriori.ned(vf,:);  cam.msv.A=D;  cam.msv.Ax=D(:,1);  cam.msv.Ay=D(:,2);  cam.msv.Az=D(:,3);
    cam.msv.nf = nf; %number of frames
    cam.msv.ntp = ntp1; %number of tie points
    cam.msv.ova = ova;
    cam.msv.ovb = ovb;
    %DEFINE PERMUTATIONS
    x = 1:nf;
    j = tril(x(ova,:), -1); j=j(j~=0); cam.msv.j=uint16(j);
    k = tril(x(ova,:)',-1); k=k(k~=0); cam.msv.k=uint16(k); cam.msv.njk=numel(j);
    %DEFINE VECTOR ORIGINS
    BAx=D(j,1)-D(k,1);  BAy=D(j,2)-D(k,2);  BAz=D(j,3)-D(k,3);  
    cam.msv.BAx=BAx(:,ovb);  cam.msv.BAy=BAy(:,ovb);  cam.msv.BAz=BAz(:,ovb);
    cam.msv.zx = [a.upx(i,vf); a.upy(i,vf)];

    %GET TP XYZ
    [~, C0(i,:), rs] = fcnminrpyrpyrpy(cam,a,cam.aposteriori.rpy(vf,:)'*d2r);    msvr(i) = mean(sqrt(rs(:)))*1000;
end

end


%LOCAL FUNCTIONS ----------------------------------------------------------
function [res, C0, rs] = fcnminrpyrpyrpy(cam,a,x) %[rpyrpyrpy...] CENTERED
%LOAD CONSTANTS
j = cam.msv.j;  k = cam.msv.k;
vf = cam.msv.vf;  nf = cam.msv.nf;  ntp = cam.msv.ntp;  ovb = cam.msv.ovb;
A = cam.msv.A;  BAx=cam.msv.BAx; BAy=cam.msv.BAy; BAz=cam.msv.BAz;
rpy = reshape(x,[3 nf])';

%VECTOR ROTATIONS
r = rpy(:,1); cr=cos(r); sr=sin(r); cr=cr(:,ovb); sr=sr(:,ovb);
p = rpy(:,2); cp=cos(p); sp=sin(p); cp=cp(:,ovb); sp=sp(:,ovb);
y = rpy(:,3); cy=cos(y); sy=sin(y); cy=cy(:,ovb); sy=sy(:,ovb);
x0=cam.msv.ucx(vf,:); y0=cam.msv.ucy(vf,:); z0=cam.msv.ucz(vf,:); %x0=cam.uc.x(vf,:); y0=cam.uc.y(vf,:); z0=cam.uc.z(vf,:);
ux1=x0.*(cp.*cy) +y0.*(sr.*sp.*cy-cr.*sy)  +z0.*(cr.*sp.*cy+sr.*sy);  %[ux1,uy1,uz1] = fcnrotateW2Brpyxyz(sr,sp,sy,cr,cp,cy,x,y,z);
uy1=x0.*(cp.*sy) +y0.*(sr.*sp.*sy+cr.*cy)  +z0.*(cr.*sp.*sy-sr.*cy);
uz1=x0.*(-sp)    +y0.*(sr.*cp)             +z0.*(cr.*cp);
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
rs =  (t1.*vx-BAx-s1.*ux).^2 + (t1.*vy-BAy-s1.*uy).^2 + (t1.*vz-BAz-s1.*uz).^2 ; %sum of the squared ranges of the misclosure vectors

% %TIE POINT CENTERS
% den = cam.msv.njk*2; %denominator = number of permutations times 2
% D = sum(A)*(cam.msv.nf-1);
% C0x = (sum(t1.*vx+s1.*ux)+D(1)) / den;
% C0y = (sum(t1.*vy+s1.*uy)+D(2)) / den;
% C0z = (sum(t1.*vz+s1.*uz)+D(3)) / den;
% C0 = zeros(ntp,3); C0(:,1)=C0x; C0(:,2)=C0y; C0(:,3)=C0z;

v1 = 1-ux1.^2;      v2 = -ux1.*uy1;     v3 = -ux1.*uz1;
v4 = v2;            v5 = 1-uy1.^2;      v6 = -uy1.*uz1;
v7 = v3;            v8 = v6;            v9 = 1-uz1.^2;

S1mat = zeros(9,ntp);
S1mat(1,:) = sum(v1);  S1mat(2,:) = sum(v2);  S1mat(3,:) = sum(v3);
S1mat(4,:) = sum(v4);  S1mat(5,:) = sum(v5);  S1mat(6,:) = sum(v6);
S1mat(7,:) = sum(v7);  S1mat(8,:) = sum(v8);  S1mat(9,:) = sum(v9);

S2mat = zeros(3,ntp);  Ax = A(:,1)';  Ay = A(:,2)';  Az = A(:,3)';
S2mat(1,:) = Ax*v1 + Ay*v2 + Az*v3;
S2mat(2,:) = Ax*v4 + Ay*v5 + Az*v6;
S2mat(3,:) = Ax*v7 + Ay*v8 + Az*v9;

C0 = zeros(ntp,3);
S1m = reshape(S1mat,3,3,ntp);
for j=1:ntp; %tpi, cam.msv.ntp
    C0(j,:) = S1m(:,:,j)\S2mat(:,j);
end
%C0x=C0(:,1)';  C0y=C0(:,2)';  C0z=C0(:,3)';

%IMAGE ANGLE RESIDUALS
%ova = cam.msv.ova;  CAx = C0x(ova,:)-A(:,ovb);  CAy = C0y(ova,:)-A(:,ovb*2);  CAz = C0z(ova,:)-A(:,ovb*3);
%ts1 = cam.msv.nf*ntp - sum(sum(  (ux1.*CAx + uy1.*CAy + uz1.*CAz).^2./(CAx.*CAx+CAy.*CAy+CAz.*CAz)  ));  %sum of sqrt of angle sines  %sin = sqrt(1-ct^2)
%ts1 = asin(sqrt(abs(  1 - (ux1.*CAx + uy1.*CAy + uz1.*CAz).^2./(CAx.*CAx+CAy.*CAy+CAz.*CAz)  )));  %sum of sqrt of angle sines  %sin = sqrt(1-ct^2)
%ts1 = 1 - (ux1.*CAx + uy1.*CAy + uz1.*CAz).^2./(CAx.*CAx+CAy.*CAy+CAz.*CAz);  %sum of sqrt of angle sines  %sin = sqrt(1-ct^2)
%res = ts1;

%IMAGE PIXEL RESIDUALS
%zx = [a.upx(:,vf); a.upy(:,vf)];
res = zeros(ntp*2,nf); %residuals
% for i=1:nf
%     ui = [C0x-A(i,1); C0y-A(i,2); C0z-A(i,3)]';
%     uic = fcnrotateB2W(r(i),p(i),y(i),ui);
%     z1 = camcc2pixel(cam, uic);
%     res(:,i) = cam.msv.zx(:,i)-z1(:); 
% end
end