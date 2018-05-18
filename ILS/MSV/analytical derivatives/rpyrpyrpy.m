clc; clear all
syms r p y x y z BAx BAy BAz vx vy vz real

el = asin(-cc(:,3)./sqrt(cc(:,1).^2 + cc(:,2).^2 + cc(:,3).^2));
az = atan2(cc(:,2),cc(:,1));


cr=cos(r); sr=sin(r);
cp=cos(p); sp=sin(p);
cy=cos(y); sy=sin(y);
%x=cam.u.x(vf,:); y=cam.u.y(vf,:); z=cam.u.z(vf,:);
ux=x.*(cp.*cy) +y.*(sr.*sp.*cy-cr.*sy)  +z.*(cr.*sp.*cy+sr.*sy);  %[ux1,uy1,uz1] = fcnrotateW2Brpyxyz(sr,sp,sy,cr,cp,cy,x,y,z);
uy=x.*(cp.*sy) +y.*(sr.*sp.*sy+cr.*cy)  +z.*(cr.*sp.*sy-sr.*cy);
uz=x.*(-sp)    +y.*(sr.*cp)             +z.*(cr.*cp);

%VECTOR INTERCEPTS
d = ux.*vx + uy.*vy + uz.*vz;
e = ux.*BAx + uy.*BAy + uz.*BAz;
f = vx.*BAx + vy.*BAy + vz.*BAz;
g = 1 - d.*d;
s1 = (d.*f - e)./g; %multiply times u
t1 = (f - d.*e)./g; %multiply times v

%MISCLOSURE VECTOR RANGE RESIDUALS
L = sqrt( (t1.*vx-BAx-s1.*ux).^2 + (t1.*vy-BAy-s1.*uy).^2 + (t1.*vz-BAz-s1.*uz).^2 );



