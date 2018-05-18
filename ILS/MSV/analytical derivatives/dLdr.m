clc; clear all
syms r p yaw x y z BAx BAy BAz vx vy vz real

cr=cos(r); sr=sin(r);
cp=cos(p); sp=sin(p);
cy=cos(yaw); sy=sin(yaw);
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

dLdr1 = diff(L,r);

%simpledLdr1 = simple(dLdr1,'IgnoreAnalyticConstraints',true)
%dLdr2 = feval(symengine,'subsex',dLdr1,[char(cr) '=cr'])

clc
%dLdr2 = subs(dLdr1,{char(cr),char(sr),char(cp),char(sp),char(cy),char(sy)},{'cr','sr','cp','sp','cy','sy'});
dLdr3 = feval(symengine,'generate::optimize',dLdr1);
for i=1:numel(dLdr3)
   str = char(dLdr3(i));
   str = strrep(str,'==','=');
   str = strrep(str,'(1/2)','.5');
   str = strrep(str,'*','.*');
   str = strrep(str,'/','./');
   str = strrep(str,'^','.^');
   str = str(~isspace(str));
   fprintf('%s; ', str )
end
    