clc; clear all
%syms r p y x y z BAx BAy BAz u v w real
%vx vy vz = u v w
%ux uy uz = a b c
%BAx BAy BAz = A B C
syms r p y x y z A B C u v w real At 


cr=cos(r); 
sr=sin(r);
cp=cos(p); 
sp=sin(p);
cy=cos(y); 
sy=sin(y);
a=x*(cp*cy) +y*(sr*sp*cy-cr*sy)  +z*(cr*sp*cy+sr*sy);
b=x*(cp*sy) +y*(sr*sp*sy+cr*cy)  +z*(cr*sp*sy-sr*cy);
c=x*(-sp)   +y*(sr*cp)           +z*(cr*cp);

d = a*u + b*v + c*w;
e = a*A + b*B + c*C;
f = u*A + v*B + w*C;
g = 1 - d*d;
s = (d*f - e)/g;
t = (f - d*e)/g;

L = sqrt( (t*u-A-s*a)^2 + (t*v-B-s*b)^2 + (t*w-C-s*c)^2 );
dLdr = diff(L,r);