% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

clc; clear all
syms r1 r2 x1 y1 z1 x2 y2 z2 BAx BAy BAz real
syms Ax Ay Az Bx By Bz Tx Ty Tz real
syms ATx ATy ATz BTx BTy BTz real


%CC2ELAZ
TAx = Tx-Ax; %T = Tie Point; A = Aircraft
TAy = Ty-Ay;
TAz = Tz-Az;
p1 = asin(-TAz./sqrt(TAx.^2 + TAy.^2 + TAz.^2)); %pointing angle to tp1
yaw1 = fcnatan2(TAy,TAx);  %pointing angle to tp1

%ROTATE
cr1=cos(r1); sr1=sin(r1);
cp1=cos(p1); sp1=sin(p1);
%cy1=1/(TAy^2/TAx^2 + 1)^.5;   sy1=TAy/(TAx*(TAy^2/TAx^2 + 1)^.5);  
cy1=cos(yaw1); sy1=sin(yaw1);
ux=x1.*(cp1.*cy1) +y1.*(sr1.*sp1.*cy1-cr1.*sy1)  +z1.*(cr1.*sp1.*cy1+sr1.*sy1);  %[ux1,uy1,uz1] = fcnrotateW2Brpyxyz(sr,sp,sy,cr,cp,cy,x,y,z);
uy=x1.*(cp1.*sy1) +y1.*(sr1.*sp1.*sy1+cr1.*cy1)  +z1.*(cr1.*sp1.*sy1-sr1.*cy1);
uz=x1.*(-sp1)     +y1.*(sr1.*cp1)                +z1.*(cr1.*cp1);

%CC2ELAZ
TBx = Tx-Bx; %T = Tie Point; A = Aircraft
TBy = Ty-By;
TBz = Tz-Bz;
p2 = asin(-TBz./sqrt(TBx.^2 + TBy.^2 + TBz.^2)); %pointing angle to tp1
yaw2 = fcnatan2(TBy,TBx);  %pointing angle to tp1

%ROTATE
cr2=cos(r2); sr2=sin(r2);
cp2=cos(p2); sp2=sin(p2);
%cy2=1/(TBy^2/TBx^2 + 1)^.5;   sy2=TBy/(TBx*(TBy^2/TBx^2 + 1)^.5);  
cy2=cos(yaw2); sy2=sin(yaw2);
vx=x2.*(cp2.*cy2) +y2.*(sr2.*sp2.*cy2-cr2.*sy2)  +z2.*(cr2.*sp2.*cy2+sr2.*sy2);  %[ux1,uy1,uz1] = fcnrotateW2Brpyxyz(sr,sp,sy,cr,cp,cy,x,y,z);
vy=x2.*(cp2.*sy2) +y2.*(sr2.*sp2.*sy2+cr2.*cy2)  +z2.*(cr2.*sp2.*sy2-sr2.*cy2);
vz=x2.*(-sp2)     +y2.*(sr2.*cp2)                +z2.*(cr2.*cp2);

%VECTOR INTERCEPTS
d = ux.*vx + uy.*vy + uz.*vz;
e = ux.*BAx + uy.*BAy + uz.*BAz;
f = vx.*BAx + vy.*BAy + vz.*BAz;
g = 1 - d.*d;
s = (d.*f - e)./g; %multiply times u
t = (f - d.*e)./g; %multiply times v

%MISCLOSURE VECTOR RANGE RESIDUALS
L = sqrt( (t.*vx-BAx-s.*ux).^2 + (t.*vy-BAy-s.*uy).^2 + (t.*vz-BAz-s.*uz).^2 );

dLdTx = diff(L,Tx); 

dLdTy = diff(L,Ty);

dLdTz = diff(L,Tz);

%dLdr1 = diff(L,r1)


%simpledLdTz = simple(dLdTz,'IgnoreAnalyticConstraints',true)
%dLdr3 = feval(symengine,'generate::optimize',[dLdTx]);
%dLdr3 = feval(symengine,'pretty',[dLdTx]);


   str = char(dLdTx);
   str = strrep(str,'==','=');
   str = strrep(str,'(1/2)','.5');
   str = strrep(str,'*','.*');
   str = strrep(str,'/','./');
   str = strrep(str,'^','.^');
   str = str(~isspace(str));
   str = strrep(str,'Ax-Tx','ATx');  str = strrep(str,'Tx-Ax','TAx');
   str = strrep(str,'Ay-Ty','ATy');  str = strrep(str,'Ty-Ay','TAy');
   str = strrep(str,'Az-Tz','ATz');  str = strrep(str,'Tz-Az','TAz');
   str = strrep(str,'Bx-Tx','BTx');  str = strrep(str,'Tx-Bx','TBx');
   str = strrep(str,'By-Ty','BTy');  str = strrep(str,'Ty-By','TBy');
   str = strrep(str,'Bz-Tz','BTz');  str = strrep(str,'Tz-Bz','TBz');
   str = strrep(str,'cos(r1)','cr1');  str = strrep(str,'sin(r1)','sr1');
   str = strrep(str,'cos(r2)','cr2');  str = strrep(str,'sin(r2)','sr2');
   
   fprintf('%s; ', str )
   
   dLdr3 = feval(symengine,'generate::optimize',eval(str));


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