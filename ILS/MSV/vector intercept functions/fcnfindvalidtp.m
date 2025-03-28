% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function j = fcnfindvalidtp(a,vf)

fx = a.upx(:,vf);
fy = a.upy(:,vf);


j = ~any(a.state(:,vf)==0 | fx==0 | fy==0 | isnan(fx) | isnan(fy) ,  2);     


