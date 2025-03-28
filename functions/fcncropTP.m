% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function a = fcncropTP(cam,a,i)

a.upx = a.upx(i,:);  
a.upy = a.upy(i,:);  
a.state = a.state(i,:);  
a.iis = a.iis(i,:);  
a.score = a.score(i);  

if cam.syntheticVideoFlag
    a.ipned = a.ipned(i,:);  
end




