% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function ned = pixel2ned(DEM,cam,i,xy,flag)
ecef = pixel2ecef(DEM,cam,i,xy,flag);
ned = ecef2ned(DEM,ecef);
end

