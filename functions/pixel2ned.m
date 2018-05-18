function ned = pixel2ned(DEM,cam,i,xy,flag)
ecef = pixel2ecef(DEM,cam,i,xy,flag);
ned = ecef2ned(DEM,ecef);
end

