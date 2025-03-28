% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [ output_args ] = fcnGoogleStaticMaps(cam)
str1 = 'http://maps.googleapis.com/maps/api/staticmap?';
str2 = '&size=640x640&scale=1&maptype=satellite&sensor=false';

lla = ned2lla(cam.DEM,mean(cam.tpnedhat));
zoom = 16;

url = sprintf('%scenter=%.5f,%.5f&zoom=%.0f%s',str1,lla(1),lla(2),zoom,str2);

[M, Mcolor] = imread(url);
M = cast(M,'double');
width = size(M,2);
height = size(M,1);

imag = zeros(height,width,3);
for idx = 1:3
    imag(:,:,idx) = reshape(Mcolor(M(:)+1+(idx-1)*size(Mcolor,1)),height,width);
end
figure; imshow(imag)
   
%http://maps.googleapis.com/maps/api/staticmap?center=40.714728,-73.998672&zoom=12&size=640x640&scale=2&maptype=satellite&sensor=false
%http://maps.googleapis.com/maps/api/staticmap?center=40.714728,-73.998672&visible=40.714728,-73.998672|40.714728,-73.918672&size=640x640&maptype=satellite&scale=2&sensor=false
end

