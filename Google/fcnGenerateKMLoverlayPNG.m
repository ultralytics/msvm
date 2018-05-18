function s = fcnGenerateKMLoverlayPNG(pathname, picfname, lats, lngs, z)
%This function creates a PNG suitable for use as a KML overlay. The PNG is
%Plate Carree projected (cylindrical).
%figure; pcolor(lngs,lats,z); shading flat;
h1 = findobj('type','figure');

extent = [lats([1 end]) lngs([1 end])];
z = z-min(min(z));
z = z/max(max(z));
dmax = max(extent(4)-extent(3),extent(2)-extent(1)); %deg
if dmax>350 %must be full earth
    dx = dmax/1000; %space 1000 pixels along the largest dimension
else
    dx = dmax/500; %space 500 pixels along the largest dimension
end

xiv = extent(2):-dx:extent(1); %lat order is reversed because row order is reversed with imwrite
yiv = extent(3):dx:extent(4);

[xi,yi] = meshgrid(xiv, yiv);
[tr, tc] = size(lats);
if  tr==1 || tc==1 %need to meshgrid these 1d vectors
    [lats, lngs] = meshgrid(lats,lngs);
    z = z';
end

cmap = jet;
zi = interp2(lats,lngs,z,xi,yi,'linear');
imwrite(uint8(zi'*63),cmap,[pathname filesep picfname],'png')
s = fcnGenerateKMLoverlay(pathname, picfname, 'PNG Overlay', extent);

h2 = findobj('type','figure');  if numel(h2)>numel(h1);  close(h2(end));  end %blank figure may pop up during imwrite!
end


function s = fcnGenerateKMLoverlay(pathname, picfname, kmlname, extent)
dir = [pathname '/'];
dirForwardSlash = strrep(dir,'\','/');

s = '<GroundOverlay>';
s = [s 10, ['	<name>' kmlname '</name>']];
s = [s 10, '	<color>bdffffff</color>'];
s = [s 10, '	<Icon>'];
s = [s 10, ['		<href>' dirForwardSlash picfname '</href>']];
s = [s 10, '       <refreshMode>onInterval</refreshMode>'];
s = [s 10, '       <refreshInterval>1</refreshInterval>'];
s = [s 10, '		<viewBoundScale>0.75</viewBoundScale>'];
s = [s 10, '	</Icon>'];
s = [s 10, '	<LatLonBox>'];
s = [s 10, sprintf('		<north>%.6f</north>',extent(2))];
s = [s 10, sprintf('		<south>%.6f</south>',extent(1))];
s = [s 10, sprintf('		<east>%.6f</east>',extent(3))];
s = [s 10, sprintf('		<west>%.6f</west>',extent(4))];
s = [s 10, '	</LatLonBox>'];
s = [s 10, '</GroundOverlay>'];
end

% function [] = fcnGenerateKMLoverlay(pathname, kmlfname, picfname, kmlname, extent)
% %extent = [min(lat) max(lat) min(lng) max(lng)]
% fprintf(['Generating ''' kmlfname '''...'])
% dir = [pathname '/'];
% dirForwardSlash = strrep(dir,'\','/');
% 
% fid = fopen([dir kmlfname],'w');
% fprintf(fid,'<?xml version="1.0" encoding="UTF-8"?>\n');
% fprintf(fid,'<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">\n');
% fprintf(fid,'<GroundOverlay>\n');
% fprintf(fid,['	<name>' kmlname '</name>\n']);
% fprintf(fid,'	<color>bdffffff</color>\n');
% fprintf(fid,'	<Icon>\n');
% fprintf(fid,['		<href>' dirForwardSlash picfname '</href>\n']);
% fprintf(fid,'       <refreshMode>onInterval</refreshMode>');
% fprintf(fid,'       <refreshInterval>1</refreshInterval>');
% fprintf(fid,'		<viewBoundScale>0.75</viewBoundScale>\n');
% fprintf(fid,'	</Icon>\n');
% fprintf(fid,'	<LatLonBox>\n');
% fprintf(fid,'		<north>%.6f</north>\n',extent(2));
% fprintf(fid,'		<south>%.6f</south>\n',extent(1));
% fprintf(fid,'		<east>%.6f</east>\n',extent(3));
% fprintf(fid,'		<west>%.6f</west>\n',extent(4));
% fprintf(fid,'	</LatLonBox>\n');
% fprintf(fid,'</GroundOverlay>\n');
% fprintf(fid,'</kml>\n');
% fclose(fid);
% 
% fprintf('   Done.\n')
% end

