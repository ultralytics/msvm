% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function fcnMSVMimages()
%cam.dateformat = 'HH:MM:SS';
clc

[file, path] = uigetfile({'*.jpg','JPG Files (*.jpg)';'*.*','All files (*.*)'},'Select images','MultiSelect','on');
if isnumeric(path); return; end; pathnow=what; addpath(pathnow.path); addpath(fcnpathm1(path)); addpath(path);  cd(path);
if ~iscell(file); x=file; clear file; file{1} = x; end
ni = numel(file);

vidObj = VideoWriter('MSVM image video'); open(vidObj); %'Uncompressed AVI' or 'MPEG-4' or 'Motion JPEG AVI'

t0 = clock; t0=t0(4:6)';
camdone = false;
lla = zeros(ni,4); %llat
valid = false(ni,1);
for i=1:numel(file)
    fprintf('image %d\n',i)
    a = imfinfo(file{i});
    
    %if isfield(a,'GPSInfo') && all(isfield(a.GPSInfo,{'GPSLatitude','GPSLongitude','GPSAltitude','GPSTimeStamp'}))
     if isfield(a,'GPSInfo') && all(isfield(a.GPSInfo,{'GPSLatitude','GPSLongitude','GPSAltitude'}))
        valid(i) = true;

        
        b = a.GPSInfo;
        t0=t0 + [0 0 1/3/86400]'; b.GPSTimeStamp = t0; %3Hz
        if ~camdone
            pixelPitch = 1.12E-3; %mm on a side
            cam.exif = a;
            cam.width = a.Width;
            cam.height = a.Height;
            cam.focalLength = a.DigitalCamera.FocalLength / pixelPitch;
            cam.fovh = atand(cam.width/2/cam.focalLength)*2; %deg
            cam.fovv = atand(cam.height/2/cam.focalLength)*2; %deg
            camdone = true;
        end
        
        lla(i,1:3) = [dms2deg(b.GPSLatitude) dms2deg(b.GPSLongitude) b.GPSAltitude].*fcnllsign(b);
        lla(i,4) = [3600 60 1]*b.GPSTimeStamp; %seconds after midnight
        
        str = sprintf('%s%s',path,file{i}); %fid = fopen(str);  if fid==-1; fprintf('Failure, file ''%s'' not found!\n',str); end
        %raw = fread(fid, [a.Width a.Height], 'uint8=>uint8')';  fclose(fid);
        %raw = fread(fid, inf, '*uint8')';  figure; imshow(raw)
        %fclose(fid);
        
        raw = imread(str);  
        %raw0 = flipdim(imrotate(imresize(raw,scaling,'bicubic'),-90), 2);    %cla; imshow(raw); axis on; colorbar; drawnow;
        writeVideo(vidObj,raw);
    end
end
badfiles = file(~valid);
close(vidObj);

lla = lla(valid,:);  file = file(valid);  nj = numel(file); %remove invalid files
[~, j] = sort(lla(:,4),'ascend');  lla = lla(j,:);  file = file(j); %sort by time

fprintf('All done, %d of %d images loaded successfully.\n',nj,ni)
if ni~=nj
    fprintf('\nThe following images lacked proper exif GPS tags:\n'); disp(char(badfiles))
end

if exist('cam','var');  save llat.mat cam lla file;  else  fprintf('\nWARNING: No data saved.\n');  end
end


%FUNCTIONS ----------------------------------------------------------------
function x=fcnllsign(b)
x=[1 1 1];
if b.GPSLatitudeRef=='S'; x(1)=-1; end
if b.GPSLongitudeRef=='W'; x(2)=-1; end
end

