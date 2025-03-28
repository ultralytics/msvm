% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [] = fcnvideo2images(varargin)
%fcnvideo2images() or fcnvideo2images(30:105)

[filename, pathname] = uigetfile({'*.avi;*.mpeg;*.mp4;*.wmv','Video Files (*.avi, *.mpeg, *.mp4, *.wmv)';'*.*','All files (*.*)'},'Select video file','MultiSelect','off');  
if filename==0; return; end;

newpath = [pathname filename(1:end-4) ' images\'];
mkdir(newpath)

vfr=VideoReader([pathname filename]);  nf=vfr.NumberOfFrames;  cam.frames=nf;  cam.width=vfr.Width;  cam.height=vfr.Height;  cam.fps=vfr.FrameRate;

frames = 1:nf;
if nargin==1
    frames = varargin{1};  nf = numel(frames);
end

fprintf('Converting ''%s'' into %.0f images...\n',filename,nf); startclock=clock;
for i = 1:nf
    j = frames(i);
    fprintf('image %.0f/%.0f\n',i,nf);
    rgb = read(vfr,j);
    imwrite(rgb,sprintf('%s %.0f.jpg',newpath,j))
end
fprintf('Done in %.1fs\n\nNew images placed in ''%s''\n',etime(clock,startclock),newpath);
