% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [] = fcnplotFeatureTracks(a,cam,originalFrame,finalFrame)
ha = fig(2,2,1.5,2,[2, 3, 1.9, 0.9, 1.9, 1.2]);

sca(ha(1))
v1=a.state==1;
pcolor(double(~v1))
shading flat; axis normal on tight
xlabel('video frame');  ylabel('tie point');  title(sprintf('Tie Point Lifetimes\nmean lifespan = %.1f frames (%.0f full-length tracks)',sum3(v1)/size(v1,1),sum(a.state(:,1)==1 & a.state(:,end)==1)));
colormap(gca,gray); box off

sca(ha(2))
bar(1:cam.frames,double(sum(v1,1)),1,'b','edgecolor',[.7 .7 1]);
xlabel('video frame');  ylabel('active tie points');  title('Active Tie Points')
fcntight(ha(1:2),'x')

sca(ha(3))
if cam.syntheticVideoFlag
    rmse = fcngettrackingerrors(a,cam);
end

sca(ha(4))
[X,Y] = meshgrid([1 cam.width],[1 cam.height]); Z=ones(size(X));
surf(X,Y,Z*1,'FaceColor','Texture','EdgeColor','none','CData',double(originalFrame),'FaceAlpha',.5);
surf(X,Y,Z*cam.frames,'FaceColor','Texture','EdgeColor','none','CData',double(finalFrame),'FaceAlpha',.5);
axis equal vis3d tight; box on; grid off
xyzlabel('x (pixels)','y (pixels)','frame');

length = sum(v1,2);
x=a.upx;
y=a.upy;
z=ones(size(x(:,1)))*(1:cam.frames); z(~v1)=nan;
ml = 2; %minimum length
v2=find(length>ml);  nv = numel(v2);  nanv=nan(1,nv);
x=x(v2,:)';  y=y(v2,:)';  z=z(v2,:)';  length=length(v2);

flagColorTracks = false;
colormap(gca,jet)
if flagColorTracks
    cx = linspace(1,max(length),size(colormap,1));
    for i=1:nv
        color = interp1(cx,colormap,length(i));
        plot3(x(:,i),y(:,i),z(:,i),'-','linewidth',.2,'color',color)
    end
else
    plot3(fcncol(cat(1,x,nanv)),fcncol(cat(1,y,nanv)),fcncol(cat(1,z,nanv)),'-','linewidth',.2,'color',[.7 .7 .7])
end

plot3(x(1,:),y(1,:),z(1,:),'.','linewidth',.2,'color','b','markersize',10)
plot3(x(cam.frames,:),y(cam.frames,:),z(cam.frames,:),'.','linewidth',.2,'color','r','markersize',10)
set(gca,'clim',fcnminmax(length))

colorbar
%fcncolorbar(1,' frames')
set(gca,'zdir','reverse'); axis ij
title(sprintf('Tie Point Lifetimes, max=%.0f frames\n(showing all %.0f tracks longer than %.0f frames)',max(length),numel(length),ml))
daspect([1 1 cam.frames/((cam.width+cam.height)/2)])
set(gca,'cameraviewangle',9.5)
view(-24,34);

