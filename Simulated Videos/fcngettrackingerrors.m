% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function rms = fcngettrackingerrors(a,cam)
%load DEMsolclose;  cam.DEM = DEM;

minlifespan = 4; %frames
valid = sum(a.state==1,2) > minlifespan; %track length > 1
a.upx = a.upx(valid,:);
a.upy = a.upy(valid,:);
a.score = a.score(valid,:);
a.state = a.state(valid,:);

x=1:cam.frames;
np = sum(valid);
nf = cam.frames;
gsd = mean(cam.true.focus.gsd);
valid = a.state == 1; %track length > 1

% %GET NED ERRORS -----------------------------------------------------------
% ecef = zeros([size(a.upx,1) 3 size(a.upx,2)]);
% for i=1:nf
%     f = [a.upx(:,i) a.upy(:,i)];
%     ecef(:,:,i) = pixel2ecef(cam.DEM,cam,i,f,'true');
% end
% diff = bsxfun(@minus,ecef(:,:,1:end),ecef(:,:,1))*1000;
% diff = reshape( sqrt(mean(diff.^2,2)),  [np nf]);
% rms.meters = fcnrms(diff(:));
% rms.pixels = rmsmeters/gsd;

%GET PIXEL ERRORS ---------------------------------------------------------
pixeldiff = zeros(np,nf);
a.ipned = zeros(np,3);
for i= 1:np
    j = find(a.state(i,:)==1,1,'first');
    a.ipned(i,:) = ecef2ned(cam.DEM, pixel2ecef(cam.DEM,cam,cam.frameID(j),[a.upx(i,j) a.upy(i,j)],'true'));
end


mpd = zeros(nf,1); %mean pixel difference
for i=1:nf
    dx = [a.upx(:,i) a.upy(:,i)] - ned2pixel(cam,cam.frameID(i),a.ipned,'true');
    pixeldiff(:,i) = sqrt(dx(:,1).^2 + dx(:,2).^2);
    mpd(i) = mean(pixeldiff(valid(:,i),i));
end
pixeldiff(~valid) = nan;
rms.pixels = fcnrms(pixeldiff(valid));
rms.meters = rms.pixels*gsd;

%PLOT ---------------------------------------------------------------------
h1=plot(x,pixeldiff);
h2=plot(x,mpd,'b','linewidth',3,'Color',[.7 .7 .7]); uistack(h2,'top');

xlabel('video frame'); ylabel('RMSE (pixels)'); 
title(sprintf('Tracking Errors\nRMSE=%.2fm, or %.1fpixels at %.2fm GSD for %.0f tracks > %.0f frames',...
    rms.meters,rms.pixels,gsd,np,minlifespan)); set(gca,'xlim',[0 nf*1.03],'ylim',[0 max3(pixeldiff)]);
legend([h1(1) h2],'RMSE','mean RMSE','Location','Best')

% %LABEL EACH TRACK ---------------------------------------------------------
% cmap = jet;
% color = interp1(linspace(min(a.score),max(a.score),size(cmap,1))',cmap,a.score);
% for i=1:np
%     text(nf,pixeldiff(i,nf),sprintf(' %.0f',i),'fontsize',8,'color',color(i,:))
%     set(h1(i),'Color',color(i,:));
% end
% box off
