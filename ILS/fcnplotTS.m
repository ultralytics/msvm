% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [hf] = fcnplotTS(ts,tsv,nmc,xstr,cam)

[h, hf] = fig(3,4);
x = tsv;
c = {'r','g','b'};
af=.2; %alpha fraction
h1=zeros(3,1);
vx=1:numel(tsv);  
x=x(vx);

axes(h(3)); j=[1 3 5]; %#ok<*MAXES>
for i=1:3;
    k=j(i);  mu=ts.sigma12mean(vx,k);  s=ts.sigma12std(vx,k); 
    h1(i)=errorbar(x,ts.xhat12mean(vx,k),ts.xhat12std(vx,k),'.-','color',c{i},'markersize',1); hold on;
    fill([x fliplr(x)],[mu+s; flipud(mu-s)],c{i},'edgecolor','none'); 
    fill([x fliplr(x)],[-mu-s; flipud(-mu+s)],c{i},'edgecolor','none');
end
legend(h1,'x','y','z'); xlabel(xstr); ylabel('pos (m)'); title('aircraft {\bfposition} error'); alpha(af); box off

axes(h(7)); j=[7 9 11]; 
for i=1:3;
    k=j(i);  mu=ts.sigma12mean(vx,k);  s=ts.sigma12std(vx,k); 
    h1(i)=errorbar(x,ts.xhat12mean(vx,k),ts.xhat12std(vx,k),'.-','color',c{i},'markersize',1); hold on;
    fill([x fliplr(x)],[mu+s; flipud(mu-s)],c{i},'edgecolor','none'); 
    fill([x fliplr(x)],[-mu-s; flipud(-mu+s)],c{i},'edgecolor','none');
end
z=max(ts.sigma12mean(vx,j)+ts.sigma12std(vx,j),[],2); y=mean(ts.xhat12std(vx,j),2); for i=2:numel(x); text(x(i),z(i),sprintf('  %.2f^o',y(i)),'horizontalalignment','left','rotation',90,'fontsize',7); end; text(x(1),z(1),sprintf('  %.2f^o',y(1)),'horizontalalignment','left','verticalalignment','top','rotation',90,'fontsize',7);   
legend(h1,'roll','pitch','yaw','location','southeast'); xlabel(xstr); ylabel('rpy (deg)'); title('aircraft {\bfrpy} error'); alpha(af); box off

axes(h(4)); j=[2 4 6];
for i=1:3;
    k=j(i);  mu=ts.sigma12mean(vx,k);  s=ts.sigma12std(vx,k); 
    h1(i)=errorbar(x,ts.xhat12mean(vx,k),ts.xhat12std(vx,k),'.-','color',c{i},'markersize',1); hold on;
    fill([x fliplr(x)],[mu+s; flipud(mu-s)],c{i},'edgecolor','none'); 
    fill([x fliplr(x)],[-mu-s; flipud(-mu+s)],c{i},'edgecolor','none');
end
legend(h1,'x','y','z'); xlabel(xstr); ylabel('vel (m/s)'); title('aircraft {\bfvelocity} error'); alpha(af); box off

axes(h(8)); j=[8 10 12];
for i=1:3;
    k=j(i);  mu=ts.sigma12mean(vx,k);  s=ts.sigma12std(vx,k); 
    h1(i)=errorbar(x,ts.xhat12mean(vx,k),ts.xhat12std(vx,k),'.-','color',c{i},'markersize',1); hold on;
    fill([x fliplr(x)],[mu+s; flipud(mu-s)],c{i},'edgecolor','none'); 
    fill([x fliplr(x)],[-mu-s; flipud(-mu+s)],c{i},'edgecolor','none');
end
legend(h1,'roll','pitch','yaw'); xlabel(xstr); ylabel('rpy rate (deg/s)'); title('aircraft {\bfrpy} rate error'); alpha(af); axis tight; box off

axes(h(11)); j=[1 2 3];
for i=1:3;
    k=j(i);  mu=ts.sigma3mean(vx,k);  s=ts.sigma3std(vx,k);
    h1(i)=errorbar(x,ts.xhat3mean(vx,k),ts.xhat3std(vx,k),'.-','color',c{i},'markersize',1); hold on;
    fill([x fliplr(x)],[mu+s; flipud(mu-s)],c{i},'edgecolor','none'); 
    fill([x fliplr(x)],[-mu-s; flipud(-mu+s)],c{i},'edgecolor','none');
end
z=max(ts.sigma3mean(vx,j)+ts.sigma3std(vx,j),[],2); y=mean(ts.xhat3std(vx,j),2); for i=2:numel(x); text(x(i),z(i),sprintf('  %.1fm',y(i)),'horizontalalignment','left','rotation',90,'fontsize',7); end; text(x(1),z(1),sprintf('  %.1fm',y(1)),'horizontalalignment','left','verticalalignment','top','rotation',90,'fontsize',7);   
legend(h1,'x','y','z','location','southeast'); xlabel(xstr); ylabel('pos (m)'); title('{\bftie point} error'); alpha(af); box off

axes(h(6));  j=[4 5 6];
for i=1:3;
    k=j(i);  mu=ts.initRMSEmean(vx,k);  s=ts.initRMSEstd(vx,k);
    h1(i)=errorbar(x,mu,s,'.-','color',c{i},'markersize',1); hold on;
end
legend(h1,'roll','pitch','yaw','Location','best'); xlabel(xstr); ylabel('rpy (deg)'); title('{\bfinitial} aircraft rpy RSME'); axis tight; box off
z=max(ts.initRMSEstd(vx,j)+ts.initRMSEmean(vx,j),[],2); y=mean(ts.initRMSEmean(vx,j),2); for i=2:numel(x); text(x(i),z(i),sprintf('  %.2f^o',y(i)),'horizontalalignment','left','rotation',90,'fontsize',7); end; text(x(1),z(1),sprintf('  %.2f^o',y(1)),'horizontalalignment','left','verticalalignment','top','rotation',90,'fontsize',7);   

axes(h(10));  j=[1 2 3];
for i=1:3;
    k=j(i);  mu=ts.initRMSEmean(vx,k);  s=ts.initRMSEstd(vx,k);
    h1(i)=errorbar(x,mu,s,'.-','color',c{i},'markersize',1); hold on;
end
legend(h1,'x','y','z','location','best'); xlabel(xstr); ylabel('pos (m)'); title('{\bfinitial} tie point RSME'); axis tight; box off
z=max(ts.initRMSEstd(vx,j)+ts.initRMSEmean(vx,j),[],2); y=mean(ts.initRMSEmean(vx,j),2); for i=2:numel(x); text(x(i),z(i),sprintf('  %.1fm',y(i)),'horizontalalignment','left','rotation',90,'fontsize',7); end; text(x(1),z(1),sprintf('  %.1fm',y(1)),'horizontalalignment','left','verticalalignment','top','rotation',90,'fontsize',7);   

axes(h(12))
y=ts.ilssuccessfraction(vx)*100;  bar(x(vx),y,.5,'b','edgecolor','none');
xlabel(xstr); ylabel('ILS Success Rate'); title(sprintf('ILS convergence rate (based on %.0f MCs)',nmc)); axis tight; set(gca,'ylim',[min(y)-2 100.5]); box off
ytick=unique(round(linspace(floor(min(y)-2),100,10)));  set(gca,'ytick',ytick); for i=1:numel(x); text(x(i),y(i),sprintf(' %.4g%%',y(i)),'horizontalalignment','left','rotation',90,'fontsize',7);  end

axes(h(9))
y=ts.timeperMC(vx,:);  hb=bar(x(vx),y,1,'b','edgecolor','none'); set(hb(1),'facecolor','r')
xlabel(xstr); ylabel('time (s)'); title('MSV and ILS time per MC'); axis tight; set(gca,'ylim',[0 max(y(:))*1.2]); box off; legend('MSV','ILS','location','best')
z=max(y,[],2); y=sum(y,2); for i=1:numel(x); text(x(i),z(i),sprintf('  %.1fs',y(i)),'horizontalalignment','left','rotation',90,'fontsize',7); end;

axes(h(1))
axis off

axes(h(2))
axis off

axes(h(5))
axis off
cam.meanGSD = mean(cam.true.focus.gsd);
a = fieldnames(cam);
b = struct2cell(cam);
units = {'','(pix)','(pix)','(deg)','(deg)','(km)','(s)','','','','(m)','(deg)','(deg)','(m)','(deg)','(s)','','(m)'};
str=[];
j=0;
for i=1:numel(a)
    if isnumeric(b{i}) && numel(b{i})==1
        j=j+1;
        str = [str sprintf('%16s =%9.5g %-10s',a{i},b{i},units{j})];
        str = [str sprintf('\n')];
    end
end
text(-.18,-.1,str,'FontName','MonoSpaced','FontSize',8.5,'VerticalAlignment','bottom','HorizontalAlignment','Left')

annotation(gcf,'line',[1 1]*.492,[0 1],'linewidth',5,'color',[.8 .8 .8])
annotation(gcf,'textbox',[.2 .916 .1 .1],'String','MSV','color',[.8 .8 .8],'edgecolor','none','fontsize',40,'fontweight','bold')
annotation(gcf,'textbox',[.71 .916 .1 .1],'String','ILS','color',[.8 .8 .8],'edgecolor','none','fontsize',40,'fontweight','bold')

fcntight(h,'x')
end

