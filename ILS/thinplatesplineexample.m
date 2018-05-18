

nxy = 31;
xy = 2*(rand(2,nxy)-.5); 
vals = sum(xy.^2);
%vals = rand(1,nxy);
noisyvals = vals + (rand(size(vals))-.5)/5;

st = tpaps(xy,vals,.99); fnplt(st), hold on
avals = fnval(st,xy);
plot3(xy(1,:),xy(2,:),vals,'wo','markerfacecolor','k')


%quiver3(xy(1,:),xy(2,:),avals,zeros(1,nxy),zeros(1,nxy), noisyvals-avals,'r'), hold off