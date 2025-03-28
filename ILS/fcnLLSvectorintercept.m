% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function tphat = fcnLLSvectorintercept(cam,a)
activetp = a.upx~=0 | a.upy~=0;

nframes = cam.frames; 
ntp = find(any(a.upy~=0 & a.upx~=0,2),1,'last'); %number of tie points
tphat = zeros(ntp,3);
for tpi = 1:ntp
    nf = max(ceil(nframes/10),2); %number of vectors to intercept
    frames = round(linspace(1,nframes,nf));
    zned1 = zeros(nf,3);  zhat1 = zned1;
    
    camned1 = cam.apriori.ned(frames,:);
    H = zeros(nf*3,3); %90x1 = 90x3 * 3x1, z=Hx
    xhat = mean(cam.apriori.ned(frames,:)); xhat(3)=0; %initial guess
    for iterations = 1:3
        for i=1:nf
            j = frames(i);
            vt = find(activetp(:,j)); %valid tp number indices
            z1 = [a.upx(vt,j), a.upy(vt,j)];
            zned1(i,:) = pixel2neduvec(cam,j,z1(tpi,:),'aposteriori');
            
            [r1, dx1] = fcnrange(camned1(i,:),xhat);
            zhat1(i,:) = dx1./r1;
            
            % syms tx cx ty cy tz cz
            % zhatx = (tx-cx)/sqrt((tx-cx)^2 + (ty-cy)^2 + (tz-cz)^2);
            % zhaty = (ty-cy)/sqrt((tx-cx)^2 + (ty-cy)^2 + (tz-cz)^2);
            % zhatz = (tz-cz)/sqrt((tx-cx)^2 + (ty-cy)^2 + (tz-cz)^2);
            % dzdx1 = simplify( [diff(zhatx,tx), diff(zhatx,ty), diff(zhatx,tz); diff(zhaty,tx), diff(zhaty,ty), diff(zhaty,tz); diff(zhatz,tx), diff(zhatz,ty), diff(zhatz,tz)] );
            
            tx=xhat(1);      ty=xhat(2);      tz=xhat(3); %tie point pos xhat
            cx=camned1(i,1); cy=camned1(i,2); cz=camned1(i,3); %camera pos
            dx = (cx - tx);
            dy = (cy - ty);
            dz = (cz - tz);
            k1 = (dx^2 + dy^2 + dz^2)^(3/2);
            
            dzdx =   [(cy^2 - 2*cy*ty + cz^2 - 2*cz*tz + ty^2 + tz^2),                                                   -(dx*dy),                                            -(dx*dz)
                -(dx*dy),    (cx^2 - 2*cx*tx + cz^2 - 2*cz*tz + tx^2 + tz^2),                                            -(dy*dz)
                -(dx*dz),                                           -(dy*dz),     (cx^2 - 2*cx*tx + cy^2 - 2*cy*ty + tx^2 + ty^2)]./k1;
            
            H((i-1)*3+[1 2 3],1:3) = dzdx;
        end
        residuals = reshape((zned1-zhat1)',[nf*3 1]);
        %R = eye(nf*3);
        %HR = H'/R;
        HR = H'; %assumes R = I matrix
        xhat = xhat + ((HR*H)\HR*residuals)';
    end
    tphat(tpi,:) = xhat;
end
fprintf('Vector Intercept rms error (m) = %.2fm\n\n',fcnrms(tphat-a.ipned)*1000)

% figure
% quiver3(camned1(:,1),camned1(:,2),camned1(:,3),zned1(:,1),zned1(:,2),zned1(:,3),10,'b'); box on; axis equal; axis vis3d; set(gca,'zdir','reverse'); xlabel('x'); ylabel('y'); zlabel('z'); hold on;
% plot3(xhat(1),xhat(2),xhat(3),'k.','markersize',20)
end










