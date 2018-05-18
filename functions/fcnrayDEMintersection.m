function x = fcnrayDEMintersection(DEM,orig,dest)
np = size(dest,1);  x=zeros(np,3);

%RAY INTERSECTIONS --------------------------------------------------------
faces = DEM.delaunayfaces; %delaunay(DEM.lng,DEM.lat);
vertices = DEM.ecef;  % vertices stored as Nx3 matrix
vert1 = vertices(faces(:,1),:);
vert2 = vertices(faces(:,2),:);
vert3 = vertices(faces(:,3),:);

dir = bsxfun(@minus,dest,orig);
[i, r]=TriangleRayIntersection(orig, dir, vert1, vert2, vert3);  %fprintf('Number of: faces=%i, points=%i, intersections=%i; time=%f sec\n', size(faces,1), size(vertices,1), sum(i), toc);

if any(i>0)
    x(i,1) = orig(1)+r(i).*dir(i,1);
    x(i,2) = orig(2)+r(i).*dir(i,2);
    x(i,3) = orig(3)+r(i).*dir(i,3);
end

if any(i==0)
    fprintf('WARNING: POINTS NOT ON DEM, EXTRAPOLATING. UPDATE DEM FOR BETTER ACCURACY\n')
    DEM = extrapolateDEM(DEM);
    
    %RAY INTERSECTIONS --------------------------------------------------------
    faces = DEM.delaunayfaces; %delaunay(DEM.lng,DEM.lat);
    vertices = DEM.ecef;  % vertices stored as Nx3 matrix
    vert1 = vertices(faces(:,1),:);
    vert2 = vertices(faces(:,2),:);
    vert3 = vertices(faces(:,3),:);
    
    dir = bsxfun(@minus,dest,orig);
    [i, r]=TriangleRayIntersection(orig, dir, vert1, vert2, vert3);  %fprintf('Number of: faces=%i, points=%i, intersections=%i; time=%f sec\n', size(faces,1), size(vertices,1), sum(i), toc);
    
    x=zeros(np,3);
    x(i,1) = orig(1)+r(i).*dir(i,1);
    x(i,2) = orig(2)+r(i).*dir(i,2);
    x(i,3) = orig(3)+r(i).*dir(i,3);
end

% fig;
% trimesh(DEM.delaunayfaces,DEM.ned(:,1),DEM.ned(:,2),DEM.ned(:,3),DEM.nedz(:)); hold on
% o = repmat(orig,[np 1]);
% plot3(orig(:,1),orig(:,2),orig(:,3),'o','MarkerSize',10)
% plot3(dest(:,1),dest(:,2),dest(:,3),'k.','MarkerSize',2)
% plot3([o(:,1) x(:,1)]',[o(:,2) x(:,2)]',[o(:,3) x(:,3)]','-b'); axis equal vis3d
% plot3(x(:,1),x(:,2),x(:,3),'.','MarkerSize',10)
end




function DEM = extrapolateDEM(DEM)
ni = numel(DEM.lat(:,1));
DEM.lat(:,1)=DEM.lat(:,1)-1;
DEM.lat(:,ni)=DEM.lat(:,ni)+1;
DEM.lng(1,:)=DEM.lng(1,:)-1;
DEM.lng(ni,:)=DEM.lng(ni,:)+1;

%ORIGINAL
DEM.latv=DEM.lat(:);
DEM.lngv=DEM.lng(:);
DEM.altv=DEM.alt(:);

DEM.lla = [DEM.latv DEM.lngv DEM.altv];

i=round(ni/2);  
centerlla=[DEM.lat(i,i) DEM.lng(i,i) 0];  
DEM.centerecef=lla2ecef(centerlla); %ok

DEM.DCM_ECEF2NED = fcnLLA2DCM_ECEF2NED(centerlla*d2r); %NOT OK!!!

DEM.ecef = lla2ecef(DEM.lla); %ok
DEM.ned = [DEM.ecef(:,1)-DEM.centerecef(1), DEM.ecef(:,2)-DEM.centerecef(2), DEM.ecef(:,3)-DEM.centerecef(3)]*DEM.DCM_ECEF2NED'; %NOT OK
DEM.nedx = reshape(DEM.ned(:,1),[ni ni]);
DEM.nedy = reshape(DEM.ned(:,2),[ni ni]);
DEM.nedz = reshape(DEM.ned(:,3),[ni ni]);

DEM.delaunayfaces = delaunay(DEM.lng,DEM.lat);

DEM.F=scatteredInterpolant(DEM.lngv, DEM.latv, DEM.altv); %ok
DEM.Fned=scatteredInterpolant(DEM.ned(:,1), DEM.ned(:,2), DEM.altv); %ok
end