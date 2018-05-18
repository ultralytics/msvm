function cam = fcngetcamcornersLLA(cam,msvned,nlsned)
if ~isfield(cam,'pathname') || isempty(cam.pathname)
    try
        cam.pathname = fcnfile2folder([cam.filename ' Features.mp4']);
    catch
        cam.pathname = fcnfile2folder([cam.filename ' Features.avi']);
    end
end
EGM = load('EGM96single.mat');

s2day = 1/86400;
dateformat = 'yyyy-mm-ddTHH:MM:SSZ';
vf = 1:cam.frames;
xy = [1 1
    cam.width 1
    cam.width cam.height
    1 cam.height]; %clockwise from top left


%AIRCRAFT PATH
cam.google.kml.aircraftpath = aircraftpath(cam,vf);


%DEM
cam.google.kml.DEM = DEM(cam);
    
if ~all(cam.aposteriori.rpy(:)==0) %estimator has been run!
    if cam.syntheticVideoFlag
        str_true_aposteriori = 'true';
        s = cam.true;
    else
        str_true_aposteriori = 'aposteriori';
        s = cam.aposteriori;
    end
    
    
    %FOCUS VECTORS
    zv = zeros(cam.frames,4);  s.corners.lat = zv;  s.corners.lon = zv;  s.corners.alt = zv;
    s.focus.lla = zeros(cam.frames,3);  xy = [xy; cam.true.focus.pixel];
    for i=vf
        lla = pixel2lla(cam.DEM,cam,i,xy,str_true_aposteriori);
        s.corners.lat(i,:) = lla(1:4,1)';
        s.corners.lon(i,:) = lla(1:4,2)';
        s.corners.alt(i,:) = lla(1:4,3)';
        s.focus.lla(i,:) = lla(5,:);
    end
    s.corners.plotlat = s.corners.lat(:,[1:4 1])';
    s.corners.plotlon = s.corners.lon(:,[1:4 1])';
    s.corners.plotalt = s.corners.alt(:,[1:4 1])';
    
    if cam.syntheticVideoFlag
        cam.true = s;
    else
        cam.aposteriori = s;
        s.lla = lla2llag(s.lla,EGM);
        s.focus.lla = lla2llag(s.focus.lla,EGM);
        s.corners.alt(1,:) = lla2ag([s.corners.lat(1,:)' s.corners.lon(1,:)' s.corners.alt(1,:)']);
    end

    
    %MSV RED POINTS
    if exist('msvned','var') && ~isempty(msvned)
        lla = ned2lla(cam.DEM,msvned);  if ~cam.syntheticVideoFlag;  lla=lla2llag(lla,EGM);  end %MSV SOLUTION
        cam.google.kml.tphatMSV = ge_folder('MSV Points',ge_point(lla(:,2),lla(:,1),lla(:,3),'extrude',1,'iconScale',.4,'iconURL','http://maps.google.com/mapfiles/kml/paddle/1-lv.png','altitudeMode','absolute'));
    end
    
    
    %NLS GREEN POINTS
    if exist('nlsned','var') && ~isempty(nlsned)
        lla = ned2lla(cam.DEM,nlsned);  if ~cam.syntheticVideoFlag;  lla=lla2llag(lla,EGM);  end %ILS SOLUTION
        cam.google.kml.tphatILS = ge_folder('NLS Points',ge_point(lla(:,2),lla(:,1),lla(:,3),'extrude',1,'iconScale',.4,'iconURL','http://maps.google.com/mapfiles/kml/paddle/grn-blank-lv.png','altitudeMode','absolute'));
    end
    
    
    %FOCUS POINT
    str1 = ge_plot(s.focus.lla(:,2),s.focus.lla(:,1),'lineWidth',2,'lineColor',[dec2hex(round(.7*255)) fcnhexcolor([.3 1 .3])],'name','Ground Focus Points');
    str2 = ge_plot3([s.lla(1,2); s.focus.lla(1,2)]',[s.lla(1,1); s.focus.lla(1,1)]',[s.lla(1,3); s.focus.lla(1,3)]','lineWidth',2,'lineColor',[dec2hex(round(.7*255)) fcnhexcolor([.3 1 .3])],'altitudeMode','absolute','name','Focus Vectors');
    cam.google.kml.focuspoint = ge_folder([str_true_aposteriori ' Focus Points'],[str1 str2]);
    
    
    %GROUND IMAGE OUTLINE
    cam.google.kml.imageoutlineonground = ge_plot(s.corners.plotlon,s.corners.plotlat,          'lineWidth',1,'lineColor',[dec2hex(round(.5*255)) fcnhexcolor([1 1 1])],'name','Image Projections on Ground');
    str1 =                                ge_plot(s.corners.plotlon(:,1),s.corners.plotlat(:,1),'lineWidth',2,'lineColor',[dec2hex(round(1*255))  fcnhexcolor([1 1 1])],'name','Image 1 Projection on Ground');

    
    %FRAME 1 IMAGE PYRAMID    
    ov=ones(1,4);  i=1;
    str2 = ge_plot3([s.lla(ov*i,2) s.corners.lon(i,:)']',[s.lla(ov*i,1) s.corners.lat(i,:)']',[s.lla(ov*i,3) s.corners.alt(i,:)']','lineWidth',2,'altitudeMode','absolute','name','Image 1 Projection Sides');%,'timeStamp',date1);
    cam.google.kml.frame1cameracone = ge_folder('Image 1 Projection',[str1 str2]);

    plottimeflag = 1;
    if plottimeflag
        str = [];
        for i=1:cam.frames
            t = datestr(cam.google.startdate+cam.true.t(i)*s2day, dateformat);
            
            %x = [[s.corners.lat(i,:)' s.corners.lon(i,:)' s.corners.alt(i,:)']; s.lla(i,:)];  j = [5 1 2 5 2 3 5 4 3 5 1 4 5];  x = x(j,:);
            %str1 = ge_plot3(x(:,2),x(:,1),x(:,3),'lineWidth',2,'altitudeMode','absolute','polyColor',[dec2hex(round(.3*255)) fcnhexcolor([1 1 1])],'timeStamp',date1);
            %str1 = ge_plot3([s.lla(ov*i,2) s.corners.lon(i,:)']',[s.lla(ov*i,1) s.corners.lat(i,:)']',[s.lla(ov*i,3) s.corners.alt(i,:)']','lineWidth',2,'altitudeMode','absolute','timeStamp',t);
            str2 = ge_plot3([s.lla(i,2); s.focus.lla(i,2)]',[s.lla(i,1); s.focus.lla(i,1)]',[s.lla(i,3); s.focus.lla(i,3)]','lineWidth',2,'lineColor',['FF' fcnhexcolor([1 1 1])],'altitudeMode','absolute','timeStamp',t);
            %str3 = ge_plot(s.corners.plotlon(:,i),s.corners.plotlat(:,i),'lineWidth',2,'lineColor',[dec2hex(round(1*255))  fcnhexcolor([1 1 1])],'timeStamp',t);
            str = [str str2];
        end
        cam.google.kml.timevaryingcameras = ge_folder([str_true_aposteriori ' LOS'],str);
    end
    
    if isfield(cam.aposteriori,'rpySFM')
        cam2=cam;  cam2.aposteriori.rpy=cam.aposteriori.rpySFM;  s.lla=cam.aposteriori.lla;

        str = [];
        for i=1:cam.frames
            s.focus.lla(i,:) = pixel2lla(cam.DEM,cam2,i,cam.true.focus.pixel,'aposteriori');
            if ~cam.syntheticVideoFlag;  s.focus.lla(i,:)=lla2llag(s.focus.lla(i,:),EGM);  end

            t = datestr(cam.google.startdate+cam.true.t(i)*s2day, dateformat);
            %x = [[s.corners.lat(i,:)' s.corners.lon(i,:)' s.corners.alt(i,:)']; s.lla(i,:)];  j = [5 1 2 5 2 3 5 4 3 5 1 4 5];  x = x(j,:);
            %str1 = ge_plot3(x(:,2),x(:,1),x(:,3),'lineWidth',2,'altitudeMode','absolute','polyColor',[dec2hex(round(.3*255)) fcnhexcolor([1 1 1])],'timeStamp',date1);
            %str1 = ge_plot3([s.lla(ov*i,2) s.corners.lon(i,:)']',[s.lla(ov*i,1) s.corners.lat(i,:)']',[s.lla(ov*i,3) s.corners.alt(i,:)']','lineWidth',2,'altitudeMode','absolute','timeStamp',t);
            str2 = ge_plot3([s.lla(i,2); s.focus.lla(i,2)]',[s.lla(i,1); s.focus.lla(i,1)]',[s.lla(i,3); s.focus.lla(i,3)]','lineWidth',2,'lineColor',['FF' fcnhexcolor([1 0 0])],'altitudeMode','absolute','timeStamp',t);
            %str3 = ge_plot(s.corners.plotlon(:,i),s.corners.plotlat(:,i),'lineWidth',2,'lineColor',[dec2hex(round(1*255))  fcnhexcolor([1 1 1])],'timeStamp',t);
            str = [str str2];
        end
        cam.google.kml.SFMlos = ge_folder('SFM LOS',str);
    end

    if isfield(cam.aposteriori,'rpyMSV')
        cam2=cam;  cam2.aposteriori.rpy=cam.aposteriori.rpyMSV;  s.lla=cam.aposteriori.lla;

        str = [];
        for i=1:cam.frames
            s.focus.lla(i,:) = pixel2lla(cam.DEM,cam2,i,cam.true.focus.pixel,'aposteriori');
            if ~cam.syntheticVideoFlag;  s.focus.lla(i,:)=lla2llag(s.focus.lla(i,:),EGM);  end
            
            t = datestr(cam.google.startdate+cam.true.t(i)*s2day, dateformat);
            %x = [[s.corners.lat(i,:)' s.corners.lon(i,:)' s.corners.alt(i,:)']; s.lla(i,:)];  j = [5 1 2 5 2 3 5 4 3 5 1 4 5];  x = x(j,:);
            %str1 = ge_plot3(x(:,2),x(:,1),x(:,3),'lineWidth',2,'altitudeMode','absolute','polyColor',[dec2hex(round(.3*255)) fcnhexcolor([1 1 1])],'timeStamp',date1);
            %str1 = ge_plot3([s.lla(ov*i,2) s.corners.lon(i,:)']',[s.lla(ov*i,1) s.corners.lat(i,:)']',[s.lla(ov*i,3) s.corners.alt(i,:)']','lineWidth',2,'altitudeMode','absolute','timeStamp',t);
            str2 = ge_plot3([s.lla(i,2); s.focus.lla(i,2)]',[s.lla(i,1); s.focus.lla(i,1)]',[s.lla(i,3); s.focus.lla(i,3)]','lineWidth',2,'lineColor',['FF' fcnhexcolor([0 0 1])],'altitudeMode','absolute','timeStamp',t);
            %str3 = ge_plot(s.corners.plotlon(:,i),s.corners.plotlat(:,i),'lineWidth',2,'lineColor',[dec2hex(round(1*255))  fcnhexcolor([1 1 1])],'timeStamp',t);
            str = [str str2];
        end
        cam.google.kml.MSVlos = ge_folder('MSV LOS',str);
    end
    
end


%SAVE TO KML FILE ---------------------------------------------------------
fields = fieldnames(cam.google.kml);  str=[];
for i=1:numel(fields);  str = [str ' cam.google.kml.' fields{i}];  end;  str = ['[' str ']']; %#ok<*AGROW>
ge_output([cam.pathname cam.filename '.kml'],eval(str));
end

function str = aircraftpath(cam,vf)
if cam.syntheticVideoFlag;
    lla=cam.true.lla; 
else
    lla=lla2llag(cam.apriori.lla);
end;  
naf = numel(lla(:,1));


str = ge_plot3(lla(vf,2),lla(vf,1),lla(vf,3),'extrude',1,'lineWidth',1,'altitudeMode','absolute','name','Aircraft GPS');
if naf~=cam.frames
    str = [str, ge_plot3(lla(:,2),lla(:,1),lla(:,3),'extrude',1,'lineWidth',1,'altitudeMode','absolute','lineColor',[dec2hex(round(.25*255)) fcnhexcolor([1 1 1])])];
end
end

function str = DEM(cam)
n = numel(cam.DEM.lat(:,1)); ov = ones(n,1);
x = cam.DEM.lng(:,1);
y = cam.DEM.lat(1,:)';

x1 = x([1 end 1],ov);
y1 = y(:,ones(3,1))';
y2 = y([1 end 1],ov);
x2 = x(:,ones(3,1))';
lng = [x1 x2];  lat = [y1 y2];
str = ge_plot(lng,lat,'tessellate',1,'lineWidth',1,'lineColor',[dec2hex(round(.25*255)) fcnhexcolor([.7 .7 .7])],'name','GE DEM');
end