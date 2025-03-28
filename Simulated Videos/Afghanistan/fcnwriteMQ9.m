% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [] = fcnwriteMQ9(lla, rpy, heading, tilt, range)
daename  = [cd filesep 'GEfiles' filesep 'MQ9.dae'];
kmlfname = [cd filesep 'GEfiles' filesep 'MQ9.kml'];
fid = fopen(kmlfname,'w');

fprintf(fid,'<?xml version="1.0" encoding="UTF-8"?>\n');
fprintf(fid,'<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">\n');
fprintf(fid,'<Document>\n');
fprintf(fid,'<name>MQ9.kml</name>\n');
fprintf(fid,'<Style id="sn_ylw-pushpin">\n');
fprintf(fid,'</Style>\n');
fprintf(fid,'<Placemark>\n');
fprintf(fid,'  <name>MQ9 Name</name>\n');
fprintf(fid,'	<LookAt>\n');
fprintf(fid,'       <altitudeMode>relativeToGround</altitudeMode>\n');
fprintf(fid,'		<longitude>%.12f</longitude>\n',lla(2));
fprintf(fid,'		<latitude>%.12f</latitude>\n',lla(1));
fprintf(fid,'		<altitude>%.12f</altitude>\n',lla(3));
fprintf(fid,'		<heading>%.12f</heading>\n',heading);
fprintf(fid,'		<tilt>%.12f</tilt>\n',tilt);
fprintf(fid,'		<range>%.12f</range>\n',range);
fprintf(fid,'   </LookAt>\n');
fprintf(fid,'   <styleUrl>#sn_ylw-pushpin</styleUrl>\n');
fprintf(fid,'   <Model id="model_4">\n');
fprintf(fid,'   <altitudeMode>relativeToGround</altitudeMode>\n');
fprintf(fid,'   <Location>\n');
fprintf(fid,'       <longitude>%.12f</longitude>\n',lla(2)); %deg
fprintf(fid,'       <latitude>%.12f</latitude>\n',lla(1)); %deg
fprintf(fid,'       <altitude>%.12f</altitude>\n',lla(3)); %m
fprintf(fid,'   </Location>\n');
fprintf(fid,'   <Orientation>\n');
fprintf(fid,'       <heading>%.12f</heading>\n',rpy(3)); %deg
fprintf(fid,'       <tilt>%.12f</tilt>\n',rpy(2)); %deg
fprintf(fid,'       <roll>%.12f</roll>\n',rpy(1)); %deg
fprintf(fid,'   </Orientation>\n');
fprintf(fid,'   <Scale>\n');
fprintf(fid,'      <x>1</x>\n');
fprintf(fid,'      <y>1</y>\n');
fprintf(fid,'      <z>1</z>\n');
fprintf(fid,'   </Scale>\n');
fprintf(fid,'   <Link>\n');
fprintf(fid,'       <href>%s</href>\n',daename);
fprintf(fid,'   </Link>\n');
fprintf(fid,'   <ResourceMap>\n');
fprintf(fid,'   </ResourceMap>\n');
fprintf(fid,'   </Model>\n');
fprintf(fid,'</Placemark>\n');
fprintf(fid,'</Document>\n');
fprintf(fid,'</kml>\n');

fclose(fid);
end

