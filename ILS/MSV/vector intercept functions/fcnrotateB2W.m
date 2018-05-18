function varargout = fcnrotateB2W(varargin)

switch nargin
    case 2 %(rpy,x)
        rpy = varargin{1};
        r = rpy(:,1);  
        p = rpy(:,2);  
        yaw = rpy(:,3);
        x = varargin{2}; 
    case 4 %(r,p,yaw,x)
        r = varargin{1};
        p = varargin{2};
        yaw = varargin{3};
        x = varargin{4};
end


sr=sin(r); sp=sin(p); sy=sin(yaw);
cr=cos(r); cp=cos(p); cy=cos(yaw);

switch nargout
    case 1 %(x2)
        x2 = zeros(size(x));
        x2(:,1)=x(:,1).*(cp.*cy)+x(:,2).*(cp.*sy )+x(:,3).*(-sp);
        x2(:,2)=x(:,1).*(sr.*sp.*cy-cr.*sy)+x(:,2).*(sr.*sp.*sy+cr.*cy)+x(:,3).*(sr.*cp);
        x2(:,3)=x(:,1).*(cr.*sp.*cy+sr.*sy)+x(:,2).*(cr.*sp.*sy-sr.*cy)+x(:,3).*(cr.*cp);
        varargout{1} = x2;
    case 3 %(x2,y2,z2)
        %varargout = cell(1,3);
        varargout{1}=x(:,1).*(cp.*cy)+x(:,2).*(cp.*sy )+x(:,3).*(-sp);
        varargout{2}=x(:,1).*(sr.*sp.*cy-cr.*sy)+x(:,2).*(sr.*sp.*sy+cr.*cy)+x(:,3).*(sr.*cp);
        varargout{3}=x(:,1).*(cr.*sp.*cy+sr.*sy)+x(:,2).*(cr.*sp.*sy-sr.*cy)+x(:,3).*(cr.*cp);
end





% function x2 = fcnrotateB2W(r,p,y,x)
% sr=sin(r); sp=sin(p); sy=sin(y);
% cr=cos(r); cp=cos(p); cy=cos(y);
% x2 = zeros(size(x));
% x2(:,1)=x(:,1).*(cp.*cy)+x(:,2).*(cp.*sy )+x(:,3).*(-sp);
% x2(:,2)=x(:,1).*(sr.*sp.*cy-cr.*sy)+x(:,2).*(sr.*sp.*sy+cr.*cy)+x(:,3).*(sr.*cp);
% x2(:,3)=x(:,1).*(cr.*sp.*cy+sr.*sy)+x(:,2).*(cr.*sp.*sy-sr.*cy)+x(:,3).*(cr.*cp);
% end