% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

% getCorrectCameraMatrix - Check which of the 4 solutions from the
% essential matrix is the correct one.
%
%
% Given the Essential matrix and its decomposition into 4 possible solutions
% the appropriate solution is found using a voting mechanism in which each
% point correspondence votes for one solution based on the position of the
% reconstructed point, which must be in front of both cameras. Camera one
% is assumed to be [I|0].
%
%
% Input  - PXcam  -> (3x4xn) Possible camera solutions
%        - E      -> (3x3xn) Essential matrix that generated the Pcam
%        - K1     -> (3x3) Camera calibration of image 1
%        - K2     -> (3x3) Camera calibration of image 2
%        - X      -> (3*2xn) homogeneous points in images 2 and 1
%
% Output - P      -> (3x4) Correct camera matrix (rotation and translation)
%        - voting -> (4x1) Votes for each solution
%
%
%
% Author: Isaac Esteban
% IAS, University of Amsterdam
% TNO Defense, Security and Safety
% isaac@fit3d.info
% isaac.esteban@tno.nl
% Copyright TNO - 2010

function [P,voting] = getCorrectCameraMatrix(P2_4,K1,K2,X1,X2)
np=size(X1,2);

% The first camera matrix is taken P = [I|0] and the other
P1 = [eye(3) zeros(3,1)];
P = K1*P1;

% For every correspondence 
Depth1=zeros(4,np);  Depth2=Depth1;
for i=1:np
    
    % Two matching points in image coordinates (x in image 1 and xp in image 2)
    x1 = X1(1:3,i);
    x2 = X2(1:3,i);
    
    % The first camera matrix is taken P = [I|0] and the other
    xhat1 = K1\x1;
    xhat2 = K2\x2;
    
    A1n=sqrt(sum(xhat1(1).^2+1));  A2n=sqrt(sum(xhat1(2).^2+1)); %for normalizing
    A1v = [-1    0  xhat1(1)    0]/A1n;
    A2v = [ 0   -1  xhat1(2)    0]/A2n;
    
    % For each camera matrix (PXcam), reproject the pair of points in 3D and determine the depth in 3D of the point
    for p=1:4
        P2 = P2_4(:,:,p);
        
        % Build and Normalize A
        A3v = P2(3,:)*xhat2(1)-P2(1,:);
        A4v = P2(3,:)*xhat2(2)-P2(2,:);
        A = [A1v;  A2v;  A3v;  A4v];
        %Anorm = [A1v;  A2v;  A3v/sqrt(sum(A3v.^2));  A4v/sqrt(sum(A4v.^2))];
        
        % Obtain the 3D point
        [~,~,V]=svd(A);
        X3 = V(:,4);
        T = X3(4);
        
        % Check depth on first camera
        w = X3(3);
        Depth1(p,i) = w/T;
        
        % Check depth on second camera
        xi = P2*X3;
        w = xi(3);
        m3n = sqrt(sum(P2(3,1:3).^2));
        Depth2(p,i) = (sign(det(P2(:,1:3)))*w)./(T*m3n);
    end
end
[~,i] = max(sum(Depth1>0 & Depth2>0,2));
P = P2_4(:,:,i);


% function [P,voting] = getCorrectCameraMatrix(P2_4,K1,K2,X1,X2)
% np = 3; %size(X1,2);
% X3_4 = zeros(4,np,4);
% 
% % The first camera matrix is taken P = [I|0] and the other
% P1 = [eye(3) zeros(3,1)];
% P = K1*P1;
% 
% % For every correspondence
% for i=1:3 %1:np
%     
%     % Two matching points in image coordinates (x in image 1 and xp in image 2)
%     x1 = X1(1:3,i);
%     x2 = X2(1:3,i);
%     
%     % The first camera matrix is taken P = [I|0] and the other
%     xhat1 = K1\x1;
%     xhat2 = K2\x2;
%     
%     A1n=sqrt(sum(xhat1(1).^2+1));  A2n=sqrt(sum(xhat1(2).^2+1)); %for normalizing
%     A1v = [-1    0  xhat1(1)    0]/A1n;
%     A2v = [ 0   -1  xhat1(2)    0]/A2n;
%     
%     % For each camera matrix (PXcam), reproject the pair of points in 3D and determine the depth in 3D of the point
%     for p=1:4
%         P2 = P2_4(:,:,p);
%         
%         % Build and Normalize A
%         A3v = P2(3,:)*xhat2(1)-P2(1,:);
%         A4v = P2(3,:)*xhat2(2)-P2(2,:);
%         A = [A1v;  A2v;  A3v;  A4v];
%         %Anorm = [A1v;  A2v;  A3v/sqrt(sum(A3v.^2));  A4v/sqrt(sum(A4v.^2))];
%         
%         % Obtain the 3D point
%         %[~,~,V]=svd(Anorm);  X3=V(:,end);
%         [eigvecs,~]=eig(A); X3=eigvecs(:,4);
%         X3_4(:,i,p) = X3;
%     end
% end
% 
% Depth1=zeros(4,np);  Depth2=Depth1;
% for p=1:4    
%     P2 = P2_4(:,:,p);
%     X3 = X3_4(:,:,p);
%     T = X3(end,:);
%     
%     % Check depth on first camera
%     xi = P1*X3;
%     w = xi(3,:);
%     m3n = sqrt(sum(P1(3,1:3).^2));
%     Depth1(p,:) = (sign(det(P1(:,1:3)))*w)./(T*m3n);
%     
%     % Check depth on second camera
%     xi = P2*X3;
%     w = xi(3,:);
%     m3n = sqrt(sum(P2(3,1:3).^2));
%     Depth2(p,:) = (sign(det(P2(:,1:3)))*w)./(T*m3n);
% end
% 
% voting = sum(Depth1>0 & Depth2>0,2);
% [~, i] = max(voting);
% P = P2_4(:,:,i);


% function [P,voting] = getCorrectCameraMatrix(P2_4,K1,K2,X1,X2)
% np = 3; %size(X1,2);
% X3_4 = zeros(4,np,4);
% 
% % The first camera matrix is taken P = [I|0] and the other
% P1 = [eye(3) zeros(3,1)];
% P = K1*P1;
% 
% % For every correspondence
% for i=1:3 %1:np
%     
%     % Two matching points in image coordinates (x in image 1 and xp in image 2)
%     x1 = X1(1:3,i);
%     x2 = X2(1:3,i);
%     
%     % The first camera matrix is taken P = [I|0] and the other
%     xhat1 = K1\x1;
%     xhat2 = K2\x2;
%     
%     % For each camera matrix (PXcam), reproject the pair of points in 3D and determine the depth in 3D of the point
%     Depth = zeros(4,2);
%     for p=1:4
%         P2 = P2_4(:,:,p);
%         
%         % We build the matrix A
%         A = [P1(3,:)*xhat1(1)-P1(1,:)
%             P1(3,:)*xhat1(2)-P1(2,:)
%             P2(3,:)*xhat2(1)-P2(1,:)
%             P2(3,:)*xhat2(2)-P2(2,:)];
%         
%         % Normalize A
%         A1n=sqrt(sum(A(1,:).^2));  A2n=sqrt(sum(A(2,:).^2));  A3n=sqrt(sum(A(3,:).^2));  A4n=sqrt(sum(A(4,:).^2));
%         Anorm = [A(1,:)/A1n;  A(2,:)/A2n;  A(3,:)/A3n;  A(4,:)/A4n];
%         
%         % Obtain the 3D point
%         [U,S,V] = svd(Anorm);
%         X3 = V(:,end);
%         X3_4(:,i,p) = X3;
%     end
% end
% 
% Depth1=zeros(4,np);  Depth2=Depth1;
% for p=1:4    
%     P2 = P2_4(:,:,p);
%     X3 = X3_4(:,:,p);
%     T = X3(end,:);
%     
%     % Check depth on first camera
%     xi = P1*X3;
%     w = xi(3,:);
%     m3n = sqrt(sum(P1(3,1:3).^2));
%     Depth1(p,:) = (sign(det(P1(:,1:3)))*w)./(T*m3n);
%     
%     % Check depth on second camera
%     xi = P2*X3;
%     w = xi(3,:);
%     m3n = sqrt(sum(P2(3,1:3).^2));
%     Depth2(p,:) = (sign(det(P2(:,1:3)))*w)./(T*m3n);
% end
% 
% voting = sum(Depth1>0 & Depth2>0,2);
% [~, i] = max(voting);
% P = P2_4(:,:,i);