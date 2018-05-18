function P = fcnE2P(E,method)
P = zeros(3,4,4);
%[eigvec,eigval]=eigs(E)
if nargin==1;  method = 'Horn';  end

switch method
    case 'Zisserman'
        % ZISSERMAN METHOD ---------------------------------------------------------
        [U,S,V] = svd(E);
        
        % Matrix W
        W = [0  -1   0
             1   0   0
             0   0   1];
        
        % Compute 4 possible solutions (p259) and scale the translation vectors
        R1 = U*W*V';
        R2 = U*W'*V';
        t = U(:,3)/max(abs(U(:,3)));

    case 'Horn'
        %HORN METHOD --------------------------------------------------------------
        % based on Horn (1990) paper. This method is apparently more stable than the one presented by Zisserman.
        
        % Obtain bxbT
        EEp = E*E';
        B = .5*( EEp(1,1) + EEp(2,2) + EEp(3,3) )*eye(3,3) - EEp;
        
        % We obtain b by selecting the largest row and dividing by the square root of the diagonal. This is for numerical accuracy.
        if B(1,1)>B(2,2) && B(1,1)>B(3,3)
            t = B(1,:)'/B(1,1).^.5;
        elseif B(2,2)>B(1,1) && B(2,2)>B(3,3)
            t = B(2,:)'/B(2,2).^.5;
        else
            t = B(3,:)'/B(3,3).^.5;
        end
        rt=sum(t.*t).^.5;
        %USE A DIFFERENT TRANSLATION----------------
        %t = fcnvec2uvec(t0)'*rt;
        
        % Compute cofactors of E
        cofE = [fcncross1(E(:,2),E(:,3)), fcncross1(E(:,3),E(:,1)), fcncross1(E(:,1),E(:,2))];
        
        % Compute two matricex B
        BE = [   0    -t(3)    t(2)
            t(3)      0    -t(1)
            -t(2)    t(1)      0    ] * E;
        b1dot = sum(t.*t); %dot(b1,b1)
        
        % Compute both R and scale baseline
        R1 = (cofE-BE)/b1dot;
        R2 = (cofE+BE)/b1dot;
        t = t/max(abs(t));
end

%BUILD 4 POSSIBLE SOLUTIONS -----------------------------------------------
P(:,:,1) = [R1,  t];
P(:,:,2) = [R1, -t];
P(:,:,3) = [R2,  t];
P(:,:,4) = [R2, -t];
end