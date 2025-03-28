% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function E = fcnF2E(F,K)
%Fundamental Matrix F to Essential Matrix E
E = K'*F*K;

[U,S,V] = svd(E);
m = S(1,1)/2 + S(2,2)/2;
S = diag([m m 0]);
E = U*S*V';
%[Up,Sp,Vp] = svd(E);
end

