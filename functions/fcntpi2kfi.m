% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function kfi = fcntpi2kfi(tpi)
%tie point index to kalman filter index
kfi = ((tpi-1)*3+1:tpi*3)+12;


end

