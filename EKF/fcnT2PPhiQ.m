function ekf = fcnT2PPhiQ(cam,ekf,dt)
m2km = 1/1000;

%position and velocity ----------------------------------------------------
B1a = 1/108; %seconds
B2a = 1/162; %seconds
[s2a, P0a] = KF2_simple_baseline_april_10_2012_init_cov(1/B1a, 1/B2a, cam.ilsxy1s*m2km); %meters
%s2a = .18*m2km; %km/s^2
s1a = s2a/2; %km/s
[Phi61, Q61] = fcngetPhiQ(s1a,s2a,B1a,B2a,dt);

%vertical position and velocity -------------------------------------------
[s2az, P0az] = KF2_simple_baseline_april_10_2012_init_cov(1/B1a, 1/B2a, cam.ilsz1s*m2km); %meters
s1az = s2az/2; %km/s
[Phi61z, Q61z] = fcngetPhiQ(s1az,s2az,B1a,B2a,dt);
Phi61(5:6,5:6) = Phi61z(5:6,5:6);
Q61(5:6,5:6) = Q61z(5:6,5:6);

%rotational position and velocity -----------------------------------------
B1b = 1/36*100; %seconds
B2b = 1/54*100; %seconds
% B1b = 1/(cam.ilsTrpy*.4); %seconds
% B2b = 1/(cam.ilsTrpy*.6); %seconds
[s2b, P0b] = KF2_simple_baseline_april_10_2012_init_cov(1/B1b, 1/B2b, cam.ilsrpy1s); %deg
%s2b = 1/100; %deg/s^2
s1b = s2b/2; %deg/s
[Phi62, Q62] = fcngetPhiQ(s1b,s2b,B1b,B2b,dt);

%combine ------------------------------------------------------------------
zm2=zeros(2);
zm6=zeros(6);
Phi = [ Phi61     zm6
        zm6       Phi62  ];
Q =   [ Q61       zm6
        zm6       Q62    ];

ekf.Phi12 = Phi;
ekf.Q12 = Q;
ekf.P0_12 = [P0a     zm2     zm2     zm2     zm2     zm2
             zm2     P0a     zm2     zm2     zm2     zm2
             zm2     zm2     P0az    zm2     zm2     zm2
             zm2     zm2     zm2     P0b     zm2     zm2
             zm2     zm2     zm2     zm2     P0b     zm2
             zm2     zm2     zm2     zm2     zm2     P0b];
end


function [Phi61, Q61] = fcngetPhiQ(s1,s2,B1,B2,dt)
zm = zeros(2);
Phi2 = [ exp(-B1*dt)                      (B2-B1)^-1*(exp(-B1*dt)-exp(-B2*dt))
        0                                 exp(-B2*dt)                                 ];
Phi61 = [ Phi2    zm      zm
          zm      Phi2    zm
          zm      zm      Phi2 ];

Q2 = zm;
Q2(1,1) = s2^2*(B2-B1)^-2*(.5*B1^-1*(1-exp(-2*B1*dt))+.5*B2^-1*(1-exp(-2*B2*dt))-2*(B1+B2)^-1*(1-exp(-(B1+B2)*dt)))  +  s1^2*.5*B1^-1*(1-exp(-2*B1*dt));
Q2(1,2) = s2^2*(B2-B1)^-1*((B1+B2)^-1*(1-exp(-(B1+B2)*dt)) - .5*B2^-1*(1-exp(-2*B2*dt)));
Q2(2,1) = Q2(1,2);
Q2(2,2) = s2^2*.5*B2^-1*(1-exp(-2*B2*dt));
Q61 = [ Q2    zm    zm
        zm    Q2    zm 
        zm    zm    Q2 ];  
end


