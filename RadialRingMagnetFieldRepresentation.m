function RadialRingMagnetFieldRepresentation

% Currently De-bugging as of 10/09/18
% - Elliptic integrals FStarCalc and PiStarCalc working properly--verified
% through Wolfram test functions.
%
% - Hz component works, not Hr yet.
%
% - Beta function seems to be correct--have reviewed it multiple times

clear
close all
clc

% Based on Ravaud 2008

%% Constants
mu0 = 4*pi*1e-7;
sigmaStar = 1;      % "surface magnetic pole density"
rin = 0.025;         % Inner radius of magnet ring
rout = 0.028;        % Outer radius of magnet ring
h = 0.003;           % Depth of magnet ring

u2 = 0.999;
u1 = -u2;

%% Looping variables
r = 0.0255;
z = 0.0015;

rValues = [0.001,0.0045,0.009,0.0135,0.018,0.0225,0.023,0.0235,0.024,0.0249,0.027,0.0315];

iterVar = 0;
for r = rValues
iterVar = iterVar + 1;
%% Azimuthal component
Htheta = 0;     % Symmetry of the ring magnet leads to no theta component

%% Axial component

% Hz = sigmaStar/(4*pi*mu0) *
% (-4*rin*KStar[4*r*rin/(r^2+rin^2-2*r*rin+z^2)]/sqrt(r^2+rin^2-2*r*rin+z^2)
% +4*rin*KStar[4*r*rin/(r^2+rin^2-2*r*rin+(z-h)^2)]/sqrt(r^2+rin^2-2*r*rin+(z-h)^2)

Hz(iterVar) = calcHz(r,z);

%% Radial component

%

Hr(iterVar) = calcHr(r,z);


end
figure(1)
plot(rValues,Hz,'bo-')
title('Hz')
figure(2)
plot(rValues,Hr,'bo-')
title('Hr')

%% FUNCTIONS

% This is the function according to the paper UNCHANGED
%     function Hz = calcHz(r,z)
%         Hz = sigmaStar/(4*pi*mu0)*(-4*rin*KStarCalc(4*r*rin/(r^2+rin^2-2*r*rin+z^2))/...
%             sqrt(r^2+rin^2-2*r*rin+z^2)+(4*rin*KStarCalc(4*r*rin/(r^2+rin^2-2*r*rin+(z-h)^2)))...
%             /sqrt(r^2+rin^2-2*r*rin+(z-h)^2)); % (21) from Ravaud 2008
%     end
% DETERMINED THAT THERE IS AN EXTRA MINUS SIGN IN THE PAPER IN KSTAR CALCS

% This is the function according to the paper with changes that correct
% what I suspect are mistakes made by the authors... namely, 1 minus
% sign inside the KStarCalc rather than 2 which cancel.
    function Hz = calcHz(r,z)
        Hz = sigmaStar/(4*pi*mu0)*(-4*rin*KStarCalc(-4*r*rin/(r^2+rin^2-2*r*rin+z^2))/...
            sqrt(r^2+rin^2-2*r*rin+z^2)+(4*rin*KStarCalc(-4*r*rin/(r^2+rin^2-2*r*rin+(z-h)^2)))...
            /sqrt(r^2+rin^2-2*r*rin+(z-h)^2)); % (21) from Ravaud 2008
    end



    function KStar = KStarCalc(m)
        KStar = FStarCalc(pi/2,m);  % (22) from Ravaud 2008
    end
        
        
        

    function FStar = FStarCalc(phi,m)
        ellipticIntegral = @(theta) 1./sqrt(1-m*sin(theta).^2);
        phi = real(phi);
        FStar = integral(ellipticIntegral,0,phi);   % (23) from Ravud 2008
    end




    function Hr = calcHr(r,z)
        Hr = sigmaStar/(2*pi*mu0)*(betaCalc(r,z,u1)-betaCalc(r,z,u2));
    end




    function beta = betaCalc(r,z,u)
        a1 = rin*r*z;
        b1 = -rin^2*z;
        c = r^2+rin^2;
        d = -2*r*rin;
        e1 = z^2;
        a2 = -rin*r*(z-h);
        b2 = rin^2*(z-h);
        e2 = (z-h)^2;
        beta = 2*1i * (1+u) * sqrt(d*(-1+u)/(c+e1+d*u)) * (-(a1*d+b1*(c+e1))) *...
            FStarCalc(1i*asinh(sqrt(-c+d-e1)/sqrt(c+e1+d*u)),(c+d+e1)/(c-d+e1))/...
            (d*sqrt(-c+d-e1)*e1*sqrt(d*(1+u)/(c+e1+d*u))*sqrt(1-u^2))+...
            0+...
            2*1i * (1+u) * sqrt(d*(-1+u)/(c+e1+d*u)) * (b1*c-a1*d) * PiStarCalc(...
            e1/(c-d+e1),1i*asinh(sqrt(-c+d+e1)/(c+e1+d*u)),(c+d+e1)/(c-d+e1))/...
            (d*sqrt(-c+d-e1)*e1*sqrt(d*(1+u)/(c+e1+d*u))*sqrt(1-u^2))+...
            0+...
            2*1i * (1+u) * sqrt(d*(-1+u)/(c+e2+d*u)) * (-(a2*d+b2*(c+e2))) * FStarCalc(...
            1i*asinh(sqrt(-c+d-e2)/sqrt(c+e2+d*u)),(c+d+e2)/(c-d+e2))/...
            (d*sqrt(-c+d-e2)*e2*sqrt(d*(1+u)/(c+e2+d*u))*sqrt(1-u^2))+...
            0+...
            2*1i * (1+u) * sqrt(d*(-1+u)/(c+e2+d*u)) * (b2*c-a2*d) * PiStarCalc(...
            e2/(c-d+e2),1i*asinh(sqrt(-c+d+e2)/(c+e2+d*u)),(c+d+e2)/(c-d+e2))/...
            (d*sqrt(-c+d-e2)*e2*sqrt(d*(1+u)/(c+e2+d*u))*sqrt(1-u^2));
        beta = real(beta);
        
    end




    function PiStar = PiStarCalc(n,phi,m)
        ellipticIntegral = @(theta) 1./((1-n*sin(theta).^2).*sqrt(1-m*sin(theta).^2));
        PiStar = integral(ellipticIntegral,0,phi);
    end









end




















