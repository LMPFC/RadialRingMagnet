function RadialRingMagnetFieldRepresentation_101218

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

% Based on Ravaud 2008, Babic 2008

%% Constants
mu0 = 4*pi*1e-7;
sigmaStar = 1;      % "surface magnetic pole density"
rin = 0.025;         % Inner radius of magnet ring
rout = 0.035;        % Outer radius of magnet ring
h = 0.01;           % Depth of magnet ring
obsDist = 0.02;
dSpace = 0.00031;

u2 = 0.999999999;
u1 = -u2;

% Original paper settings for magnet:
% rin = 0.025;
% rout = 0.028;
% h = 0.003;
% obsDist = 0.003;
% dSpace = 0.0003;

%% Looping variables

% rValues = 0.022:0.0003:0.031;
% zValues = -0.003:0.0003:0.006;
rValues = rin-obsDist:dSpace:rout+obsDist;
zValues = -obsDist:dSpace:h+obsDist;

rSize = length(rValues);
zSize = length(zValues);

rSpace = ones(zSize,1)*rValues;
zSpace = (ones(rSize,1)*zValues)';

HzNew = zeros(zSize,rSize);
HrNew = zeros(zSize,rSize);
HthetaNew = zeros(zSize,rSize);

rIterVar = 0;
for r = rValues
    rIterVar = rIterVar + 1;
    zIterVar = 0;
    for z = zValues
        zIterVar = zIterVar + 1;
        if r < rin-0.0001 || r > rout+0.0001 || z < -0.0001 || z > h+0.0001

        %% Azimuthal component
        Htheta = 0;     % Symmetry of the ring magnet leads to no theta component
        HthetaNew = 0;
        %% Axial component
        
        % Hz = sigmaStar/(4*pi*mu0) *
        % (-4*rin*KStar[4*r*rin/(r^2+rin^2-2*r*rin+z^2)]/sqrt(r^2+rin^2-2*r*rin+z^2)
        % +4*rin*KStar[4*r*rin/(r^2+rin^2-2*r*rin+(z-h)^2)]/sqrt(r^2+rin^2-2*r*rin+(z-h)^2)
        
        % Hz(zIterVar,rIterVar) = calcHz(r,z);
        HzNew(zIterVar,rIterVar) = calcHzNew(r,z);  % Using method from: Babic 2008
        
        %% Radial component
        
        %
        
        % Hr(zIterVar,rIterVar) = calcHr(r,z);      % Using method from: Ravaud 2008
        HrNew(zIterVar,rIterVar) = calcHrNew(r,z);  % Using method from: Babic 2008
        
        end
    end
end

quiver(rSpace,zSpace,HrNew,HzNew)
% figure(1)
% plot(rValues,Hz,'bo-')
% title('Hz')
% figure(2)
% plot(rValues,Hr,'bo-')
% title('Hr')
% figure(3)
% plot(rValues,HzNew,'bo-')
% title('Hz_{New}')
% figure(4)
% plot(rValues,HrNew,'bo-')
% title('Hr_{New}')



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
            2*1i * (1+u) * sqrt(d*(-1+u)/(c+e1+d*u)) * (b1*c-a1*d) * PiStarCalc(...
            e1/(c-d+e1),1i*asinh(sqrt(-c+d+e1)/(c+e1+d*u)),(c+d+e1)/(c-d+e1))/...
            (d*sqrt(-c+d-e1)*e1*sqrt(d*(1+u)/(c+e1+d*u))*sqrt(1-u^2))+...
            2*1i * (1+u) * sqrt(d*(-1+u)/(c+e2+d*u)) * (-(a2*d+b2*(c+e2))) * FStarCalc(...
            1i*asinh(sqrt(-c+d-e2)/sqrt(c+e2+d*u)),(c+d+e2)/(c-d+e2))/...
            (d*sqrt(-c+d-e2)*e2*sqrt(d*(1+u)/(c+e2+d*u))*sqrt(1-u^2))+...
            2*1i * (1+u) * sqrt(d*(-1+u)/(c+e2+d*u)) * (b2*c-a2*d) * PiStarCalc(...
            e2/(c-d+e2),1i*asinh(sqrt(-c+d+e2)/(c+e2+d*u)),(c+d+e2)/(c-d+e2))/...
            (d*sqrt(-c+d-e2)*e2*sqrt(d*(1+u)/(c+e2+d*u))*sqrt(1-u^2));
    end
    function PiStar = PiStarCalc(n,phi,m)
        ellipticIntegral = @(theta) 1./((1-n*sin(theta).^2).*sqrt(1-m*sin(theta).^2));
        PiStar = integral(ellipticIntegral,0,phi);
    end
    function HzNew = calcHzNew(r,z)     % Relatively confident that this works.
        HzNewPlus = sigmaStar/(2*pi*mu0)*(HzPlusSumCalc(r,z,1)+HzPlusSumCalc(r,z,2));
        HzNewMinus = -sigmaStar/(2*pi*mu0)*(HzMinusSumCalc(r,z,1)+HzMinusSumCalc(r,z,2));
        HzNew = HzNewPlus + HzNewMinus;
    end
    function HzPlusComp = HzPlusSumCalc(r,z,n)
        t(1) = z-h;
        t(2) = z;
        kPlus(1) = sqrt(4*r*rin/((r+rin)^2+t(1)^2));
        kPlus(2) = sqrt(4*r*rin/((r+rin)^2+t(2)^2));
        HzPlusComp = (-1)^(n-1)*kPlus(n)*sqrt(rin/r)*KStarCalc(kPlus(n));
    end
    function HzMinusComp = HzMinusSumCalc(r,z,n)
        t(1) = z-h;
        t(2) = z;
        kMinus(1) = sqrt(4*r*rout/((r+rout)^2+t(1)^2));
        kMinus(2) = sqrt(4*r*rout/((r+rout)^2+t(2)^2));
        HzMinusComp = (-1)^(n-1)*kMinus(n)*sqrt(rout/r)*KStarCalc(kMinus(n));
    end
    function HrNew = calcHrNew(r,z)
        HrNewPlus = -sigmaStar/(4*pi*mu0)*(HrPlusSumCalc(r,z,1)+HrPlusSumCalc(r,z,2));
        HrNewMinus = sigmaStar/(4*pi*mu0)*(HrMinusSumCalc(r,z,1)+HrMinusSumCalc(r,z,2));
        HrNew = HrNewPlus + HrNewMinus;
    end
    function HrPlusComp = HrPlusSumCalc(r,z,n)
        t(1) = z-h;
        t(2) = z;
        kPlus(1) = sqrt(4*r*rin/((r+rin)^2+t(1)^2));
        kPlus(2) = sqrt(4*r*rin/((r+rin)^2+t(2)^2));
        kMinus(1) = sqrt(4*r*rout/((r+rout)^2+t(1)^2));
        kMinus(2) = sqrt(4*r*rout/((r+rout)^2+t(2)^2));
        hPlus = 4*r*rin/(r+rin)^2;
        hMinus = 4*r*rout/(r+rout)^2;
        epsPlus(1) = asin(sqrt((1-hPlus)/(1-kPlus(1)^2)));
        epsPlus(2) = asin(sqrt((1-hPlus)/(1-kPlus(2)^2)));
        epsMinus(1) = asin(sqrt((1-hMinus)/(1-kMinus(2)^2)));
        epsMinus(2) = asin(sqrt((1-hMinus)/(1-kMinus(2)^2)));
        
        % "Suitable for regular or singular cases" --doesn't work
%         HrPlusComp = (-1)^(n-1)*sqrt(rin/r)*(2*t(n)*kPlus(n)/(r+rin)*KStarCalc(kPlus(n))+...
%             pi*sqrt(rin/r)*sign(r-rin)*sign(t(n))*(1-Lambda0Calc(epsPlus(n),kPlus(n))));
        
        
        % Regular Hr
        HrPlusComp = (-1)^(n-1)*t(n)*kPlus(n)/r*sqrt(rin/r)*(KStarCalc(kPlus(n))+...
            (r-rin)/(r+rin)*PiStarCalc(hPlus,pi/2,kPlus(n)));
    end
    function HrMinusComp = HrMinusSumCalc(r,z,n)
        t(1) = z-h;
        t(2) = z;
        kPlus(1) = sqrt(4*r*rin/((r+rin)^2+t(1)^2));
        kPlus(2) = sqrt(4*r*rin/((r+rin)^2+t(2)^2));
        kMinus(1) = sqrt(4*r*rout/((r+rout)^2+t(1)^2));
        kMinus(2) = sqrt(4*r*rout/((r+rout)^2+t(2)^2));
        hPlus = 4*r*rin/(r+rin)^2;
        hMinus = 4*r*rout/(r+rout)^2;
        epsPlus(1) = asin(sqrt((1-hPlus)/(1-kPlus(1)^2)));
        epsPlus(2) = asin(sqrt((1-hPlus)/(1-kPlus(2)^2)));
        epsMinus(1) = asin(sqrt((1-hMinus)/(1-kMinus(2)^2)));
        epsMinus(2) = asin(sqrt((1-hMinus)/(1-kMinus(2)^2)));
        
        % "Suitable for regular or singular cases" --doesn't work.
%         HrMinusComp = (-1)^(n-1)*sqrt(rout/r)*(2*t(n)*kMinus(n)/(r+rout)*KStarCalc(kMinus(n))+...
%             pi*sqrt(rout/r)*sign(r-rout)*sign(t(n))*(1-Lambda0Calc(epsMinus(n),kMinus(n))));
        
        
        % Regular Hr
        HrMinusComp = (-1)^(n-1)*t(n)*kMinus(n)/r*sqrt(rout/r)*(KStarCalc(kMinus(n))+...
            (r-rout)/(r+rout)*PiStarCalc(hMinus,pi/2,kMinus(n)));
        
    end
    function Lambda0 = Lambda0Calc(phi,m)
        Lambda0 = FStarCalc(phi,1-m)/(KStarCalc(1-m))+2/pi*KStarCalc(m)*ZetaCalc(phi,1-m);
    end
    function Zeta = ZetaCalc(phi,m)
        Zeta = EStarCalc(phi,m)-(EStarCalc(pi/2,m)*FStarCalc(phi,m))/KStarCalc(m);
    end
    function EStar = EStarCalc(phi,m)
        ellipticIntegral = @(theta) sqrt(1-m.^2.*sin(theta).^2);
        EStar = integral(ellipticIntegral,0,phi);
    end


end




















