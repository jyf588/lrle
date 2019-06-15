function [FM, lMtilde, FMactFL, FMactFV, FMpas, cos_alpha] = HillModel_RigidTendon_my(a,lMT,vMT,Parameters)

ActiveFVParameters = [-0.215812499592629	-32.5966241739010	-1.12412150864786	0.912644059915004];
PassiveFLParameters = [-0.995172050006169 53.5981500331442];
Faparam = [0.814483478343008 1.05503342897057 0.162384573599574 0.0633034484654646 0.433004984392647 0.716775413397760 -0.0299471169706956 0.200356847296188];

Fmax = Parameters.maxIsoForce;
lMopt = Parameters.optFL;
lTs = Parameters.tendonSL;
alphaopt = Parameters.alpha;
vMmax  = Parameters.maxConVel;

%%%% multiply vMax(10.0) with lMopt
vMmax = vMmax .* lMopt;

% Hill-type muscle model: geometric relationships
w = lMopt.*sin(alphaopt);
lM = sqrt((lMT-lTs).^2+w.^2); % Rigid Tendon: lT = lTs
lMtilde = lM./lMopt;
lT = lTs*ones(size(lMtilde));
cos_alpha = (lMT-lT)./lM;

b11 = Faparam(1);
b21 = Faparam(2);
b31 = Faparam(3);
b41 = Faparam(4);
b12 = Faparam(5);
b22 = Faparam(6);
b32 = Faparam(7);
b42 = Faparam(8);

b13 = 0.1;
b23 = 1;
b33 = 0.5*sqrt(0.5);
b43 = 0;
num3 = lMtilde-b23;
den3 = b33+b43*lMtilde;
FMtilde3 = b13*exp(-0.5*num3.^2./den3.^2);

num1 = lMtilde-b21;
den1 = b31+b41*lMtilde;
FMtilde1 = b11*exp(-0.5*num1.^2./den1.^2);

num2 = lMtilde-b22;
den2 = b32+b42*lMtilde;
FMtilde2 = b12*exp(-0.5*num2.^2./den2.^2);

FMactFL = FMtilde1+FMtilde2+FMtilde3;

% cos_alpha(cos_alpha<0) = 0.0;   %%%%
% Active muscle force-velocity relationship
vMtilde = (vMT./vMmax).*cos_alpha;
e1 = 1.475*ActiveFVParameters(1);
e2 = 0.25*ActiveFVParameters(2);
e3 = ActiveFVParameters(3) + 0.75;
e4 = ActiveFVParameters(4) - 0.027;

FMactFV = e1*log((e2*vMtilde+e3)+sqrt((e2*vMtilde+e3).^2+1))+e4;

% Active muscle force
FMact = a.*FMactFL.*FMactFV;

% Passive muscle force-length characteristic
e0 = 0.6;
kpe = 4;
t5 = exp(kpe * (lMtilde - 1.0) / e0);
FMpas = (t5  - PassiveFLParameters(1)) / PassiveFLParameters(2);

% Muscle force
FM_norm = (FMact+FMpas);
FM = Fmax*FM_norm;

return