function f = swing_DC_ObjFun(X,auxdata)
%
% Computes the objective function value for the  simpleLowerLimb
%	X = current set of optimization parameters
%   auxdata = extra parameters and pointers to
%             instantiated OpenSim objects
%   f = sum of squared muscle activation integrals
% extract the nesessary auxiliary data

% Import the OpenSim modeling classes
import org.opensim.modeling.*

% del = sqrt(eps);    % findiff step size

N        = auxdata.N;
% osimModel = auxdata.model;
% h        = auxdata.h;
Nstates  = auxdata.Nstates;
Ncontrols = auxdata.Ncontrols;
% Ncoord   = auxdata.Ncoord;
tor_limit = auxdata.tor_limit;

useMeta   = auxdata.useMeta;
Wall_E = auxdata.Wall_E;
Ball_E = auxdata.Ball_E;

states=zeros(N,Nstates); %pre-allocate size
for i = 1:Nstates
    states(:,i) = X(N*(i-1)+1:N*i,1); %column: state; row: nodes (time steps)
end

controls = zeros(N,Ncontrols); %pre-allocate size
for i = 1:Ncontrols
    controls(:,i) = X(Nstates*N + N*(i-1)+1:Nstates*N + N*i,1); %column: controls; row: nodes (time steps)
end

if (useMeta)
    f = 0;
    for i = 1:N
        % hip, knee, ankle
        y = NN_FD([states(i,[3:5 8:10]) controls(i,2:4)]',Wall_E,Ball_E,tor_limit);
        f = f + y / 100 * 3;      % max wanted energy value 100 (1000)
    end
    f = f + sum(abs(controls(:,1))) / 2;
%     f = f + sum(abs(controls(:,1)));
else
    controls_scaled = controls;
    controls_scaled(:,1) = controls(:,1) / 2;   % shoulder torques are smaller
    f = sum(abs(controls_scaled(:)));
%     f = sum(abs(controls(:)));
end

% f = f * 100;s
global hisy
hisy = [hisy; reshape(X,1,[])];

end

function y = NN_FD(x,Wall,Ball,tor_limit)
% pass Wall and Ball as ausdata.Wall(Ball), where Wall(Ball) is a cell array
N = length(Wall);
x(7:9) = x(7:9).*tor_limit;
x(5:6) = x(5:6)/7.5;
x(4) = x(4)/2.5;
x(7:9) = x(7:9)/100.0;
for i = 1:N
   W = Wall{i};
   B = Ball{i};
   wxb = double(W)*x + double(B);
   if (i<N)
       x = exp(wxb.*(wxb<=0))-1 + wxb.*(wxb>0);
   else
       x = wxb;
   end
end
y = x;
end

