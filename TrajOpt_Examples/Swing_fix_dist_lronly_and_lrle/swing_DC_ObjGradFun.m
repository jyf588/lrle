function gradf = swing_DC_ObjGradFun(X,auxdata)

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

gradf = zeros(length(X),1);
if (useMeta)
    for i = 1:N
        % dydx is 1 by 9
        dydx = NN_BP([states(i,[3:5 8:10]) controls(i,2:4)]',Wall_E,Ball_E,tor_limit);
        gradf(([3:5 8:10 12:14]-1)*N + i,:) = dydx' / 100 * 3;
        gradf((11-1)*N + i,:) = sign(X((11-1)*N + i)) / 2.0;
%         gradf((11-1)*N + i,:) = sign(X((11-1)*N + i));
    end
else
    gradf((Nstates*N+1):end) = sign(X((Nstates*N+1):end));
    gradf((11-1)*N + (1:N),:) = gradf((11-1)*N + (1:N),:) / 2.0;
end

% gradf = gradf * 100;
gradf = sparse(gradf);
end

function dydx = NN_BP(x,Wall,Ball,tor_limit)
% pass Wall and Ball as ausdata.Wall(Ball), where Wall(Ball) is a cell array
N = length(Wall);
dydx = eye(length(x));
x(7:9) = x(7:9).*tor_limit;
x(5:6) = x(5:6)/7.5;
x(4) = x(4)/2.5;
x(7:9) = x(7:9)/100.0;
for i = 1:N
   W = Wall{i};
   B = Ball{i};
   wxb = double(W)*x + double(B);
   if (i<N)
       dydx = repmat(exp(min(0,wxb)),1,size(W,2)) .* double(W) * dydx;
       x = exp(wxb.*(wxb<=0))-1 + wxb.*(wxb>0);
   else
       dydx = double(W) * dydx;
       x = wxb;
   end
end
dydx(:,5:6) = dydx(:,5:6)/7.5;
dydx(:,4) = dydx(:,4)/2.5;
dydx(:,7:9) = dydx(:,7:9)/100.0;
dydx(:,7:9) = dydx(:,7:9) .* repmat(tor_limit', size(dydx,1), 1);
end