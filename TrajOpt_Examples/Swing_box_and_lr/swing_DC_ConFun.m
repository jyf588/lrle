function c = swing_DC_ConFun(X,auxdata)
%
% Computes the nonlinear equality and inequality contraints
%	X       = current set of optimization parameters
%	auxdata = extra parameters and pointers to
%            instantiated OpenSim objects

% g = 9.80665;

% Import the OpenSim modeling classes
import org.opensim.modeling.*

% extract the nesessary auxiliary data
N         = auxdata.N;
h         = auxdata.h;
Nstates   = auxdata.Nstates;
Ncontrols = auxdata.Ncontrols;
Ncoord   = auxdata.Ncoord;
osimModel = auxdata.model;
dc_time   = auxdata.time;
init_q = auxdata.init_q;
tor_limit = auxdata.tor_limit;

addNN = auxdata.addNN;
Wall = auxdata.Wall;
Ball = auxdata.Ball;

osimState = osimModel.updWorkingState(); 

% Here, we extract all of the states and controls by node (time step) from X.
% This will result in two matrices, states and controls, that have a
% number of rows equal to the number of nodes and a number of columns
% equal to the number of states/controls, as appropriate.

states=zeros(N,Nstates); %pre-allocate size
for i = 1:Nstates
    states(:,i) = X(N*(i-1)+1:N*i,1); %column: state; row: nodes (time steps)
end

controls = zeros(N,Ncontrols); %pre-allocate size
for i = 1:Ncontrols
    controls(:,i) = X(Nstates*N + N*(i-1)+1:Nstates*N + N*i,1); %column: controls; row: nodes (time steps)
end


% ----------- Compute the equality constraints, ceq ---------------------------

ceq = [];

% Compute the constraint violation using backward Euler method
states_dot = zeros(N-1,Nstates);
for i = 1:Nstates
   for j = 1:N-1
      states_dot(j,i) = (states(j+1,i)-states(j,i))/h;
   end
end

% Get state derivatives from OpenSim
x_dot = zeros(N-1,Nstates);
for i = 1:N-1
   x_dot(i,:) = computeOpenSimModelXdot(states(i+1,:)',controls(i+1,:)',dc_time(i+1),osimModel, osimState)';
end

% Evaluate the constraint violoation at every time step except the last node, N
ceq_temp= zeros(N-1,Nstates);
for i = 1:Nstates 
    ceq_temp(:,i) = states_dot(1:end,i) - x_dot(1:end,i);
end

ceq_dyn = zeros((N-1)*Nstates, 1);  % Pre-allocate some space
% Re-arrange the constraint violations in one long vector
for i = 1:Nstates
    ceq_dyn((N-1)*(i-1)+1:(N-1)*i,:) = ceq_temp(:,i);
end
ceq = [ceq; ceq_dyn];

ceq_task = zeros(Nstates, 1);
% Task constraints: start from fixed state
for i = 1:Ncoord
    ceq_task(2*i-1,1) = states(1,i) - init_q(i);
    ceq_task(2*i,1) = states(1,i + Nstates/2);
end
ceq = [ceq; ceq_task];

% Nonlinear inequality contraints (cn):
cn = [];
% % add GRF constraints (Fy >= 0)
% cn_grf = zeros(N-1,1);
% for i = 1:N-1
%     q_ddq = [states(i,1:Nstates/2) states_dot(i,(Nstates/2+1):Nstates)];
%     com_acc_y = calcCOMaccy(q_ddq);
%     cn_grf(i,1) = -com_acc_y-g;    % com_ddy >= -9.85
% end
% % append to existing constr
% cn = [cn;cn_grf];
% 
if (addNN)
    % add my NN constraints: each time step -1<NN_{123}(x,u)<1
    cn_mynn = zeros(3*N,1);
    %[N123(t=1), (t=2), (t=3),...]
    for i = 1:N
        % hip, knee, ankle
%         y = NN_FD_tanh([states(i,[2:4 6:8]) controls(i,:)]',Wall,Ball,tor_limit);
        y = NN_FD([states(i,[3:5 8:10]) controls(i,2:4)]',Wall,Ball,tor_limit);
        cn_mynn(3*(i-1)+1 : 3*i) = y;
    end
    cn = [cn;cn_mynn];
end

c = [ceq; cn];
% 
% % inlined function
% function accy = calcCOMaccy(q_ddq)
%     for l = 0:length(q_ddq)-1
%        osimState.updY().set(l, q_ddq(l+1)); % use the second last vel and pos
%     end
% %     osimModel.computeStateVariableDerivatives(osimState); 
%     accy = osimModel.calcMassCenterVelocity(osimState).get(1);
% end

end

function y = NN_FD(x,Wall,Ball,tor_limit)
% pass Wall and Ball as ausdata.Wall(Ball), where Wall(Ball) is a cell array
N = length(Wall);
x(7:9) = x(7:9).*tor_limit;
% x(4:5) = x(4:5)/7.5;
% x(6) = x(6)/2.5;
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

% function y = NN_FD_tanh(x,Wall,Ball,tor_limit)
% % pass Wall and Ball as ausdata.Wall(Ball), where Wall(Ball) is a cell array
% N = length(Wall);
% x(7:9) = x(7:9).*tor_limit;
% x(4:5) = x(4:5)/7.5;
% x(6) = x(6)/2.5;
% x(7:9) = x(7:9)/100.0;
% for i = 1:N
%    W = Wall{i};
%    B = Ball{i};
%    wxb = double(W)*x + double(B);
%    if (i<N)
%        x = tanh(wxb);
%    else
%        x = wxb;
%    end
% end
% y = x;
% end