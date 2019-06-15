function c = jumperLeftRevoluteFootStatic_DC_ConFun(X,auxdata)
%
% Computes the nonlinear equality and inequality contraints
%	X       = current set of optimization parameters
%	auxdata = extra parameters and pointers to
%            instantiated OpenSim objects

g = 9.80665;
ctrl_q = 2:4;
ctrl_dq = 6:8;

% Import the OpenSim modeling classes
import org.opensim.modeling.*

% extract the nesessary auxiliary data
N         = auxdata.N;
h         = auxdata.h;
Nstates   = auxdata.Nstates;
Ncontrols = auxdata.Ncontrols;
osimModel = auxdata.model;
dc_time   = auxdata.time;
init_q = auxdata.init_q;
tor_limit = auxdata.tor_limit;


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
for i = 1:Nstates/2
    ceq_task(2*i-1,1) = states(1,i) - init_q(i);
    ceq_task(2*i,1) = states(1,i + Nstates/2);
end
ceq = [ceq; ceq_task];

% Nonlinear inequality contraints (cn):-------------------------------
cn = [];
% add GRF constraints (Fy >= 0)
cn_grf = zeros(N-1,1);
for i = 1:N-1
    q_ddq = [states(i,1:Nstates/2) states_dot(i,(Nstates/2+1):Nstates)];
    com_acc_y = calcCOMaccy(q_ddq);
    cn_grf(i,1) = -com_acc_y-g;    % com_ddy >= -9.85
end
% append to existing constr
cn = [cn;cn_grf];

c = [ceq; cn];

% add muscle equliaty constraints
musModel = auxdata.musModel;
DOF = auxdata.DOF;
MusclesOI = auxdata.MusclesOI;
CoordinatesOI = auxdata.CoordinatesOI;
musState = musModel.updWorkingState(); 
ceq_mus = zeros(Ncontrols*N,1);
M = length(MusclesOI);
act = zeros(N,M);
for i = 1:M
    act(:,i) = X((Ncontrols+Nstates)*N + (i-1)*N+1 : (Ncontrols+Nstates)*N + i*N);
end

for i = 1:N
    [tau,~] = ComputeMuscleTorques(musModel, musState, MusclesOI, CoordinatesOI, DOF, act(i,:)', states(i,ctrl_q)', states(i,ctrl_dq)');
    ceq_mus((i-1)*Ncontrols+1 : i*Ncontrols) = tau - controls(i,:)'.*(tor_limit');  % f(a,l,dl) - tau == 0
end
c = [c; ceq_mus];

% inlined function
function accy = calcCOMaccy(q_ddq)
    for l = 0:length(q_ddq)-1
       osimState.updY().set(l, q_ddq(l+1)); % use the second last vel and pos
    end
%     osimModel.computeStateVariableDerivatives(osimState); 
    accy = osimModel.calcMassCenterVelocity(osimState).get(1);
end

end