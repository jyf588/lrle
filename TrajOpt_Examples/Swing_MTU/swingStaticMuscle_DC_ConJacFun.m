function jac = swingStaticMuscle_DC_ConJacFun(X, auxdata)

% g = 9.80665;
del = sqrt(eps);    % step size
% the leg q dq's
ctrl_q = 3:5;
ctrl_dq = 8:10;
ctrl_states = [ctrl_q ctrl_dq];
% the leg tau
ctrl_tau = 2:4;


% Import the OpenSim modeling classes
import org.opensim.modeling.*

% extract the nesessary auxiliary data
N         = auxdata.N;
h         = auxdata.h;
Nstates   = auxdata.Nstates;
Ncontrols = auxdata.Ncontrols;
Ncoord    = auxdata.Ncoord;
osimModel = auxdata.model;
dc_time   = auxdata.time;
% init_q = auxdata.init_q;
tor_limit = auxdata.tor_limit;

parModels = auxdata.c;
parModels_mus = auxdata.c2;

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

gradceq = [];

% ----------- Compute the equality constraints, ceq ---------------------------

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

% first fill in the entries dC(i,j)/d q(k,j+1) k = 1..8
si_mat = zeros(N-1, Nstates*Nstates);
sj_mat = zeros(N-1, Nstates*Nstates);
sv_mat = zeros(N-1, Nstates*Nstates);
parfor j = 1:(N-1)
    parState_i = parModels.Value.updWorkingState();
    si_t = zeros(1,Nstates*Nstates);
    sj_t = zeros(1,Nstates*Nstates);
    sv_t = zeros(1,Nstates*Nstates);
    sp_t = 1;
    for k = 1:Nstates
        % for all i = 1..8, they are all related to dyneqn f
        pstate = states(j+1,:)';
        pstate(k) = pstate(k) + del; % forward diff
        pxdot = computeOpenSimModelXdot(pstate, controls(j+1,:)',dc_time(j+1),parModels.Value,parState_i)';
        pceq = (pxdot - x_dot(j,:)) / del * (-1); % *-1 since ceq = state_dot - xdot
        % there are 6 pceq elements here, insert them.
        for i = 1:Nstates
            sj_t(sp_t) = (i-1)*(N-1)+j; % sj is the # of ceq
            si_t(sp_t) = (k-1)*N + (j+1); % si is the # of state x=(q,u)
            sv_t(sp_t) = pceq(i);
            if (i == k)
                % special case: dC(i,j)/d q(i,j+1) has additional 1/h
                sv_t(sp_t) = sv_t(sp_t) + 1.0/h;
            end
            sp_t = sp_t +1;
        end
    end
    si_mat(j,:) = si_t;
    sj_mat(j,:) = sj_t;
    sv_mat(j,:) = sv_t;
end
si = si_mat(:)';
sj = sj_mat(:)';
sv = sv_mat(:)';

% similarly, fill in dC(i,j)/d u(k,j+1) k = 1..3
si_mat = zeros(N-1, Nstates*Ncontrols);
sj_mat = zeros(N-1, Nstates*Ncontrols);
sv_mat = zeros(N-1, Nstates*Ncontrols);
parfor j = 1:(N-1)
    parState_i = parModels.Value.updWorkingState();
    si_t = zeros(1,Nstates*Ncontrols);
    sj_t = zeros(1,Nstates*Ncontrols);
    sv_t = zeros(1,Nstates*Ncontrols);
    sp_t = 1;
    for k = 1:Ncontrols
        pctrl = controls(j+1,:)';
        pctrl(k) = pctrl(k) + del;
        pxdot = computeOpenSimModelXdot(states(j+1,:)', pctrl, dc_time(j+1),parModels.Value,parState_i)';
        pceq = (pxdot - x_dot(j,:)) / del * (-1);
        % there are 8 pceq elements here, insert them.
        for i = 1:Nstates
            sj_t(sp_t) = (i-1)*(N-1)+j; % sj is the # of ceq
            si_t(sp_t) = (Nstates + k-1)*N + (j+1); % si is the # of state x=(q,u)
            sv_t(sp_t) = pceq(i);
            sp_t = sp_t +1;
        end
    end
    si_mat(j,:) = si_t
    sj_mat(j,:) = sj_t;
    sv_mat(j,:) = sv_t;
end
si = [si si_mat(:)'];
sj = [sj sj_mat(:)'];
sv = [sv sv_mat(:)'];
sp = (N-1)*(Nstates+Ncontrols)*Nstates + 1;

% we have filled in 24*(8+3)*8, we are left with 24*8*1
for i = 1:Nstates
    for j = 1:(N-1)
        sj(sp) = (i-1)*(N-1)+j;
        si(sp) = (i-1)*N+j;
        sv(sp) = -1.0/h;
        sp = sp+1;
    end
end

gradceq_dyn = sparse(sj, si, sv, (N-1)*Nstates, length(X));  % ipopt transpose from fmincon
gradceq = [gradceq; gradceq_dyn];


sj = zeros(1,Nstates);
si = zeros(1,Nstates);
sv = zeros(1,Nstates);
sp = 1;
% Task constraints: start from fixed state
for i = 1:Ncoord
    sj(sp) = 2*i-1;
    si(sp) = (i-1)*N + 1;
    sv(sp) = 1.0;
    sp = sp+1;

    sj(sp) = 2*i;
    si(sp) = (i + Nstates/2 - 1)*N + 1;
    sv(sp) = 1.0;
    sp = sp+1;
end

gradceq_task = sparse(sj, si, sv, Nstates, length(X)); % ipopt transpose from fmincon
gradceq = [gradceq; gradceq_task];

% % Nonlinear inequality contraints (cn):-----------------------------
% gradcn = [];
% 
% % add GRF constraints (Fy >= 0)
% cn_grf = zeros(N-1,1);
% for i = 1:N-1
%     q_ddq = [states(i,1:Nstates/2) states_dot(i,(Nstates/2+1):Nstates)];
%     com_acc_y = calcCOMaccy(q_ddq);
%     cn_grf(i,1) = -com_acc_y-g;
% end
% 
% % add gradient for GRF constraints (Fy >= 0) gradcn2   
% len = (Nstates/2*3)*(N-1); % dq_i dq_{i+1} q_i
% SI = zeros(1,len);
% SJ = zeros(1,len);
% SV = zeros(1,len);
% SP = 1;
% 
% for i = 1:N-1
%     SJ(SP:SP+(Nstates/2*3)-1) = i;
% 
%     SI(SP:SP+Nstates/2-1) = (0 : Nstates/2-1) * N + i; % qi
%     q_ddq = [states(i,1:Nstates/2) states_dot(i,(Nstates/2+1):Nstates)];
%     for j = 1:Nstates/2
%         pq_ddq = q_ddq;
%         pq_ddq(j) = pq_ddq(j) + del;
%         pcom_acc_y = calcCOMaccy(pq_ddq);
%         SV(SP+j-1) = ((-pcom_acc_y-g) - cn_grf(i,1)) / del;
%     end
% 
%     SI(SP+Nstates/2 : SP+Nstates-1) = (Nstates/2 : Nstates-1) * N + i; % dqi
%     SI(SP+Nstates : SP+(Nstates/2*3)-1) = (Nstates/2 : Nstates-1) * N + i+1; % dq{i+1}
%     q_0 = [states(i,1:Nstates/2) zeros(1,Nstates/2)];
%     for j = 1:Nstates/2
%         pq_0 = q_0;
%         pq_0(j+Nstates/2) = 1/h; % note additional sign change from -pcom_acc_y-g
%         dcdqi = calcCOMaccy(pq_0);
%         SV(SP+Nstates/2+j-1) = dcdqi;   % dqi
%         SV(SP+Nstates+j-1) = -dcdqi; % dq{i+1}
%     end
% 
%     SP = SP+(Nstates/2*3);
% end
% gradcn_grf = sparse(SJ, SI, SV, N-1, length(X)); % ipopt transpose from fmincon
% gradcn = [gradcn; gradcn_grf];
% 
% jac = [gradceq; gradcn];

% add muscle equality constraints
musModel = auxdata.musModel;
DOF = auxdata.DOF;
nDOF = length(DOF);
MusclesOI = auxdata.MusclesOI;
CoordinatesOI = auxdata.CoordinatesOI;
musState = musModel.updWorkingState(); 
ceq_mus_mat = zeros(nDOF,N);
M = length(MusclesOI);
act = zeros(N,M);
for i = 1:M
    act(:,i) = X((Ncontrols+Nstates)*N + (i-1)*N+1 : (Ncontrols+Nstates)*N + i*N);
end
len = N*nDOF*M;
si = zeros(1,len);
sj = zeros(1,len);
sv = zeros(1,len);
sp = 1;
for i = 1:N
    % actually, not all act related to every tau; simplify sparsity calc here
    [tau,dTaudA] = ComputeMuscleTorques(musModel, musState, MusclesOI, CoordinatesOI, DOF, act(i,:)', states(i,ctrl_q)', states(i,ctrl_dq)');
    ceq_mus_mat(:,i) = tau - controls(i,ctrl_tau)'.*(tor_limit');
    for j = 1 : nDOF
        si(sp:sp+M-1) = ((Ncontrols+Nstates)+(0:M-1))*N + i;
        sj(sp:sp+M-1) = nDOF*(i-1) + j;
        sv(sp:sp+M-1) = dTaudA(j,:)';
        sp = sp+M;
    end
end

si_tau = zeros(1,nDOF*N);
sj_tau = zeros(1,nDOF*N);
sv_tau = zeros(1,nDOF*N);
for i = 1:N
    for j = 1:nDOF
        si_tau((i-1)*nDOF+j) = (Nstates+ctrl_tau(j)-1)*N + i;
        sj_tau((i-1)*nDOF+j) = nDOF*(i-1) + j;
        sv_tau((i-1)*nDOF+j) = -tor_limit(j);
    end
end
si = [si si_tau];
sj = [sj sj_tau];
sv = [sv sv_tau];

% first fill in the entries dC(j,controls)/d (q/dq)(j,states)
si_mat = zeros(N, nDOF*(nDOF*2)); % Ncontrols*2: all the states(qdq) that are controlled
sj_mat = zeros(N, nDOF*(nDOF*2));
sv_mat = zeros(N, nDOF*(nDOF*2));
parfor j = 1:N
    parState_mus_i = parModels_mus.Value.updWorkingState();
    si_t = zeros(1,nDOF*(nDOF*2));
    sj_t = zeros(1,nDOF*(nDOF*2));
    sv_t = zeros(1,nDOF*(nDOF*2));
    sp_t = 1;
    for k = 1:nDOF*2
        pstate = states(j,ctrl_states)';
        pstate(k) = pstate(k) + del; % forward diff
        [ptau,~] = ComputeMuscleTorques(parModels_mus.Value, parState_mus_i, MusclesOI, CoordinatesOI, DOF, act(j,:)', pstate(1:nDOF), pstate(nDOF+1:end));
        pceq_mus = (ptau - controls(j,ctrl_tau)'.*(tor_limit') - ceq_mus_mat(:,j)) / del; 
        % there are Ncontorls pceq elements here, insert them.
        for i = 1:nDOF
            sj_t(sp_t) = (j-1)*nDOF + i; % sj is the # of ceq
            si_t(sp_t) = (ctrl_states(k)-1)*N + j; % si is the # of state x=(q,u)
            sv_t(sp_t) = pceq_mus(i);
            sp_t = sp_t +1;
        end
    end
    si_mat(j,:) = si_t;
    sj_mat(j,:) = sj_t;
    sv_mat(j,:) = sv_t;
end
si = [si si_mat(:)'];
sj = [sj sj_mat(:)'];
sv = [sv sv_mat(:)'];
gradceq_mus = sparse(sj, si, sv, nDOF*N, length(X)); % ipopt transpose from fmincon

jac = [gradceq; gradceq_mus];

% % inlined function
% function accy = calcCOMaccy(q_ddq)
%     for l = 0:length(q_ddq)-1
%        osimState.updY().set(l, q_ddq(l+1)); % use the second last vel and pos
%     end
% %     osimModel.computeStateVariableDerivatives(osimState); 
%     accy = osimModel.calcMassCenterVelocity(osimState).get(1);
% end

end