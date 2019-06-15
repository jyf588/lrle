function jac = swing_DC_ConJacFun(X, auxdata)

% g = 9.80665;
del = sqrt(eps);    % step size

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

addNN = auxdata.addNN;
Wall = auxdata.Wall;
Ball = auxdata.Ball;


parModels = auxdata.c;

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

gradceq_dyn = sparse(sj, si, sv, (N-1)*Nstates, (Nstates+Ncontrols)*N);  % ipopt transpose from fmincon
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

% add final vx*vy constraint
for i = 0:Nstates-1
   osimState.updY().set(i, X((i+1)*N-1,1)); % use the second last vel and pos
end
fvx = osimModel.calcMassCenterVelocity(osimState).get(0);
fvy = osimModel.calcMassCenterVelocity(osimState).get(1);
f = calcObj(fvx, fvy) + 4.0;

sv = zeros(1, Nstates);
sj = zeros(1, Nstates);
si = zeros(1, Nstates);
for i = 1:Nstates
    sj(i) = 1;
    si(i) = i*N-1;
    osimState.updY().set(i-1, X(i*N-1,1)+del);
    pfvx = osimModel.calcMassCenterVelocity(osimState).get(0);
    pfvy = osimModel.calcMassCenterVelocity(osimState).get(1);
    pf = calcObj(pfvx, pfvy) + 4.0;
    sv(i) = (pf - f)/del;
    osimState.updY().set(i-1, X(i*N-1,1));
end
gradf = sparse(sj, si, sv, 1, length(X));
gradceq = [gradceq; gradf];

% Nonlinear inequality contraints (c):------------------------------
gradcn = [];
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
% gradcn_grf = sparse(SJ, SI, SV, N-1, (Nstates+Ncontrols)*N); % ipopt transpose from fmincon
% gradcn = [gradcn; gradcn_grf];
% 
if (addNN)
    % add my NN constraints: each time step -1<NN_{123}(x,u)<1
    Nsc = 3*3;   % *3 means u,dq,q
    len = 3*N*Nsc;
    SI = zeros(1,len);
    SJ = zeros(1,len);
    SV = zeros(1,len);
    SP = 1;

    %[N123-1(t=1), (t=2), (t=3),...]
    for i = 1:N
        % hip, knee, ankle
%         dydx = NN_BP_tanh([states(i,[2:4 6:8]) controls(i,:)]',Wall,Ball,tor_limit);
        dydx = NN_BP([states(i,[3:5 8:10]) controls(i,2:4)]',Wall,Ball,tor_limit);
        for j = 1 : 3
            SI(SP:SP+Nsc-1) = ([3:5 8:10 12:14]-1)*N + i;
            SJ(SP:SP+Nsc-1) = 3*(i-1) + j;
            SV(SP:SP+Nsc-1) = dydx(j,:)';
            SP = SP+Nsc;
        end
    end
    gradcn_mynn = sparse(SJ,SI,SV, 3*N, (Nstates+Ncontrols)*N);
    gradcn = [gradcn; gradcn_mynn];
end

jac = [gradceq; gradcn];

% % inlined function
% function accy = calcCOMaccy(q_ddq)
%     for l = 0:length(q_ddq)-1
%        osimState.updY().set(l, q_ddq(l+1)); % use the second last vel and pos
%     end
%     accy = osimModel.calcMassCenterVelocity(osimState).get(1);
% end

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

% function dydx = NN_BP_tanh(x,Wall,Ball,tor_limit)
% % pass Wall and Ball as ausdata.Wall(Ball), where Wall(Ball) is a cell array
% N = length(Wall);
% dydx = eye(length(x));
% x(7:9) = x(7:9).*tor_limit;
% x(5:6) = x(5:6)/7.5;
% x(4) = x(4)/2.5;
% x(7:9) = x(7:9)/100.0;
% for i = 1:N
%    W = Wall{i};
%    B = Ball{i};
%    wxb = double(W)*x + double(B);
%    if (i<N)
%        dydx = repmat(1-tanh(wxb).^2,1,size(W,2)) .* double(W) * dydx;
%        x = tanh(wxb);
%    else
%        dydx = double(W) * dydx;
%        x = wxb;
%    end
% end
% dydx(:,5:6) = dydx(:,5:6)/7.5;
% dydx(:,4) = dydx(:,4)/2.5;
% dydx(:,7:9) = dydx(:,7:9)/100.0;
% dydx(:,7:9) = dydx(:,7:9) .* repmat(tor_limit', size(dydx,1), 1);
% end

function obj = calcObj(vx,vy)
    if (vx < 0 && vy > 0)
        % this is what we want
        obj = vx * vy; % since we are minimizing; up to a scale
    end
    % shortest l2 distance to desired convex region
    if (vx < 0 && vy <= 0)
        obj = -vy;
    end
    if (vx >= 0 && vy > 0)
        obj = +vx;
    end
    if (vx >= 0 && vy <= 0)
        obj = sqrt(vx^2+vy^2);
    end
end