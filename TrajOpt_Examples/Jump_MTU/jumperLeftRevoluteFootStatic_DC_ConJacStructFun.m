function jacStruct = jumperLeftRevoluteFootStatic_DC_ConJacStructFun(auxdata)

% extract the nesessary auxiliary data
N         = auxdata.N;
h         = auxdata.h;
Nstates   = auxdata.Nstates;
Ncontrols = auxdata.Ncontrols;

MusclesOI = auxdata.MusclesOI;
M = length(MusclesOI);
X = ones((Ncontrols+Nstates+M)*N,1); % dummy

ctrl_q = 2:4;
ctrl_dq = 6:8;
ctrl_states = [ctrl_q ctrl_dq];

gradceq = [];

% first fill in the entries dC(i,j)/d q(k,j+1) k = 1..8
si_mat = zeros(N-1, Nstates*Nstates);
sj_mat = zeros(N-1, Nstates*Nstates);
for j = 1:(N-1)
    si_t = zeros(1,Nstates*Nstates);
    sj_t = zeros(1,Nstates*Nstates);
    sp_t = 1;
    for k = 1:Nstates
        for i = 1:Nstates
            sj_t(sp_t) = (i-1)*(N-1)+j; % sj is the # of ceq
            si_t(sp_t) = (k-1)*N + (j+1); % si is the # of state x=(q,u)
            sp_t = sp_t +1;
        end
    end
    si_mat(j,:) = si_t;
    sj_mat(j,:) = sj_t;
end
si = si_mat(:)';
sj = sj_mat(:)';

% similarly, fill in dC(i,j)/d u(k,j+1) k = 1..3
si_mat = zeros(N-1, Nstates*Ncontrols);
sj_mat = zeros(N-1, Nstates*Ncontrols);
for j = 1:(N-1)
    si_t = zeros(1,Nstates*Ncontrols);
    sj_t = zeros(1,Nstates*Ncontrols);
    sp_t = 1;
    for k = 1:Ncontrols
        for i = 1:Nstates
            sj_t(sp_t) = (i-1)*(N-1)+j; % sj is the # of ceq
            si_t(sp_t) = (Nstates + k-1)*N + (j+1); % si is the # of state x=(q,u)
            sp_t = sp_t +1;
        end
    end
    si_mat(j,:) = si_t;
    sj_mat(j,:) = sj_t;
end
si = [si si_mat(:)'];
sj = [sj sj_mat(:)'];
sp = (N-1)*(Nstates+Ncontrols)*Nstates + 1;

% we have filled in 24*(8+3)*8, we are left with 24*8*1
for i = 1:Nstates
    for j = 1:(N-1)
        sj(sp) = (i-1)*(N-1)+j;
        si(sp) = (i-1)*N+j;
        sp = sp+1;
    end
end

gradceq_dyn = sparse(sj, si, ones(1, length(si)), (N-1)*Nstates, length(X));  % ipopt transpose from fmincon
gradceq = [gradceq; gradceq_dyn];

sj = zeros(1,Nstates);
si = zeros(1,Nstates);
sp = 1;
% Task constraints: start from fixed state
for i = 1:Nstates/2
    sj(sp) = 2*i-1;
    si(sp) = (i-1)*N + 1;
    sp = sp+1;

    sj(sp) = 2*i;
    si(sp) = (i + Nstates/2 - 1)*N + 1;
    sp = sp+1;
end

gradceq_task = sparse(sj, si, ones(1, length(si)), Nstates, length(X)); % ipopt transpose from fmincon
gradceq = [gradceq; gradceq_task];

% Nonlinear inequality contraints (c):---------------------------------
gradcn = [];

% add gradient for GRF constraints (Fy >= 0) gradcn2   
len = (Nstates/2*3)*(N-1); % dq_i dq_{i+1} q_i
SI = zeros(1,len);
SJ = zeros(1,len);
SP = 1;
for i = 1:N-1
    SJ(SP:SP+(Nstates/2*3)-1) = i;
    SI(SP:SP+Nstates/2-1) = (0 : Nstates/2-1) * N + i; % qi
    SI(SP+Nstates/2 : SP+Nstates-1) = (Nstates/2 : Nstates-1) * N + i; % dqi
    SI(SP+Nstates : SP+(Nstates/2*3)-1) = (Nstates/2 : Nstates-1) * N + i+1; % dq{i+1}
    SP = SP+(Nstates/2*3);
end
gradcn2 = sparse(SJ, SI, ones(1,length(SI)), N-1, length(X)); % ipopt transpose from fmincon
gradcn = [gradcn; gradcn2];

jacStruct = [gradceq; gradcn];

% add muscle equality constraints
len = N*Ncontrols*M;
si = zeros(1,len);
sj = zeros(1,len);
sp = 1;
for i = 1:N
    % actually, not all a related to every tau; simplify sparsity calc here
    for j = 1 : Ncontrols
        si(sp:sp+M-1) = ((Ncontrols+Nstates)+(0:M-1))*N + i;
        sj(sp:sp+M-1) = Ncontrols*(i-1) + j;
        sp = sp+M;
    end
end

si_tau = zeros(1,Ncontrols*N);
sj_tau = zeros(1,Ncontrols*N);
for i = 1:N
    for j = 1:Ncontrols
        si_tau((i-1)*Ncontrols+j) = (Nstates+j-1)*N + i;
        sj_tau((i-1)*Ncontrols+j) = Ncontrols*(i-1) + j;
    end
end
si = [si si_tau];
sj = [sj sj_tau];

% first fill in the entries dC(j,controls)/d (q/dq)(j,states)
si_mat = zeros(N, Ncontrols*(Ncontrols*2)); % Ncontrols*2: all the states(qdq) that are controlled
sj_mat = zeros(N, Ncontrols*(Ncontrols*2));
for j = 1:N
    si_t = zeros(1,Ncontrols*(Ncontrols*2));
    sj_t = zeros(1,Ncontrols*(Ncontrols*2));
    sp_t = 1;
    for k = 1:Ncontrols*2
        for i = 1:Ncontrols
            sj_t(sp_t) = (j-1)*Ncontrols + i; % sj is the # of ceq
            si_t(sp_t) = (ctrl_states(k)-1)*N + j; % si is the # of state x=(q,u)
            sp_t = sp_t +1;
        end
    end
    si_mat(j,:) = si_t;
    sj_mat(j,:) = sj_t;
end
si = [si si_mat(:)'];
sj = [sj sj_mat(:)'];

gradceq_mus = sparse(sj, si, ones(1,length(si)), Ncontrols*N, length(X)); % ipopt transpose from fmincon

jacStruct = [jacStruct; gradceq_mus];

end