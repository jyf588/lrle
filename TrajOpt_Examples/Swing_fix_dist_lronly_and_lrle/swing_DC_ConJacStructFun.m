function jacStruct = swing_DC_ConJacStructFun(auxdata)

% extract the nesessary auxiliary data
N         = auxdata.N;
h         = auxdata.h;
Nstates   = auxdata.Nstates;
Ncontrols = auxdata.Ncontrols;
Ncoord    = auxdata.Ncoord;

addNN = auxdata.addNN;


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

gradceq_dyn = sparse(sj, si, ones(1, length(si)), (N-1)*Nstates, (Nstates+Ncontrols)*N);  % ipopt transpose from fmincon
gradceq = [gradceq; gradceq_dyn];

sj = zeros(1,Nstates);
si = zeros(1,Nstates);
sp = 1;
% Task constraints: start from fixed state
for i = 1:Ncoord
    sj(sp) = 2*i-1;
    si(sp) = (i-1)*N + 1;
    sp = sp+1;

    sj(sp) = 2*i;
    si(sp) = (i + Nstates/2 - 1)*N + 1;
    sp = sp+1;
end

gradceq_task = sparse(sj, si, ones(1, length(si)), Nstates, (Nstates+Ncontrols)*N); % ipopt transpose from fmincon
gradceq = [gradceq; gradceq_task];

% add final vx*vy constraint
sj = zeros(1, Nstates);
si = zeros(1, Nstates);
for i = 1:Nstates
    sj(i) = 1;
    si(i) = i*N-1;
end
gradf = sparse(sj, si, ones(1, length(si)), 1, (Nstates+Ncontrols)*N);
gradceq = [gradceq; gradf];

% Nonlinear inequality contraints (c):
gradcn = [];
% 
% % add gradient for GRF constraints (Fy >= 0) gradcn2   
% len = (Nstates/2*3)*(N-1); % dq_i dq_{i+1} q_i
% SI = zeros(1,len);
% SJ = zeros(1,len);
% SP = 1;
% for i = 1:N-1
%     SJ(SP:SP+(Nstates/2*3)-1) = i;
%     SI(SP:SP+Nstates/2-1) = (0 : Nstates/2-1) * N + i; % qi
%     SI(SP+Nstates/2 : SP+Nstates-1) = (Nstates/2 : Nstates-1) * N + i; % dqi
%     SI(SP+Nstates : SP+(Nstates/2*3)-1) = (Nstates/2 : Nstates-1) * N + i+1; % dq{i+1}
%     SP = SP+(Nstates/2*3);
% end
% gradcn2 = sparse(SJ, SI, ones(1,length(SI)), N-1, (Nstates+Ncontrols)*N); % ipopt transpose from fmincon
% gradcn = [gradcn; gradcn2];
% 
if (addNN)
    % add my NN constraints: each time step -1<NN_{123}(x,u)<1
    Nsc = 3*3;   % *3 means u,dq,q
    len = 3*N*Nsc;
    SI = zeros(1,len);
    SJ = zeros(1,len);
    SP = 1;
    %[N123-1(t=1), (t=2), (t=3),...]
    for i = 1:N
        % hip, knee, ankle
        for j = 1 : 3
            SI(SP:SP+Nsc-1) = ([3:5 8:10 12:14]-1)*N + i;
            SJ(SP:SP+Nsc-1) = 3*(i-1) + j;
            SP = SP+Nsc;
        end
    end
    gradcn_mynn = sparse(SJ,SI, ones(1,length(SI)), 3*N, (Nstates+Ncontrols)*N);
    gradcn = [gradcn; gradcn_mynn];
end

jacStruct = [gradceq; gradcn];

end