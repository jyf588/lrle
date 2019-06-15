function [Torque_Manageable,act,eT,MusclesOI] = TorqueLimitClassifier(Model_OS, DOF, q, qdot, tau)
% This function evaluates whether the muscles can produce lower limb joint
% torques for a given kinematic state of the model.
%
% Inputs:
% modelName Name of the osim-model. The osim-model is supposed to be in the
%           current working directory. If this is not the case, edit
%           MAIN_path.
% DOF   cell array of strings of length nDof containing the names of the
%       degrees of freedom that are considered.
% q     N x nDof matrix of generalized coordinates (radians, meters)
%       N is the number of samples, nDof is the number of degrees of
%       freedom
% qdot  N x nDof matrix of generalized velocities
% tau   N x nDof matrix of torques (Newton meter)
%
% Outputs:
% Torque_Manageable     vector of length N containing 1 for feasible samples 
%                       and 0 for infeasible samples
% act                   N x M matrix with muscle activations (bounded
%                       between 0 and 1)
%                       M is the number of muscles actuating the degrees of
%                       freedom under consideration.
% eT                    N x nDof matrix with reserve torques. These are
%                       near zero when tau is feasible.
% MusclesOI             Structure holding information on the used muscles

import org.opensim.modeling.*
MAIN_path = pwd;
load ActiveFVParameters
load PassiveFLParameters
load Faparam

N = size(q,1);

nDof = length(DOF);
% Calculate all kinds of information based on the model and the DOFs you
% are interested in (OI = Of Interest)
[MusclesOI,CoordinatesOI] = getMuscleStructure(Model_OS,DOF);
M = length(MusclesOI);

% Initialize some matrices
act = ones(N,M); % Activations
FMltilde = ones(N,M); % Normalized Force-Length dependence
FMvtilde = ones(N,M); % Normalized Force-Velocity dependence
Fpe = ones(N,M);      % Passive (elastic) force
cos_alpha = ones(N,M); % cosine of the muscles pennation angle

kinematics = [q, qdot];

% Now we calculate based on the kinematic state the Muscle-Tendon Length & Velocity and the Moment-Arm Matrix
[LMT,VMT,dM] = get_LMT_vMT_dM(Model_OS,MusclesOI,CoordinatesOI,kinematics,DOF);

% We apply the Hill Model under the assumption of a rigid tendon and
% calculate the normalized quantities of passive and active muscle force
% components. The active component will be scaled with the muscle
% activation and added to the passive force to produce the final muscle force.
for m = 1:M
    [~, ~, FMltilde(:,m), FMvtilde(:,m), Fpe(:,m), cos_alpha(:,m)] = HillModel_RigidTendon(act(:,m),LMT(:,m),VMT(:,m),MusclesOI(m),ActiveFVParameters,PassiveFLParameters,Faparam);
end

FMo = ones(size(act,1),1)*[MusclesOI(:).maxIsoForce];
Fpas = FMo.*Fpe.*cos_alpha;
Fact = FMo.*FMltilde.*FMvtilde.*cos_alpha;


% Add optimal force for reserve torques
Topt = 1;
Fpas = [Fpas zeros(N,nDof)];  % The minimal force of our actuators (muscle force or torque)
Fact = [Fact Topt*ones(N,nDof)];

% Initial guess on the solution (activations of muscles and reserves
% minimal), this is our optimization vector
x0 = repmat(0.01*ones(M+nDof,1),N,1);


% The bounds
options.lb = repmat([zeros(M,1); -1500*ones(nDof,1)],N,1);  % Lower bounds of the optimization variables
options.ub = repmat([ones(M,1); 1500*ones(nDof,1)],N,1);    % Upper bounds of the optimization variables
options.cl = repmat(zeros(nDof,1),N,1); % Lower bounds on the constraints (comply with the torques that are provided)
options.cu = repmat(zeros(nDof,1),N,1); % Upper bounds on the constraints (comply with the torques that are provided)

I = N*(M+nDof);
% Set up the auxiliary data.
options.auxdata = { M N nDof reshape(Fact', I, 1)  reshape(Fpas', I, 1) ...
     tau dM};


options.ipopt.mu_strategy      = 'adaptive';
options.ipopt.max_iter         = 1500;
options.ipopt.tol              = 1e-5;
options.ipopt.hessian_approximation = 'limited-memory';
% options.ipopt.derivative_test = 'first-order';
options.ipopt.print_level = 2;

options.ipopt.acceptable_compl_inf_tol = 2.0;
options.ipopt.acceptable_constr_viol_tol = 1e-2;
options.ipopt.acceptable_obj_change_tol = 1e-3;
options.ipopt.acceptable_tol = 1e-1;
options.ipopt.acceptable_iter = 10;

% The callback functions.
funcs.objective         = @objective_SO;
funcs.constraints       = @constraints_SO;
funcs.gradient          = @gradient_SO;
funcs.jacobian          = @jacobian_SO;
funcs.jacobianstructure = @jacobianstructure_SO;

[x,~] = ipopt_auxdata(x0,funcs,options);

x_opt = reshape(x, M+nDof, N)';

act = x_opt(:,1:M);
eT = x_opt(:, M+1:M+nDof)*Topt;

ReserveActuatorSum = sum(abs(eT),2);
for i = 1:size(ReserveActuatorSum,1)
    if ReserveActuatorSum(i) < 1
        ReserveActuatorSum(i) = 1;
    else 
        ReserveActuatorSum(i) = 0;
    end
end


Torque_Manageable = ReserveActuatorSum;
end
% ------------------------------------------------------------------
function f = objective_SO (x, auxdata)
%f     = 0.5 * sum(x.^2);
[M, N, nDof, Fmax, Fpas, ID_data, MomentArm] = deal(auxdata{:});
x_opt = reshape(x, M+nDof, N)';
f     = 0.5 * (sum(sum(x_opt(:,1:M).^2)) + 10*sum(sum(x_opt(:,M+1:end).^2)));  % We want to minimize the optimization vector (activations of muscles and reserves)
% NOTE: We especially want to minimize the reserves to see whether the
% muscles can deliver the desired moment since the activations of the
% reserves are the torques itself and the activations of the muscles scale
% with an optimal force we don't need to put an extremely high weight on the term
% minimizing the reserves.
end

% ------------------------------------------------------------------
function c = constraints_SO (x, auxdata)
[M, N, nDof, Fmax, Fpas, ID_data, MomentArm] = deal(auxdata{:});

F = Fmax .* x + Fpas;   % Here we calculate the force of a muscle: active force component scaled with the activation + the passive component.

c = zeros(nDof*N,1);

for k = 1:nDof
    F_matrix = reshape(F, M+nDof, N)';
    MomentArm_matrix = reshape(MomentArm(:,k), M+nDof, N)';
    c(k:nDof:end) = sum(F_matrix.*MomentArm_matrix, 2) - ID_data(:,k);   % Inverse Dynamics Constraint - deliver the desired torque/moment
end
end
% ------------------------------------------------------------------
function g = gradient_SO (x, auxdata)
%g = x;
[M, N, nDof, Fmax, Fpas, ID_data, MomentArm] = deal(auxdata{:});
x_opt = reshape(x, M+nDof, N)';
gtemp = [x_opt(:,1:M) 10*x_opt(:,M+1:end)];   % Gradient of the objective to the variables
I = N*(M+nDof);
g = reshape(gtemp',I,1);
end
% ------------------------------------------------------------------
function J = jacobianstructure_SO (auxdata)
[M, N, nDof, Fmax, Fpas, ID_data, MomentArm] = deal(auxdata{:});

nA = M + nDof; % number of actuators
J = zeros(nDof*N,(nA)*N);
for i = 1:N
    for k = 1:nDof
        J(k+(i-1)*nDof,nA*(i-1)+1:nA*i) = MomentArm((i-1)*nA+1:i*nA,k)';  
    end
end

J = sparse(J);
end

% ------------------------------------------------------------------
function J = jacobian_SO (x, auxdata)
[M, N, nDof, Fmax, Fpas, ID_data, MomentArm] = deal(auxdata{:});

nA = M + nDof; % number of actuators
J = zeros(nDof*N,(nA)*N);
for i = 1:N
    for k = 1:nDof
        J(k+(i-1)*nDof,nA*(i-1)+1:nA*i) = Fmax((i-1)*nA+1:i*nA)'.*MomentArm((i-1)*nA+1:i*nA,k)'; % Partial derivatives of the constraints to the optimization variables
    end
end

J = sparse(J);
end


