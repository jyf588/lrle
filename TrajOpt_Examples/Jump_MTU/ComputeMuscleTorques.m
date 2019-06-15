function [tau, dTaudA] = ComputeMuscleTorques(Model_OS, state, MusclesOI, CoordinatesOI, DOF, a, q, qdot)
% just make sure state is initialized before first call of
% ComputeMuscleTorques
% Inputs:
% DOF   cell array of strings of length nDof containing the names of the
%       degrees of freedom that are considered.
% a     length(MuscleOI)*1 activations
% q     nDof*1 matrix of generalized coordinates (radians, meters)
% qdot  nDof*1 matrix of generalized velocities
% Outputs:
% tau   nDof*1 matrix of torques
% dTaudM  f0*fl*fv*cos(alpha) for each muscle * dM(M,nDoF) nDoF*M dTaudM

import org.opensim.modeling.*

M = length(MusclesOI);
nDof = length(DOF);

kinematics = [q; qdot];

% Now we calculate based on the kinematic state the Muscle-Tendon Length & Velocity and the Moment-Arm Matrix
[LMT,VMT,dM] = get_LMT_vMT_dM_my(Model_OS,state, MusclesOI,CoordinatesOI,kinematics,DOF);

% Initialize some matrices
FMltilde = ones(M,1); % Normalized Force-Length dependence
FMvtilde = ones(M,1); % Normalized Force-Velocity dependence
Fpe = ones(M,1);      % Passive (elastic) force
cos_alpha = ones(M,1); % cosine of the muscles pennation angle

% We apply the Hill Model under the assumption of a rigid tendon and
% calculate the normalized quantities of passive and active muscle force
% components. The active component will be scaled with the muscle
% activation and added to the passive force to produce the final muscle force.
for m = 1:M
    [~, ~, FMltilde(m), FMvtilde(m), Fpe(m), cos_alpha(m)] = HillModel_RigidTendon_my(1.0, LMT(m),VMT(m),MusclesOI(m));
end
FMo = reshape([MusclesOI(:).maxIsoForce], M, 1);
Fpas = FMo.*Fpe.*cos_alpha;
Fact = FMo.*FMltilde.*FMvtilde.*cos_alpha;

FM = a.* Fact + Fpas;
tau = dM.'*FM;
dTaudA = (dM.').* repmat(Fact.',nDof,1); 
end



