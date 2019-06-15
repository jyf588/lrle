clear 

load('R4E_2D.mat')

% Specify the osim-model name
% This model has 23 dofs and 92 muscles.
% Alternatively, use gait2354_simbody.osim
modelName = 'gait2392_simbody.osim';
% Set up to evaluate left hip flexion, hip
% adduction, hip rotation, knee flexion and ankle plantairflexion torques.
DOF = {'hip_flexion_l' 'hip_adduction_l' 'hip_rotation_l' 'knee_angle_l' 'ankle_angle_l'};
batchsize = 200;   % number of samples per batch

qdqtau_total = valid_qdqtau;
% qdqtau_total = qdqtau_total(1:1600,:);

repeat = floor(size(qdqtau_total,1)/batchsize);
qdqtau_total = qdqtau_total(1:repeat*batchsize,:);

p_qdqtau = parallel.pool.Constant(qdqtau_total);

import org.opensim.modeling.*
MAIN_path = pwd;
% Model_OS = Model(fullfile(MAIN_path,modelName));
c = parallel.pool.Constant(@() org.opensim.modeling.Model(fullfile(MAIN_path,modelName)) );
spmd
   c.Value.initSystem();
end
auxdata.c = c;

% generate repeat*batchsize samples and store them in 2 matrices
In = zeros(repeat, batchsize, 9);
Out = zeros(repeat, batchsize, 1);
Out2 = zeros(repeat, batchsize, 3);
parfor i = 1:repeat
    q = p_qdqtau.Value((i-1)*batchsize+1 : i*batchsize,[1 1 1 2 3]);
    q_dot = p_qdqtau.Value((i-1)*batchsize+1 : i*batchsize,[4 4 4 5 6]);
    tau = p_qdqtau.Value((i-1)*batchsize+1 : i*batchsize,[7 7 7 8 9]);
    q(:,[2 3]) = 0.0;
    q_dot(:,[2 3]) = 0.0;
    tau(:,[2 3]) = 0.0;

    % Call the classifier, it returns a vector which indicates whether the
    % chosen sample is feasible, activations that are needed and the values of
    % the reserveactuators.
    [Torque_Manageable,Activations,ReserveActuators,MusclesOI,EnergyValue] = TorqueLimitClassifierEnergyObj(c.Value,DOF,q,q_dot,tau);
    In(i,:,:) = [q(:,[1 4 5]) q_dot(:,[1 4 5]) tau(:,[1 4 5])];
    Out(i,:,:) = EnergyValue;
    Out2(i,:,:) = ReserveActuators(:,[1 4 5]);
end

InR = reshape(In,[repeat*batchsize,9]);
OutE = reshape(Out,[repeat*batchsize,1]);
OutR = reshape(Out2,[repeat*batchsize,3]);
save('InOutRE_ValidR4E_2D_q-1224&-2501&-1211_dq5&15&15_tau200250250_2392.mat','InR','OutR','OutE')

