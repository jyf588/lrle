clear 

load('R4E_3D.mat')

% Specify the osim-model name
% This model has 23 dofs and 92 muscles.
% Alternatively, use gait2354_simbody.osim
modelName = 'gait2392_simbody.osim';
% Set up to evaluate left hip flexion, hip
% adduction, hip rotation, knee flexion and ankle plantairflexion torques.
DOF = {'hip_flexion_l' 'hip_adduction_l' 'hip_rotation_l' 'knee_angle_l' 'ankle_angle_l'};
batchsize = 200;   % number of samples per batch

%TODO: Needs loading
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

% calc E for the repeat*batchsize input samples
In = zeros(repeat, batchsize, 15);
Out = zeros(repeat, batchsize, 1);
Out2 = zeros(repeat, batchsize, 5);
% Acts = zeros(repeat, batchsize, 43);
parfor i = 1:repeat
    q = p_qdqtau.Value((i-1)*batchsize+1 : i*batchsize,1:5);
    q_dot = p_qdqtau.Value((i-1)*batchsize+1 : i*batchsize,6:10);
    tau = p_qdqtau.Value((i-1)*batchsize+1 : i*batchsize,11:15);
    % Call the classifier, it returns a vector which indicates whether the
    % chosen sample is feasible, activations that are needed and the values of
    % the reserveactuators.
    [Torque_Manageable,Activations,ReserveActuators,MusclesOI,EnergyValue] = TorqueLimitClassifierEnergyObj(c.Value,DOF,q,q_dot,tau);
    In(i,:,:) = [q q_dot tau];
    Out(i,:,:) = EnergyValue;
    Out2(i,:,:) = ReserveActuators;
%     Acts(i,:,:) = Activations;  % not really used
end

InR = reshape(In,[repeat*batchsize,15]);
OutE = reshape(Out,[repeat*batchsize,1]);
OutR = reshape(Out2,[repeat*batchsize,5]);
% ActsE = reshape(Acts,[repeat*batchsize,43]);
save('InOutRE_ValidR4E_3D_q-0915&-0606&-0606&-2101&-0909_dq555&10&15_tau250&80&80&250&250_2392_simple_volumn_E','InR','OutR','OutE')

