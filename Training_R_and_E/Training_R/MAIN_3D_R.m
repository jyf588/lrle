clear 

% Create some dummy input data to test
% Specify the osim-model name
% This model has 23 dofs and 92 muscles.
% Alternatively, use gait2354_simbody.osim
modelName = 'gait2392_simbody.osim';
% Set up to evaluate left hip flexion, hip
% adduction, hip rotation, knee flexion and ankle plantairflexion torques.
DOF = {'hip_flexion_l' 'hip_adduction_l' 'hip_rotation_l' 'knee_angle_l' 'ankle_angle_l'};
batchsize = 200;   % number of samples per ipopt batch, does not really matter

import org.opensim.modeling.*
MAIN_path = pwd;
c = parallel.pool.Constant(@() org.opensim.modeling.Model(fullfile(MAIN_path,modelName)) );
spmd
   c.Value.initSystem();
end
auxdata.c = c;

repeat = 12000;
% generate repeat*batchsize samples
In = zeros(repeat, batchsize, 15);
% Out = zeros(repeat, batchsize, 1);
Out2 = zeros(repeat, batchsize, 5);
parfor i = 1:repeat
    q = [2.4*rand(batchsize,1) - 0.9*ones(batchsize,1) 1.2*rand(batchsize,1) - 0.6*ones(batchsize,1) ...
    1.2*rand(batchsize,1) - 0.6*ones(batchsize,1) 2.2*rand(batchsize,1) - 2.1*ones(batchsize,1)...
    1.8*rand(batchsize,1) - 0.9*ones(batchsize,1)];
    q_dot = [10*rand(batchsize,3) - 5*ones(batchsize,3), 20*rand(batchsize,1) - 10*ones(batchsize,1), 30*rand(batchsize,1) - 15*ones(batchsize,1)];
    tau = [500*rand(batchsize,1) - 250*ones(batchsize,1), 160*rand(batchsize,2) - 80*ones(batchsize,2), 500*rand(batchsize,2) - 250*ones(batchsize,2)];

    % Call the classifier, it returns a vector which indicates whether the
    % chosen sample is feasible, activations that are needed and the values of
    % the reserveactuators.
    [Torque_Manageable,Activations,ReserveActuators,MusclesOI] = TorqueLimitClassifier(c.Value,DOF,q,q_dot,tau);
    In(i,:,:) = [q q_dot tau];
%     Out(i,:,:) = Torque_Manageable;
    Out2(i,:,:) = ReserveActuators;
end

InR = reshape(In,[repeat*batchsize,15]);
OutR = reshape(Out2,[repeat*batchsize,5]);
save('InOutR_2_4M_3D_q-0915&-0606&-0606&-2101&-0909_dq555&10&15_tau250&80&80&250&250_2392','InR','OutR')

% Note: In case of infeasible torque (indicated by 0) the activations
% should max out, i.e. reach 1 (at least for some muscles) and the reserve
% actuator torques are non-trivial.