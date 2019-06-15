clear 

% Create some dummy input data to test
% Specify the osim-model name
% This model has 23 dofs and 92 muscles.
% Alternatively, use gait2354_simbody.osim
modelName = 'gait2392_simbody.osim';
% Set up to evaluate left hip flexion, hip
% adduction, hip rotation, knee flexion and ankle plantairflexion torques.
DOF = {'hip_flexion_l' 'hip_adduction_l' 'hip_rotation_l' 'knee_angle_l' 'ankle_angle_l'};
batchsize = 300;   % number of samples per ipopt batch, does not really matter

import org.opensim.modeling.*
MAIN_path = pwd;
c = parallel.pool.Constant(@() org.opensim.modeling.Model(fullfile(MAIN_path,modelName)) );
spmd
   c.Value.initSystem();
end
auxdata.c = c;

repeat = 1000;
% generate repeat*batchsize samples
In = zeros(repeat, batchsize, 9);
% Out = zeros(repeat, batchsize, 1);
Out2 = zeros(repeat, batchsize, 3);
parfor i = 1:repeat
    q = [3.6*rand(batchsize,1) - 1.2*ones(batchsize,1), zeros(batchsize,1), ...
    zeros(batchsize,1), 2.6*rand(batchsize,1) - 2.5*ones(batchsize,1), ...
    2.3*rand(batchsize,1) - 1.2*ones(batchsize,1)];
    q_dot = [10*rand(batchsize,1) - 5*ones(batchsize,1), zeros(batchsize,2), 30*rand(batchsize,2) - 15*ones(batchsize,2)];
    tau = [400*rand(batchsize,1) - 200*ones(batchsize,1), zeros(batchsize,2), 500*rand(batchsize,2) - 250*ones(batchsize,2)];

    % Call the classifier, it returns a vector which indicates whether the
    % chosen sample is feasible, activations that are needed and the values of
    % the reserveactuators.
    [Torque_Manageable,Activations,ReserveActuators,MusclesOI] = TorqueLimitClassifier(c.Value,DOF,q,q_dot,tau);
    In(i,:,:) = [q(:,[1 4 5]) q_dot(:,[1 4 5]) tau(:,[1 4 5])];
%     Out(i,:,:) = Torque_Manageable;
    Out2(i,:,:) = ReserveActuators(:,[1 4 5]);
end

InR = reshape(In,[repeat*batchsize,9]);
OutR = reshape(Out2,[repeat*batchsize,3]);
save('InOutR_2D_300k_q-1224&-2501&-1211_dq5&15&15_tau200250250_2392','InR','OutR')

% Note: In case of infeasible torque (indicated by 0) the activations
% should max out, i.e. reach 1 (at least for some muscles) and the reserve
% actuator torques are non-trivial.