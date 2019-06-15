function [LMT,VMT,dM] = get_LMT_vMT_dM_my(Model_OS, state, MusclesOI,CoordinatesOI,kinematics,DOF)
% Find  Muscle-Tendon Length & Velocity and the Moment-Arm Matrix
% kinematics: [q;q_dot] 2nDof * 1

LMT = zeros(length(MusclesOI),1);
VMT = zeros(length(MusclesOI),1);
dM = zeros(length(MusclesOI),length(DOF));

% state = Model_OS.initSystem; % Seems not necessatry

for i = 1:length(DOF)
    stateName = cell2mat(DOF(i));
    state_value = kinematics(i);
    Model_OS.setStateVariable(state,stateName,state_value);
    stateName_u = [stateName  '_u'];
    state_value = kinematics(i+length(DOF));
    Model_OS.setStateVariable(state,stateName_u,state_value);
end
% We need to realize the multibody system (which is part of SimBody and
% is thus not part of the API exposed to Matlab). We need an opensim
% function that realizes velocity without possibly violating any of the
% constraints in the model, as this would induce an error.
Model_OS.computeStateVariableDerivatives(state);

%Model_OS.equilibrateMuscles(state); % Can't use this as some kinematic states might induce errors in the OpenSim source code.

% Now that the model is in the correct state we can read out the muscle
% information

coordinateSet = Model_OS.getCoordinateSet;
for i = 1:length(MusclesOI)
    muscle = Model_OS.getMuscles.get(MusclesOI(i).index-1);
    LMT(i) = muscle.getLength(state);
    VMT(i) = muscle.getLengtheningSpeed(state);
    for k = 1:length(DOF)
        if find(CoordinatesOI(k).index == [MusclesOI(i).coordinate(:).coordinateIndex])
            coordinate = coordinateSet.get(CoordinatesOI(k).index-1);
            dM(i,k) = muscle.computeMomentArm(state, coordinate);
        end
    end
end