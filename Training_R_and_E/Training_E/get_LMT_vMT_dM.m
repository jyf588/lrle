function [LMT,VMT,dM] = get_LMT_vMT_dM(Model_OS,MusclesOI,CoordinatesOI,kinematics,DOF)
% Find  Muscle-Tendon Length & Velocity and the Moment-Arm Matrix

LMT = zeros(size(kinematics,1),length(MusclesOI));
VMT = zeros(size(kinematics,1),length(MusclesOI));
dM = zeros((length(MusclesOI)+length(DOF))*size(kinematics,1),length(DOF));
state = Model_OS.updWorkingState();
for t = 1:size(kinematics,1)
    
    for i = 1:length(DOF)
        stateName = cell2mat(DOF(i));
        state_value = kinematics(t,i);
        Model_OS.setStateVariable(state,stateName,state_value);
        stateName = [stateName  '_u'];
        state_value = kinematics(t,i+length(DOF));
        Model_OS.setStateVariable(state,stateName,state_value);
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
        LMT(t,i) = muscle.getLength(state);
        VMT(t,i) = muscle.getLengtheningSpeed(state);
        for k = 1:length(DOF)
            if find(CoordinatesOI(k).index == [MusclesOI(i).coordinate(:).coordinateIndex])
            coordinate = coordinateSet.get(CoordinatesOI(k).index-1);
            dM((t-1)*(length(MusclesOI)+length(DOF)) + i,k) = muscle.computeMomentArm(state, coordinate);
            else
            end
        end
    end
    for k = 1:length(DOF)
        dM((t-1)*(length(MusclesOI)+length(DOF)) + length(MusclesOI) + k,k) = 1;  % We add this part to the moment arm matrix as the moment arms of the reserve actuators (which is 1 to the joint it is actuating and 0 to the others)
    end
end



