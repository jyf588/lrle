TF = 1.0;
import org.opensim.modeling.*
osimModel = Model('LeftSideModel_2D_torque_revolute_new.osim');
osimState = osimModel.initSystem();
dc_time = linspace(0,TF,101).';
dt = dc_time(2) - dc_time(1);
Nstates       = osimModel.getNumStateVariables();
StatesData = struct();
StatesData.name = [char(osimModel.getName()), '_Initial_States_halfCronch'];
StatesData.nRows = size(dc_time, 1);
StatesData.nColumns = Nstates+1; %All the states + time
StatesData.inDegrees = false;
StatesData.labels = cell(1,StatesData.nColumns); 
StatesData.labels{1}= 'time';

X_state_init = zeros(StatesData.nRows, Nstates);

% toeaux, ankle, knee, hip

start_pose = [0 48.158 -104.737 97.895] / 180 * pi;
end_pose = start_pose;
mid_pose = start_pose;

for i = 1:Nstates/2
    X_state_init(:,i) = [linspace(start_pose(i), mid_pose(i), 40).'; linspace(mid_pose(i),end_pose(i), StatesData.nRows-40).'];
    X_state_init(:,i+Nstates/2) = [0; (X_state_init(2:end,i) - X_state_init(1:end-1,i)) / dt;];
end

X_state_init = X_state_init + randn(StatesData.nRows, Nstates) * 0.002; % only pertube a little bit


for j = 2:1:StatesData.nColumns
   StatesData.labels{j} = osimModel.getStateVariableNames().getitem(j-2);
end
StatesData.data = [dc_time, X_state_init];
writeOpenSimStatesFile(StatesData)