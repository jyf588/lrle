TF = 1.0;
import org.opensim.modeling.*
osimModel = Model('LeftSideModel_2D_torque_revolute_new.osim');
osimState = osimModel.initSystem();
dc_time = linspace(0,TF,101).';
Ncontrols = osimModel.getNumControls();
ControlData = struct();
ControlData.name = [char(osimModel.getName()), '_Initial_Controls'];
ControlData.nRows = size(dc_time, 1);
ControlData.nColumns = Ncontrols+1; %All the controls + time
ControlData.inDegrees = false;
ControlData.labels = cell(1,ControlData.nColumns); 
ControlData.labels{1}= 'time';
% ControlData.labels{2}= 'hip_flexion_torq';
% ControlData.labels{3}= 'knee_angle_torq';
% ControlData.labels{4}= 'ankle_angle_torq';
X_controls_init = ones(ControlData.nRows, Ncontrols) * -0.02;
X_controls_init(:,2) = X_controls_init(:,2) * -1.0;

for j = 2:1:ControlData.nColumns
   ControlData.labels{j} = osimModel.getActuators().get(j-2).getName();
end
ControlData.data = [dc_time, X_controls_init];
writeOpenSimControlFile(ControlData)