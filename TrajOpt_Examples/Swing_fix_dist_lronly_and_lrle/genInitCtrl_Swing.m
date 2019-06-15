TF = 3.0;
import org.opensim.modeling.*
osimModel = Model('swing_model.osim');
osimState = osimModel.initSystem();
dc_time = linspace(0,TF,201).';
Ncontrols = osimModel.getNumControls();
ControlData = struct();
ControlData.name = [char(osimModel.getName()), '_Initial_Controls_3.0s'];
ControlData.nRows = size(dc_time, 1);
ControlData.nColumns = Ncontrols+1; %All the controls + time
ControlData.inDegrees = false;
ControlData.labels = cell(1,ControlData.nColumns); 
ControlData.labels{1}= 'time';

X_controls_init = ones(ControlData.nRows, Ncontrols) * 0.02;

for j = 2:1:ControlData.nColumns
   ControlData.labels{j} = osimModel.getActuators().get(j-2).getName();
end
ControlData.data = [dc_time, X_controls_init];
writeOpenSimControlFile(ControlData)