function x_dot = computeOpenSimModelXdot(states,controls,t,osimModel,osimState)
% This function sets the model to a particular state, applies the controls
% (excitations) and returns the derivaties of the state variables.
%
% Author: Brian Umberger, UMass Amherst
%
% Note: This function draws from the function OpenSimPlantFunction.m by Daniel
% Jacobs (Stanford), which is part of the OpenSim dynamic walking example
%

% Import the OpenSim modeling classes
import org.opensim.modeling.*

% Update model state with current values  
% osimState.setTime(t);
numVar = osimState.getNY();
for i = 0:numVar-1
   osimState.updY().set(i, states(i+1,1));
end

% Update the state velocity calculations
osimModel.computeStateVariableDerivatives(osimState);

% Determine number of muscles/controls in the model
% muscles = osimModel.getMuscles(); 
% nMuscles = muscles.getSize();
Ncontrols     = osimModel.getNumControls();

% Get a reference to the current model controls
modelControls = osimModel.updControls(osimState);

% Initialize a vector for the actuator controls (muscle excitations)
actuatorControls = Vector(1, 0.0);

% Set Actuator Controls and update model controls with the new values
for i = 0:Ncontrols-1
    actuatorControls.set(0, controls(i+1,1));
    osimModel.updActuators().get(i).addInControls(actuatorControls,modelControls);
end

% Set the control (excitation) values
osimModel.setControls(osimState,modelControls);

% Update the derivative calculations in the state variable
osimModel.computeStateVariableDerivatives(osimState);

% Set output variable to the new state
x_dot = zeros(numVar,1);

for i = 0:1:numVar-1
   x_dot(i+1,1) = osimState.getYDot().get(i);
end


end
