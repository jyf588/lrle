function f = jumperLeftRevoluteFootStatic_DC_ObjFun(X,auxdata)
%
% Computes the objective function value for the  simpleLowerLimb
%	X = current set of optimization parameters
%   auxdata = extra parameters and pointers to
%             instantiated OpenSim objects
%   f = sum of squared muscle activation integrals
% extract the nesessary auxiliary data

% Import the OpenSim modeling classes
import org.opensim.modeling.*

del = sqrt(eps);    % findiff step size

N        = auxdata.N;
osimModel = auxdata.model;
% h        = auxdata.h;
% Nstates  = auxdata.Nstates;
% Ncoord   = auxdata.Ncoord;

% if(osimModel.getWorkingState().getNY() == 0)
%    osimState = osimModel.initSystem();
% else
   osimState = osimModel.updWorkingState(); 
% end

numVar = osimState.getNY();
for i = 0:numVar-1
   osimState.updY().set(i, X((i+1)*N-1,1)); % use the second last vel and pos
end
% osimModel.computeStateVariableDerivatives(osimState); 

fh = osimModel.calcMassCenterPosition(osimState).get(1);
fhv = osimModel.calcMassCenterVelocity(osimState).get(1);
f = -100 * (fhv^2/2/9.80665 + fh);

global hisy
hisy = [hisy; reshape(X,1,[])];

end

