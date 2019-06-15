function gradf = jumperLeftRevoluteFootStatic_DC_ObjGradFun(X,auxdata)

import org.opensim.modeling.*

del = sqrt(eps);    % findiff step size

N        = auxdata.N;
osimModel = auxdata.model;
% h        = auxdata.h;
% Nstates  = auxdata.Nstates;
% Ncoord   = auxdata.Ncoord;

osimState = osimModel.updWorkingState(); 
numVar = osimState.getNY();
for i = 0:numVar-1
   osimState.updY().set(i, X((i+1)*N-1,1)); % use the second last vel and pos
end

fh = osimModel.calcMassCenterPosition(osimState).get(1);
fhv = osimModel.calcMassCenterVelocity(osimState).get(1);
f = -100 * (fhv^2/2/9.80665 + fh);

sv = zeros(numVar,1);
sj = zeros(numVar,1);
si = zeros(numVar,1);
for i = 1:numVar
    sj(i) = 1;
    si(i) = i*N-1;
    osimState.updY().set(i-1, X(i*N-1,1)+del);
    pfh = osimModel.calcMassCenterPosition(osimState).get(1);
    pfhv = osimModel.calcMassCenterVelocity(osimState).get(1);
    pf = -100 * (pfhv^2/2/9.80665 + pfh);
    sv(i) = (pf - f)/del;
    osimState.updY().set(i-1, X(i*N-1,1));
end
gradf = sparse(si, sj, sv, length(X), 1);

end