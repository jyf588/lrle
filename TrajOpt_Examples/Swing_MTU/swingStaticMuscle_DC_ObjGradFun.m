function gradf = swingStaticMuscle_DC_ObjGradFun(X,auxdata)

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

fvx = osimModel.calcMassCenterVelocity(osimState).get(0);
fvy = osimModel.calcMassCenterVelocity(osimState).get(1);
f = 100 * calcObj(fvx, fvy);

sv = zeros(numVar,1);
sj = zeros(numVar,1);
si = zeros(numVar,1);
for i = 1:numVar
    sj(i) = 1;
    si(i) = i*N-1;
    osimState.updY().set(i-1, X(i*N-1,1)+del);
    pfvx = osimModel.calcMassCenterVelocity(osimState).get(0);
    pfvy = osimModel.calcMassCenterVelocity(osimState).get(1);
    pf = 100 * calcObj(pfvx, pfvy);
    sv(i) = (pf - f)/del;
    osimState.updY().set(i-1, X(i*N-1,1));
end
gradf = sparse(si, sj, sv, length(X), 1);

    function obj = calcObj(vx,vy)
        if (vx < 0 && vy > 0)
            % this is what we want
            obj = vx * vy; % since we are minimizing; up to a scale
        end
        % shortest l2 distance to desired convex region
        if (vx < 0 && vy <= 0)
            obj = -vy;
        end
        if (vx >= 0 && vy > 0)
            obj = +vx;
        end
        if (vx >= 0 && vy <= 0)
            obj = sqrt(vx^2+vy^2);
        end
    end
end