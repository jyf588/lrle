% ----------------------------------------------------------------------- %
% JumperLeftRevoluteFootStatic_DC.m
%
% ----------------------------------------------------------------------- %

clear

global hisy;
hisy = [];

% TODO: different from the joint limits defined in osim 
% [toe_aux ankle knee hip]
qLB = [-0.8 -1.0 -2.5 -0.9];
qUB = [0.1 1.1 0.1 2.2];
vB = [20 15 15 5];

init_q = [0 48.158 -104.737 97.895] / 180 * pi;
auxdata.init_q = init_q;
init_a = 0.1;
tor_limit = [250 250 200]; % ankle knee hip

% Node and timing information
N = 200;                 % # of nodal points
duration = 1.0;         % time in sec to complete task
h = duration/(N-1);     % time interval between nodes
dc_time = h*(0:N-1)';   % list of time points (temporal grid)
% Import the OpenSim modeling classes
import org.opensim.modeling.*

model_name = 'LeftSideModel_2D_torque_revolute_new.osim';

% Read in the osim model
osimModel = Model(model_name); % 
% Initialize the model (this builds the system and initialize the state)
osimState = osimModel.initSystem();

c = parallel.pool.Constant(@() org.opensim.modeling.Model(model_name) );
spmd
   c.Value.initSystem();
end
auxdata.c = c;

% Get the number of states, coordinates, muscles and controls from the model;
% in this case the number of controls equals the number of muscles
Nstates       = osimModel.getNumStateVariables();
Ncontrols     = osimModel.getNumControls();
Ncoord        = osimModel.getNumCoordinates(); 
% model_muscles = osimModel.getMuscles();
% Nmuscles      = model_muscles.getSize();


% Auxiliary data to be passed to the optimizer
auxdata.model      = osimModel;
auxdata.time       = dc_time;
auxdata.N          = N;
auxdata.h          = h;
auxdata.Nstates    = Nstates;
auxdata.Ncontrols  = Ncontrols; 
auxdata.Ncoord     = Ncoord;
% auxdata.Nmuscles   = Nmuscles;
auxdata.tor_limit     = tor_limit;


if (Ncontrols ~= 3)
    disp('Ncontrols should be 3 for lower limb');
    return
end

if (Ncoord ~= 4)
    disp('Ncoord should be 4 for lower limb revolute');
    return
end

% Get the names of the states from the model
states_all = cell(Nstates,1);
for i = 1:Nstates
   states_all(i,1) = cell(osimModel.getStateVariableNames().getitem(i-1));
end

% Get the names of the controls/muscles from the model (same in this case) 
controls_all = cell(Ncontrols,1);
for i = 1:Ncontrols
   controls_all(i,1) = cell(osimModel.getActuators().get(i-1).getName());
end

%-----------------------------------------------------------------------------
% 
% In the LowerLimbTorque_Revolute example, the states are:
% 8 states (q(1,2,3,4), dq(1,2,3,4))
% the controls are:
% tau(1,2,3)
%-----------------------------------------------------------------------------

% Load the .sto file that contain the initial guess for the states
[file_input, pathname] = uigetfile({'*.sto', 'OpenSim States Files (*.sto)'}, ...
                         'Select the initial states file','MultiSelect', 'off');
temp_s = importdata(strcat(pathname,file_input)); % import states data
old_time_s = temp_s.data(:,1);         % time 
old_data_s = temp_s.data(:,2:end);     % the 8 states
old_text_s = temp_s.textdata(7,2:end); % states names

% Arrange the initial guess by nodes and states
x0_temp = zeros(N,Nstates);  % pre-allocate space
chk_counter =0;
for j = 1:size(states_all,1) 
    for k = 1:size(old_data_s,2)
        if strcmp(old_text_s(k),states_all(j)) == 1
            % interpolate initial guess to the defined temporal grid
            dc_time(dc_time > old_time_s(end) & dc_time < (old_time_s(end)+1e-6)) = old_time_s(end);
            x0_temp(:,j) = interp1(old_time_s,old_data_s(:,k),dc_time);
            chk_counter = chk_counter+1;
        end
    end
end
if (chk_counter ~= size(states_all,1))
    disp('Inconsistent number of states for inital guess!');
    disp('Check the input file...');
    return
end

X0 = zeros((Ncontrols+Nstates)*N,1);
% Arrange the initial guess into a column vector
% [ [pos(t_0) ... pos(t_N)], [vel(t_0) ... vel(t_N)], ... etc]
for i = 1:Nstates
    X0(N*(i-1)+1:N*i,1) = x0_temp(:,i);
end


% Load the file that contain the initial guess for the controls (excitations)
[file_input, pathname] = uigetfile({'*.sto', 'OpenSim Controls (excitation) Files (*.sto)'}, ...
                                   'Select the initial controls file','MultiSelect', 'off');
temp_i = importdata(strcat(pathname,file_input)); % import controls data (1st column is time)
old_time_i = temp_i.data(:,1);         % time 
old_data_i = temp_i.data(:,2:end);     % the controls
old_text_i = temp_i.textdata(7,2:end); % controls names


% Arrange the initial guess by nodes and controls
u0_temp = zeros(N,Ncontrols);
for j = 1:size(controls_all,1)
    for k = 1:size(old_data_i,2)
        if strcmp(old_text_i(k),controls_all(j)) == 1
            % interpolate to the temporal grid
            dc_time(dc_time > old_time_i(end) & dc_time < (old_time_i(end)+1e-6)) = old_time_i(end);
            u0_temp(:,j) = interp1(old_time_i,old_data_i(:,k),dc_time);
        end
    end
end

% Append the initial guess for the controls to the end of X0
for i = 1:Ncontrols
    X0(Nstates*N + N*(i-1)+1 : Nstates*N + N*i, 1) = u0_temp(:,i);
end


% Check: make sure both files (states and controls) are consistent
if (old_time_i(end) ~= old_time_s(end) || old_time_i(end) ~= duration)
    disp('Time stamp of states and controls for initial guess do not match!');
    disp('Check the input files...');
    return
end

%------------------------------
% Initialize Muscle computations (used to calculate muscle paths, lm,vm only)
%---------------------------------
mus_model_name = 'gait2392_simbody.osim';
% Read in the osim model
musModel = Model(mus_model_name);
musModel.initSystem();
% DOF = {'hip_flexion_l' 'hip_adduction_l' 'hip_rotation_l' 'knee_angle_l' 'ankle_angle_l'};
DOF = {'ankle_angle_l' 'knee_angle_l' 'hip_flexion_l'};
[MusclesOI,CoordinatesOI] = getMuscleStructure(musModel,DOF);
auxdata.musModel      = musModel;
auxdata.DOF           = DOF;
auxdata.MusclesOI     = MusclesOI;
auxdata.CoordinatesOI = CoordinatesOI;

% parallel model replicates
c2 = parallel.pool.Constant(@() org.opensim.modeling.Model(mus_model_name) );
spmd
   c2.Value.initSystem();
end
auxdata.c2 = c2;

% append activations to X0 (a1(1~N), a2(1~N), ... a_24(1~N))
X0 = [X0; init_a*ones(length(MusclesOI)*N,1)];

%-----------------------------------------------------------------
% Set-up and solve the nonlinear programming problem
%-----------------------------------------------------------------

% Bounds on the optimization parameters

for i = 1:Ncoord
    Pos_LB(N*(i-1)+1:N*i) = qLB(i);
    Pos_UB(N*(i-1)+1:N*i) = qUB(i);
end
for i = 1:Ncoord
    Vel_LB(N*(i-1)+1:N*i) = -vB(i);
    Vel_UB(N*(i-1)+1:N*i) = vB(i);
end

Con_LB(1:Ncontrols*N) = -1.0;    Con_UB(1:Ncontrols*N) = 1.0;

% append activations
Act_LB(1:length(MusclesOI)*N) = 0.0; Act_UB(1:length(MusclesOI)*N) = 1.0;

options.lb = [Pos_LB Vel_LB Con_LB Act_LB]';
options.ub = [Pos_UB Vel_UB Con_UB Act_UB]';

options.cl = [zeros((N-1)*Nstates + Nstates,1); ones(N-1,1)*(-Inf)]; % first (N-1)*Nstates + Nstates + floor(N/10*7) ceq, N-1 cn (<=0)
options.cu = [zeros((N-1)*Nstates + Nstates,1); zeros(N-1,1)]; % Upper bounds on the constraints (comply with the torques that are provided)
options.cl(end, 1) = 0; % last one should be exactly zero contact force
% add muscle equality constr here
options.cl = [options.cl; zeros(Ncontrols*N,1)];
options.cu = [options.cu; zeros(Ncontrols*N,1)];

options.auxdata = auxdata;

options.ipopt.mu_strategy      = 'adaptive';
options.ipopt.max_iter         = 250;
% options.ipopt.derivative_test  = 'first-order';
% options.ipopt.derivative_test_perturbation  = sqrt(eps);
options.ipopt.nlp_scaling_method = 'gradient-based';
options.ipopt.nlp_scaling_max_gradient = 50;
options.ipopt.tol              = 1e-1;
options.ipopt.hessian_approximation = 'limited-memory';
options.ipopt.limited_memory_max_history = 500;
% options.ipopt.print_level = 2;
% options.ipopt.print_frequency_iter = 2;
options.ipopt.acceptable_compl_inf_tol = 2.0;
options.ipopt.acceptable_constr_viol_tol = 1.0;
options.ipopt.acceptable_obj_change_tol = 1e-3;
options.ipopt.acceptable_tol = 1.0;
options.ipopt.acceptable_iter = 10;

% The callback functions.
funcs.objective         = @jumperLeftRevoluteFootStatic_DC_ObjFun;
funcs.constraints       = @jumperLeftRevoluteFootStatic_DC_ConFun;
funcs.gradient          = @jumperLeftRevoluteFootStatic_DC_ObjGradFun;
funcs.jacobian          = @jumperLeftRevoluteFootStatic_DC_ConJacFun;
funcs.jacobianstructure = @jumperLeftRevoluteFootStatic_DC_ConJacStructFun;
% funcs.iterfunc          = @iterfunc;

% start a timer
tic;

[Xopt,~] = ipopt_auxdata(X0,funcs,options);



% Define anonymous functions for the objective and constraint functions
% objfun = @(x)jumperLeftRevoluteFoot_DC_ObjFun(x,auxdata);
% confun = @(x)jumperLeftRevoluteFoot_DC_ConFun(x,auxdata);


% % Set some parameters for fmincon and then run the optimization
% % TODO; let us test Useparallel later. 
% options = optimset('algorithm','interior-point','TolFun',1e-4,'TolX',1e-4, ...
%                    'TolCon',1e-4,'FinDiffType','forward','MaxFunEvals',1e5, ...
%                    'Hessian','bfgs','display','iter','GradConstr','on','DerivativeCheck','off',... 
%                     'FinDiffRelStep',sqrt(eps));
                
% Set some parameters for fmincon and then run the optimization
% TODO; let us test Useparallel later. 
% options = optimset('algorithm','interior-point','TolFun',5e-4,'TolX',5e-4, ...
%                    'TolCon',1e-2,'FinDiffType','forward','MaxFunEvals',1e7, 'MaxIter',1e5, ...
%                    'Hessian','bfgs','display','iter','GradConstr','on', 'GradObj','on',...
%                    'DerivativeCheck','off', 'FinDiffRelStep',sqrt(eps),'OutputFcn',@myoutput);



% [Xopt,fval,exitflag,output] = fmincon(objfun,X0,[],[],[],[],lb,ub,confun,options);

% stop the timer
runtime = toc;


% Store activation values
act = reshape(Xopt((Ncontrols+Nstates)*N+1:end), N, length(MusclesOI));
save(['activations',char(date),'_1.0s_252520MusStatic2392.mat'],'act')

%------------------------------------------------------------
% Print the optimal states and controls to STO files
%------------------------------------------------------------

% Put Xopt back in 'row = time, column = states/controls' format

X_state_opt = zeros(N,Nstates); %pre-allocate size
for i = 1:Nstates
    X_state_opt(:,i) = Xopt(N*(i-1)+1:N*i,1);
end

X_controls_opt = zeros(N,Ncontrols); %pre-allocate size
for i = 1:Ncontrols
    X_controls_opt(:,i) = Xopt(Nstates*N + N*(i-1)+1:Nstates*N + N*i,1);
end

% Create data structure for the states file
StatesData = struct();
StatesData.name = [char(osimModel.getName()), '_Optimal_States_', char(date),'_1.0s_252520MusStatic2392'];
StatesData.nRows = size(dc_time, 1);
StatesData.nColumns = Nstates+1; %All the states + time
StatesData.inDegrees = false;
StatesData.labels = cell(1,StatesData.nColumns); 
StatesData.labels{1}= 'time';
for j = 2:1:StatesData.nColumns
   StatesData.labels{j} = char(states_all(j-1));
end
StatesData.data = [dc_time, X_state_opt];
writeOpenSimStatesFile(StatesData)

% Create data structure for the controls file0 
ControlData = struct();
ControlData.name = [char(osimModel.getName()), '_Optimal_Controls_', char(date),'_1.0s_252520MusStatic2392'];
ControlData.nRows = size(dc_time, 1);
ControlData.nColumns = Ncontrols+1; %All the controls + time
ControlData.inDegrees = false;
ControlData.labels = cell(1,ControlData.nColumns); 
ControlData.labels{1}= 'time';
for j = 2:1:ControlData.nColumns
   ControlData.labels{j} = char(controls_all(j-1));
end
ControlData.data = [dc_time, X_controls_opt];
writeOpenSimControlFile(ControlData)
%------------------------------------------------------------

% Print out some details about the optimization run
disp( ['Optimization: elapsed time  = ' num2str(runtime) ' s'])
% disp( ['Optimization: min fun value = ' num2str(fval)])
% disp( ['Optimization: # iterations  = ' num2str(output.iterations)])
% disp( ['Optimization: # func evals  = ' num2str(output.funcCount)])
disp('   ')

% function cont = iterfunc(itnum,obj,data)
%     global hisy
%     cont = true;
%     hisy = [hisy; reshape(data.x,1,[])];
% end

