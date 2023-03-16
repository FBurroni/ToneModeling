%% Load Table
clearvars
load("muThaiTones.mat");

%% Plot average trajectories

% Figure
f = figure("WindowState","maximized");

% Load color array
C = struct2array(load("toneColors.mat"));
M = struct2array(load("muThaiTones.mat"));

% Plot five tones
for k = 1 : 5
    plot(M(k,:),"Color",C(k,:),"LineWidth",3); hold on;
end

% Change aspect of the plot
set(gca,"FontName","Assistant","FontSize",40)
grid on
legend(["M","L","F","H","R"])

close(f)
%% Code the Task-Dynamic model

% % Choose tone
% toneIx = 1;
% 
% % Plot example of TV evolution using the tas
% t = linspace(1,5,1000);
% Y = [.2;0];
% p = [100 -4];
% 
% [t,y] = ode45(@(t,y)t
% dFunStep(t,y,p),t,Y);
% plot(t,y(:,1))
%% Create optimization problem with only stiffness and target
close all

% Choose tone
toneIx = 1;

% Set number of parameters (Target,Stiffness) and tones (Mid)
nPars = 2; nGest = 1;

% Optimization problem
p = optimvar("p",nPars,nGest,"Type","continuous","LowerBound",[-1000 -1000],"UpperBound",[1000 1000]);

% Original time
t = 1:size(M,2);

% New time
nSamp = 100;
t2 = linspace(t(1),t(end),100);

% Initial Conditions
Y = [M(toneIx,1);0];

% Prepare parameters that need to be fed to ODE
fcn = fcn2optimexpr(@ParamsToTDODE,p,t2,Y);

% Get ground truth and resample to 500 samples
groundTruth = M(toneIx,:);
groundTruth = interp1(t,groundTruth,t2);

% Define objective function
obj = sum((fcn - groundTruth).^2);

% Define optimization problem
prob = optimproblem("Objective",obj); 

% Initial guess and solve problems
p0.p = [0 0];
[psol,sumsq] = solve(prob,p0,"Solver","lsqnonlin");
pOpt = [psol.p(1),psol.p(2)];
% pOpt = [psol.p(1),-1.1];

% Get solution for optimized parameters
ySol= ParamsToTDODE(pOpt,t2,Y);

% Plot
plot(t2,groundTruth,"Color",C(toneIx,:),"LineWidth",3); hold on;
plot(t2,ySol,"o","MarkerEdgeColor",C(2,:),"MarkerFaceColor",C(2,:).*.5+[1 1 1].*.5,...
    "MarkerSize",10)

grid on
legend(["Data","Synthesized"])
set(gca,"FontName","Assistant","FontSize",20)


%% Create optimization problem with activation time

% Choose tone
toneIx = 5;

% Set number of parameters (Target,Stiffness) and tones (Mid)
nPars = 6; nGest = 2;

% Optimization problem
lb = [0 -5 t2(1)  t2(30) .1 0]';
ub = [10 -.5 t2(1) t2(60) .25 1]';
lb2 = [0 .25 t2(30) t2(end) .1 0]';
ub2 = [10 5 t2(70) t2(end) .25 1]';
p = optimvar("p",nPars,nGest,"Type","continuous","LowerBound",[lb lb2],"UpperBound",[ub,ub2]);

% Original time
t = 1:size(M,2);

% New time
nSamp = 1000;
t2 = linspace(t(1),t(end),100);

% Initial Conditions
Y = [M(toneIx,1);0];

% Prepare parameters that need to be fed to ODE
fcn = fcn2optimexpr(@ParamsToTDODE2,p,t2,Y);

% Get ground truth and resample to 500 samples
groundTruth = M(toneIx,:);
groundTruth = interp1(t,groundTruth,t2);

% Define objective function
obj = sum((fcn - groundTruth).^2);

% Define optimization problem
prob = optimproblem("Objective",obj); 

% Initial guess and solve problems
p0.p = zeros(nPars,nGest);
options = optimoptions('lsqnonlin','FiniteDifferenceStepSize',1e-8,"FiniteDifferenceType","central");
[psol,sumsq] = solve(prob,p0,"Solver","lsqnonlin","Options",options);
pOpt = struct2array(psol);

% pOpt = [.02 .01; -1.25 .1; 0 25; 24 37];

% Get solution for optimized parameters
ySol= ParamsToTDODE2(pOpt,t2,Y);

% Plot
plot(t2,groundTruth,"Color",C(toneIx,:),"LineWidth",3); hold on;
plot(t2,ySol,"o","MarkerEdgeColor",C(2,:),"MarkerFaceColor",C(2,:).*.5+[1 1 1].*.5,...
    "MarkerSize",10)

grid on
legend(["Data","Synthesized"])
set(gca,"FontName","Assistant","FontSize",20)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dYdt = tdFunStep(t,y,p)
% We need to turn the second-order ODE into a first-order lineary system
% rename velocity z' as z2, and acceleration z'' to z2'

% Position
z = y(1);

% Velocity
z2 = y(2);

% Stiffness
K = p(1);

% Target
z0=p(2);

% Damping
B = 2.*sqrt(K);

% Distance between position and target
deltaZ = z-z0;

% Inertial coefficients
M = 1;

% TV Equation
dz2dt = M .* (-B .* z2 - K .* deltaZ);

% Assemble a vector with velocity and acceleration
dYdt = [z2; dz2dt];
end

function dYdt = tdFunTime(t,y,p)
% We need to turn the second-order ODE into a first-order lineary system
% rename velocity z' as z2, and acceleration z'' to z2'

% Position
z = y(1);

% Velocity
z2 = y(2);

% Stiffness
K = p(1);
% K = ones(1,length(t)).*p(1);

% Activation
% a = zeros(1,length(t));
z0 = zeros(1,length(t))+exp(-t*p(6));

for k = 1 : size(p,2)
    z0 = z0  + getAContRamp(t,p(3,k),p(4,k),p(5,k)) .* p(2,k);
end

% Damping
B = 2.*sqrt(K);

% Target
% z0=p(2).*a;

% Distance between position and target
% deltaZ = z-z0;
% z0 = a.*z0;

% Inertial coefficients
M = 1;

% TV Equation
dz2dt = M .* (-B .* z2 - K .* (z-z0));

% Assemble a vector with velocity and acceleration
dYdt = [z2; dz2dt];
end

function a = getAStep(t,tOn,tOff)

a = zeros(1,length(t));
a(t>=tOn&t<=tOff)=1;

end


function a = getACont(t,tOn1,tOn2,tOff1,tOff2)

a = zeros(1,length(t));
a(t>tOn2&t<tOff1)=1;
a(t>tOn1&t<=tOn2)=sin(pi/2.*(t(t>tOn1&t<=tOn2)-tOn1)/(tOn2-tOn1));
a(t>=tOff1&t<=tOff2)=cos(pi/2.*(t(t>=tOff1&t<=tOff2)-tOff1)/(tOff2-tOff1));

end

function a = getAContRamp(t,tOn1,tOff2,rampDur)

tOn2 = tOn1+rampDur;
tOff1 = tOff2-rampDur;

a = zeros(1,length(t));

a(t>tOn2&t<tOff1)=1;
a(t>tOn1&t<=tOn2)=sin(pi/2.*(t(t>tOn1&t<=tOn2)-tOn1)/(tOn2-tOn1));
a(t>=tOff1&t<=tOff2)=cos(pi/2.*(t(t>=tOff1&t<=tOff2)-tOff1)/(tOff2-tOff1));

end

function solpts = ParamsToTDODE(p,t0,y0)
sol =  ode89(@(t,y)tdFunStep(t,y,p),t0,y0);
solpts = deval(sol,t0);
solpts = solpts(1,:);
end


function solpts = ParamsToTDODE2(p,t0,y0)
sol =  ode89(@(t,y)tdFunTime(t,y,p),t0,y0);
solpts = deval(sol,t0);
solpts = solpts(1,:);
end