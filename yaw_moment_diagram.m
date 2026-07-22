clearvars;clc;
kph2mps = 1000/3600;

% specify speed, delta, beta, mass
V = 54 * kph2mps; % Speed [m/s]
delta = (-12:1:12) * pi/180; % steer angle [rad]
beta = (-4.4:0.4:4.4) * pi/180; % chassis slip angle [rad]
mass = 200; % kg

% initialize tire model
mF = magicformula('Hoosier 43100 18.0x6.0-10 R20 7in Rim.mat');
mR = magicformula('Hoosier 43100 18.0x6.0-10 R20 7in Rim.mat');
mF.coeffs.SCALING_COEFFICIENTS.LKY = mF.coeffs.SCALING_COEFFICIENTS.LKY * 1.05;
mR.coeffs.SCALING_COEFFICIENTS.LKY = mR.coeffs.SCALING_COEFFICIENTS.LKY * 0.95;

% vehicle geometry
l = 1.5; % wheelbase [m]
wt_pct_front = 0.40; % percentage of weight on the front axle

% calculate for array of delta and beta
Ld = length(delta);
Lb = length(beta);
for i = Ld:-1:1
    for j = Lb:-1:1
        results = bicycle_model(mF, mR, V, delta(i), beta(j), mass, l, wt_pct_front);
        Y(i,j) = results.Y;
        N(i,j) = results.N;
    end
end

%% Plot
f = figure(2);
tl = tiledlayout(f,1,1);
ax = nexttile(tl);
ax.NextPlot = 'add';

% Non-dimensionalize values
Y = Y/(mass * 9.81);
N = N/(mass * 9.81 * l);

% Plot constant slip
for i = 1:Ld
    ps = plot(ax, Y(i,:), N(i,:), 'b');
end

% Plot constant steer
for j = 1:Lb
    pb = plot(ax, Y(:,j), N(:,j), 'r');
end

% Make the plots pretty
grid(ax, 'on')
xlabel(ax, 'Lateral Acceleration [G]')
ylabel(ax, 'Yaw Moment Coefficient [-]')
title(ax, 'Yaw Moment Diagram')
legend(ax, [ps, pb], 'Constant Slip', 'Constant Steer')