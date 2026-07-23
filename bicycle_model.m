function out = bicycle_model(mF, mR, V, delta, beta, mass, l, bPCT)
%function out = bicycle_model(m, V, delta, beta, mass, l, bPCT)
% m = magicformula object (for evaluating tire forces & moments)
% V = Vehicle Speed [m/s]
% delta = steer angle [rad]
% beta = chassis slip angle [rad]
% mass = vehicle mass [kg]
% l = wheelbase [m]
% bPCT = percentage of vehicle weight on the front axle

    r = 0; % initial yaw velocity estimate
    ay = [10,0]; % initial lateral acceleration estimate
    
    % solver
    rtol = 1e-4; % iteration tolerance
    rel = 0.7; % relaxation parameter
    e = 10; % initial error
    i = 1; % iteration counterThe
    while e > rtol
        % update
        i = i + 1;
    
        % calculate lateral force
        [YF, YR, a, b] = slips_forces(mF, mR, r, V, delta, beta, mass, l, bPCT);
        [Y, N] = Y_N(a, b, YF, YR);
    
        % calculate lateral force & yaw rate
        ay(i) = rel * Y/mass + (1-rel)*ay(i-1);
        r = ay(i)/V;
    
        % calculate error
        e = abs(ay(i) - ay(i-1));
    end
    
    % output
    out.Y = Y; % lateral force sum
    out.N = N; % moment sum
    out.ay = ay(end); % lateral acceleration
    out.r = r; % yaw rate
    out.R = V^2/ay(end); % corner radius
end

function [Y,N] = Y_N(a, b, YF, YR)
    % Sum lateral force
    Y = YF + YR;

    % Sum moments
    N = a*YF - b*YR;
end

function [YF, YR, a, b, aF, aR] = slips_forces(...
    mF, mR, r, V, delta, beta, mass, l, bPCT)

    sr = 0; % slip ratio
    ia = 0; % inclination angle
    p = 70e3; % pressure (10 psi)

    % Calculate normal forces
    [a, b, fzF, fzR] = fz(mass,l,bPCT);

    % Calculate slip ratios
    aF = beta + a*r/V - delta;
    aR = beta - b*r/V;

    % Calculate lateral forces
    YF = mF.fy(aF, sr, ia, fzF, p);
    YR = mR.fy(aR, sr, ia, fzR, p);
end

function [a, b, fzF, fzR] = fz(mass,l,bPCT)
    b = bPCT*l; % distance to front axle
    a = l-b; % distance to rear axle
    ww = mass * 9.81; % weight
    fzF = ww * b/l; % front normal force
    fzR = ww - fzF; % rear normal force
end