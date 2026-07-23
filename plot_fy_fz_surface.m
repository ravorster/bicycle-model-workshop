% Slip angle and Normal Force Arrays
sa = (-6:0.1:6) * pi/180;
fz = 200:200:2200;
Lsa = length(sa);
Lfz = length(fz);
sr = 0;
ia = 0;
p = 70e3; % 10psi

% Instantiate magic formula object
m = magicformula('Hoosier 43100 18.0x6.0-10 R20 7in Rim.mat');

% Calculate Fy
fy = zeros(Lsa, Lfz);
for i = 1:Lsa
    for j = 1:Lfz
        fy(i,j) = m.fy(sa(i), sr, ia, fz(j), p);
    end
end

% Plot as as surface
f = figure(1);
tl = tiledlayout(f, 1,2);
ax1 = nexttile(tl);
surf(ax1, fz,sa,fy)
xlabel(ax1, 'Normal Force [N]')
ylabel(ax1, 'Slip Angle [rad]')
zlabel(ax1, 'Lateral Force [N]')

ax2 = nexttile(tl);
surf(ax2, fz,sa,fy./fz)
xlabel(ax2, 'Normal Force [N]')
ylabel(ax2, 'Slip Angle [rad]')
zlabel(ax2, 'Lateral \mu=fy/fz [-]')
% view(ax2, [0,-1,0])
% zlim(ax2, [0.75, 0.85])