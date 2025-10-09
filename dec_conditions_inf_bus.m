%% Stability Conditions for GFM converter connected to an infinite bus
%% Setup
% Frequency vector
w = logspace(log10(1e-3),log10(1e4),1000)*2*pi;
f = w/2/pi;
% s operator
s = tf('s');
%% Extract converter and network models
Yc = GmDssCell{2}; % Converter Admittance
Yc = minreal(Yc(1:2,1:2)); 
[~,Ynet] = ObjYbusDss.GetDSS(ObjYbusDss); % Network Admittance
Ynet = Ynet(3:end,3:end); % Include infinity bus as ground bus

%% Define trasnformations
% Extract operating point
PF = ApparatusPowerFlow{2};
V = PF(3); phi = PF(4); v = V*exp(1i*phi);
R = [cos(phi) sin(phi); -sin(phi) cos(phi)];
P = PF(1); Q = PF(2);
i = (P-1i*Q)/v';
vd0 = real(v); vq0 = imag(v);
id0 = real(i); iq0 = imag(i);

% Define fiters
lp_w = 0.5*2*pi;
lp = lp_w/(s+lp_w);
hp_w = 0.5*2*pi;
hp = 1-hp_w/(s+hp_w);
% Virtual admittance filter
Lv = 0.1/Wbase; Rv = 0.01*1; z = abs(Rv+1i*Lv*Wbase);
Yv = [s*Lv+Rv -Wbase*Lv; Wbase*Lv s*Lv+Rv]/z; Yv = inv(Yv);
% Select trasnformation
switch 4
    case 1
        % Stardard admittance reference frame (Fig. 2)
        E = eye(2);
        C = zeros(2,2);
        F = eye(2);
    case 2
        % Stardard admittance with virtual admittance filter (Fig. 6)
        E = eye(2);
        C = zeros(2,2);
        F = inv(Yv);
    case 3
        % Power polar reference frame (Fig. 5)
        E = [vd0 vq0; vq0 -vd0];
        F = [-vq0 vd0; vd0 vq0];
        C = [id0 iq0;-iq0 id0]; 
    case 4
        % Mixed filtered transfomation (Fig. 7)
        E = [vd0 vq0; vq0 -vd0];
        F = [-vq0 vd0; vd0 vq0]*lp+inv(Yv)*inv(E)*hp;
        C = [id0 iq0;-iq0 id0]*lp;
end
%% Apply transformation
Yc = (E*Yc+C)*F;
Ynet = (E*Ynet-C)*F;
%% Compute Singular Values
sv_Yc = sigma(Yc,w);
sv_Yc_max = sv_Yc(1,:);
      
sv_Ynet = sigma(Ynet,w);
sv_Ynet_min = sv_Ynet(end,:);
%% Compute Phases - Converter
phase_Yc_max = [];
phase_Yc_min = [];

for i = 1:length(w)
    ww = w(i);
    A = evalfr(Yc,1i*ww);
    e = eig(A);
    % Verify A sectorial
    if check_sectoriality(A)
        phases = compute_phases(A);
    else
        phases = [-2*pi,2*pi];
    end
    phase_Yc_max = [phase_Yc_max phases(2)];
    phase_Yc_min = [phase_Yc_min phases(1)];
end

%% Compute Phases - Network
Znet = inv(Ynet);
phase_Znet_max = [];
phase_Znet_min = [];

for i = 1:length(w)
    ww = w(i);
    A = evalfr(Znet,1i*ww);
    % Verify A sectorial
    if check_sectoriality(A)
        phases = compute_phases(A);
    else
        phases = [-2*pi,2*pi];
    end

    phase_Znet_max = [phase_Znet_max phases(2)];
    phase_Znet_min = [phase_Znet_min phases(1)];
end
%% Compute Eigenvalues Phase Closed loop system
phase_CL_max = [];
phase_CL_min = [];
eig_cl = [];

for i = 1:length(w)
    ww = w(i);
    A = evalfr(Yc*Znet,1i*ww);
    e = eig(A);
    phase_CL_max = [phase_CL_max angle(e(1))];
    phase_CL_min = [phase_CL_min angle(e(2))];
end

%% Plots - Small gain
ax1 = subplot(3,1,1);
semilogx(f,20*log10(sv_Yc_max),'b',f,20*log10(sv_Ynet_min),'k');
grid on; 
title('Small-Gain Criteria');
ylabel('Gain (dB)')

%% Plots - Small phase
ax2 = subplot(3,1,2);
semilogx(f,phase_Yc_max,'b',f,pi-phase_Znet_max,'k'); hold on;
semilogx(f,phase_Yc_min,'b',f,pi-phase_Znet_min,'k');
semilogx(f,-pi-phase_Znet_max,'k');
semilogx(f,-pi-phase_Znet_min,'k');
x_fill = [f, fliplr(f)];
y_fill = [phase_Yc_max, fliplr(phase_Yc_min)];
fill(x_fill, y_fill, [0.7 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5)
y_fill = [pi-phase_Znet_max, fliplr(pi-phase_Znet_min)];
fill(x_fill, y_fill, [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3)
y_fill = [-pi-phase_Znet_max, fliplr(-pi-phase_Znet_min)];
fill(x_fill, y_fill, [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3)
legend('Converter', 'Grid');
title('Small-Phase Criteria');
ylabel('Phase (rad)')
grid on; ylim([-2*pi 2*pi]);

%% Plots - Comparison Small phase bounds and actual eigenvalues phase
ax3 = subplot(3,1,3); hold off;
semilogx(f,phase_Yc_max+phase_Znet_max,'g',f, ...
    pi*ones(size(phase_Znet_max)),'k'); hold on;
semilogx(f,1*phase_Yc_min+phase_Znet_min,'g',f, ...
    -pi*ones(size(phase_Znet_min)),'k');
semilogx(f,phase_CL_max,'b--',f,phase_CL_min,'b--');
legend('Phase sum', '','','', 'Actual Phase');
ylabel('Phase (rad)')
title('Check Phase Bounds');
grid on
linkaxes([ax1,ax2,ax3],'x');

xlim([f(1) f(end)])
ylim([-2*pi 2*pi]);
xlabel('f (Hz)');