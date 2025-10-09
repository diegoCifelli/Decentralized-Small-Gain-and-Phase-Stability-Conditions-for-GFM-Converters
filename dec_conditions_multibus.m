%% Decentralized Stability Conditions for a multi-converter power system
% Script to use for the IEEE14 Bus System
%% Setup
% Frequency vector
w = logspace(log10(1e-3),log10(1e4),1000)*2*pi;
f = w/2/pi;
% s operator
s = tf('s');

%% Extract converters' model and network model
appT = cell2mat(ApparatusType); % Apparatus types
appBus = find(appT<=90); % Buses with apparatus
appEmpty = find(appT>90); % Buses without apparatus
inf_bus = 1;

% Extract converters' model
Yconv = {};
for a=appBus
    c_t = GmDssCell{a};
    c_t = c_t([1 2],[1 2]);
    c_t = minreal(c_t);
    Yconv{a} = c_t;
end
% Extract converters' model
[~,Ynet] = ObjYbusDss.GetDSS(ObjYbusDss); % 
% Include inf bus as short-circuit
Yconv{inf_bus} = tf(eye(2)*1e14);
Ynet_a = include_converter_into_net_y(Ynet,Yconv{inf_bus},inf_bus);
% Kron Reduction
Ynet_red = kron_reduction_dq(Ynet_a,sort([appEmpty(1:end),inf_bus]));
Ynet_red = sminreal(Ynet_red);
appBus =setdiff(appBus,inf_bus);
%% Define trasnformations
C = {};
E = {};
F = {};
% Define fiters
lp_w = 0.5*2*pi;
lp = lp_w/(s+lp_w);
hp_w = 0.5*2*pi;
hp = 1-hp_w/(s+hp_w);
% Virtual admittance filter
Lv = 0.1/Wbase; Rv = 0.01; z = abs(Rv+1i*Lv*Wbase);
Yv = [s*Lv+Rv -Wbase*Lv; Wbase*Lv s*Lv+Rv]/z; Yv = inv(Yv);
for a=appBus
    % Get setpoint
    PF = ApparatusPowerFlow{a};
    V = PF(3); phi = PF(4); v = V*exp(1i*phi);
    P = PF(1); Q = PF(2);
    i = (P-1i*Q)/v';
    vd0 = real(v); vq0 = imag(v);
    id0 = real(i); iq0 = imag(i);
    e = [vd0 vq0; vq0 -vd0];
    E{a} = e;
    c = [id0 iq0;-iq0 id0]*lp;
    C{a} = c;
    f = [-vq0 vd0/V; vd0 vq0/V]*lp+inv(Yv)*inv(e)*hp;
    F{a} = f;   
end

%% Apply transformation to Yconv
if 1
    k = 1;
    for a=appBus
        Yconv{a} = (E{a}*Yconv{a}+C{a})*F{a};
        k = k+1;
    end
end
%% Apply transformation to Ynet
EE =[]; CC = []; FF = [];
if 1
    k = 1;
    for a=appBus
        EE = blkdiag(EE,E{a});
        FF = blkdiag(FF,F{a});
        CC = blkdiag(CC,C{a});
        k = k+1;
    end
    Ynet_red_t = (EE*Ynet_red-CC)*FF;
else
    Ynet_red_t = Ynet_red;
end

%% Compute singular values
% Converters
sv_Yc_max_v = {};
for a=appBus
    c = Yconv{a};
    sv_Yc = sigma(c,w);
    sv_Yc_max = sv_Yc(1,:);
    sv_Yc_max_v{a} = sv_Yc_max;
end    

% Network
[sv_net] = sigma(Ynet_red_t,w);
sv_Ynet_min = sv_net(end,:);
%% Phase converters
phase_Yc_max = {};
phase_Yc_min = {};
for a=appBus
    c = Yconv{a};
    p_max = [];p_min = [];
    for i = 1:length(w)
        ww = w(i);
        A = evalfr(c,1i*ww);
        % Verify A sectorial
        if check_sectoriality(A)
            phases = compute_phases(A);
        else
            phases = [-2*pi,2*pi];
        end
        p_max = [p_max phases(2)];
        p_min = [p_min phases(1)];
     end
     phase_Yc_max{a} = p_max;
     phase_Yc_min{a} = p_min;
end   

%% Phase network
Znet = inv(Ynet_red_t);
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

%% Plots - Small gain
colors = linspecer(length(appBus)); % Colors for plot
ax1 = subplot(2,1,1);
k = 1;
for a=appBus
    semilogx(w/2/pi, 20*log10(sv_Yc_max_v{a}),'Color',colors(k,:)); hold on;
    k = k+1;
end
semilogx(w/2/pi, 20*log10(sv_Ynet_min),'-k'); hold on;
grid on; 
title('Small-Gain Criteria');
ylabel('Gain (dB)')
%% Plots - Small phase
ax2 = subplot(2,1,2);
k = 1;
f = w/2/pi;
p_max_tot=phase_Yc_max{appBus(1)};
p_min_tot=phase_Yc_min{appBus(1)};
% Convert phase
for a=appBus
     p_max = phase_Yc_max{a};
     p_min = phase_Yc_min{a};
     semilogx(f,p_max*180/pi,'Color',colors(k,:)); hold on;
     semilogx(f,p_min*180/pi,'Color',colors(k,:));
%      p_max_tot = max(p_max_tot,p_max);
%      p_min_tot = min(p_min_tot,p_min);
     x_fill = [f, fliplr(f)];
     y_fill = [p_max*180/pi, fliplr(p_min*180/pi)]; hold on
     fill(x_fill, y_fill, colors(k,:), 'EdgeColor', 'none', 'FaceAlpha', 0.3)
     k = k+1;
end  
% Network phase
semilogx(f,180-phase_Znet_max*180/pi,'k',f,-180-phase_Znet_max*180/pi,'k'); 
semilogx(f,180-phase_Znet_min*180/pi,'k',f,-180-phase_Znet_min*180/pi,'k');
y_fill = [180-phase_Znet_max*180/pi, fliplr(180-phase_Znet_min*180/pi)];
fill(x_fill, y_fill, [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3)
y_fill = [-180-phase_Znet_max*180/pi, fliplr(-180-phase_Znet_min*180/pi)];
fill(x_fill, y_fill, [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3)
lng = strings(3*length(appBus)+2,1);
lng(1:3:end-2) = 'J_{C,' + string(appBus) +'}';
lng(end) = 'J_{net}';
legend(lng);
title('Phase Test')
grid on
xlim([f(1) f(end)])
xlabel('F [Hz]')

linkaxes([ax1,ax2],'x');
xlim([f(1) f(end)])
xlabel('F [Hz]')

%% Functions
function Y_net_w_c = include_converter_into_net_y(net_y,c,bus)
    net_y_dim = size(net_y,1);
    % Compute converter indexes to insert in network admittance
    a = bus*2-3+2; b = sort([a, a+1]);
    idx = b;
    % Build full matrix converter to insert in network admittance
    Yc = ss(zeros(net_y_dim,net_y_dim),zeros(net_y_dim,net_y_dim),zeros(net_y_dim,net_y_dim),[]); % Empty ss
    Yc(idx,idx) = c; % Add converters admittance
    % Compute network admittance with converters admittance
    Y_net_w_c = (net_y+Yc);   
end