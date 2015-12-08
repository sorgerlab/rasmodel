t0 = 0;
tf = 4e4;
tspan = [t0 tf];

m = run_timecourse();
m.parameters.Vem_0 = 0;
% Parameters to control negative feedback
m.parameters.k_spe = 1e-4;  % k_spe = 1e-3 and k_dspe = 50 work well too
m.parameters.k_dspe = 10;
m.parameters.k_pp2e = 10;
m.parameters.k_mer = 0.1;
m.parameters.k_mee = 10;
m.parameters.kf5 = 0.25; % Important parameter to control extent of adaptive response
m.parameters.KRAS_0 = 1e5;
m.parameters.BRAF_0 = 1e4;

% Simulate mutants
% ================
% NRAS mutation
% m.parameters.kf5 = 0.05;

% BRAF overexpression
% m.parameters.BRAF_0 = 1e6;


options = odeset('Events', @add_vemurafenib, 'RelTol', 1e-4, 'AbsTol', 1e-6);

tout = t0;
yout = m.get_initial_values;
teout  = [];
yeout = [];
ieout = [];

initial_values = m.get_initial_values();

while t0 < tf
    
    [t, y, te, ye, ie] = ode15s(@m.odes, [t0 tf], initial_values, options);
    
    nt = length(t);
    
    tout = [tout; t(2:nt)];
    yout = [yout; y(2:nt, :)];
    teout = [teout; te];
    yeout = [yeout; ye];
    ieout = [ieout; ie];
    
    initial_values = y(nt, :);
    
    ie
    if isscalar(ie) == 0
        ie = 0;
    end
    
    if ie == 1
        initial_values(8) = 2e4;
        initial_values(13) = 0;%2e5;
    end
    t0 = t(nt);
    
    if t0 >= tf
        break;
    end
end
    
y_obs = m.get_observables(yout);

plot(tout, y_obs.ERK_P/m.parameters.ERK_0, 'LineWidth', 4)
axis([1.5e4 3e4 0 1])
xlabel('time', 'Fontsize', 20)
ylabel('ERK~P', 'Fontsize', 20)

% hold on
% plot(erk_out{2, 1}, erk_out{2, 2}, 'k')
% plot(erk_out{3, 1}, erk_out{3, 2}, 'g')
% plot(erk_out{4, 1}, erk_out{4, 2}, 'b')
% plot(erk_out{1, 1}, erk_out{1, 2}, 'r')
% axis([4e4 6.5e4 0 1])
% xlabel('time (a.u)', 'Fontsize', 20)
% ylabel('ERK~P', 'Fontsize', 20)
% h_legend = legend('KRAS mut', 'BRAF(V600E) ovex', 'Trimatinib', 'BRAF(V600E)') ;
% set(h_legend,'FontSize',14)
% title('Addition of Vemurafenib (2e5 a.u) at t = 5e4 a.u', 'Fontsize', 15)
%     