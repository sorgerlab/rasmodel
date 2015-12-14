t0 = 0;
tf = 1e5;
tspan = [t0 tf];

m = test_vr();
% m = vemurafenib_resistance();
m.parameters.Vem_0 = 0;
m.parameters.kr_rg_bind_1 = 0.25;
m.parameters.EGF_0 = 5e3;
% m.parameters.kf_ee_transphos_1 = 1;

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
        initial_values(2) = 2e5;
        %initial_values(8) = 2e5;
    end
    t0 = t(nt);
    
    if t0 >= tf
        break;
    end
end
    
y_obs = m.get_observables(yout);

plot(tout, y_obs.ERK_P/m.parameters.MAPK1_0, 'LineWidth', 2)
% axis([1.5e4 3e4 0 1])
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