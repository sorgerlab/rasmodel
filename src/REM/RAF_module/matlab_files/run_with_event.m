t0 = 0;
tf = 1e5;
tspan = [t0 tf];

m = run_timecourse();
m.parameters.Vem_0 = 0;
% Parameters to control negative feedback
m.parameters.k_spe = 1e-4;  % k_spe = 1e-3 and k_dspe = 50 work well too
m.parameters.k_dspe = 10;
m.parameters.kf5 = 0.5; % Important parameter to control extent of adaptive response 

% Simulate mutants
% ================
% KRAS mutation
% m.parameters.kf5 = 0.05;

% BRAF overexpression
m.parameters.BRAF_0 = 1e5;


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
        initial_values(8) = 2e5;
    end
    t0 = t(nt);
    
    if t0 >= tf
        break;
    end
end
    
y_obs = m.get_observables(yout);

plot(tout, y_obs.ERK_P/m.parameters.ERK_0)
axis([4e4 8e4 0 1])
xlabel('time (a.u)', 'Fontsize', 20)
ylabel('ERK~P', 'Fontsize', 20)
title('Addition of Vemurafenib (2e5 a.u) at t = 5e4 a.u', 'Fontsize', 15)
    