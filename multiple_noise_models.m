% plots the probability density functions and cumulative distribution
% functions of the phase offsets for the different noise models

ntrials = 10;
n_neurons = 2;
noise_intensity = 1;
coupling_strength = 1;
t = [0:0.01:500];
noise = {'Subunit', 'FoxLuSystemSize', 'VClamp', 'MarkovChain'};

delta_phi_su = [];
delta_phi_cd = [];
delta_phi_vc = [];
delta_phi_mc = [];
for i = 1:ntrials
    disp(i);
    [offset_su, lengths_su] = multiple_phase_offsets(t, n_neurons, coupling_strength, noise_intensity, noise{1});
    [offset_cd, lengths_cd] = multiple_phase_offsets(t, n_neurons, coupling_strength, noise_intensity, noise{2});
    [offset_vc, lengths_vc] = multiple_phase_offsets(t, n_neurons, coupling_strength, noise_intensity, noise{3});
    [offset_mc, lengths_mc] = multiple_phase_offsets(t, n_neurons, coupling_strength, noise_intensity, noise{4});
    
    delta_phi_su = [delta_phi_su; offset_su(1:lengths_su, :)];
    delta_phi_cd = [delta_phi_cd; offset_cd(1:lengths_cd, :)];
    delta_phi_vc = [delta_phi_vc; offset_vc(1:lengths_vc, :)];
    delta_phi_mc = [delta_phi_mc; offset_mc(1:lengths_mc, :)];
end

figure(1);
lineStyles = linspecer(n_neurons);
hold on;

for i = 1:n_neurons-1
    dp_su = dist_from_asynchrony(delta_phi_su(:, i));
    [f_su,kd_su] = ksdensity(dp_su, 'bandwidth', 0.1); 
    plot(kd_su, f_su, 'color', 'red');
    dp_cd = dist_from_asynchrony(delta_phi_cd(:, i));
    [f_cd,kd_cd] = ksdensity(dp_cd, 'bandwidth', 0.1); 
    plot(kd_cd, f_cd, 'color', 'blue');
    dp_vc = dist_from_asynchrony(delta_phi_vc(:, i));
    [f_vc,kd_vc] = ksdensity(dp_vc, 'bandwidth', 0.1); 
    plot(kd_vc, f_vc, 'color', 'green');
    dp_mc = dist_from_asynchrony(delta_phi_mc(:, i));
    [f_mc,kd_mc] = ksdensity(dp_mc, 'bandwidth', 0.1); 
    plot(kd_mc, f_mc, 'color', 'black');
end
xlabel('Phase Offset, \Delta\Phi');
ylabel('Probability Density');
title(['Phase Offsets Between Coupled Neurons (\kappa = ', num2str(coupling_strength), ')']);
xlim([0,1]);
lgd = legend(noise);
% lgd.Position = [0.7363 0.6782 0.1000 0.2226];
legend boxoff; 

figure(2);
lineStyles = linspecer(n_neurons);
hold on;
for i = 1:n_neurons-1
    [f_su,kd_su] = ksdensity(delta_phi_su(1:lengths_su(i), i), 'bandwidth', 0.1); 
    plot(kd_su, f_su, 'color', 'red');
    [f_cd,kd_cd] = ksdensity(delta_phi_cd(1:lengths_cd(i), i), 'bandwidth', 0.1); 
    plot(kd_cd, f_cd, 'color', 'blue');
    [f_vc,kd_vc] = ksdensity(delta_phi_vc(1:lengths_vc(i), i), 'bandwidth', 0.1); 
    plot(kd_vc, f_vc, 'color', 'green');
    [f_mc,kd_mc] = ksdensity(delta_phi_mc(1:lengths_mc(i), i), 'bandwidth', 0.1); 
    plot(kd_mc, f_mc, 'color', 'black');
end
xlabel('Phase Offset, \Delta\Phi');
ylabel('Probability Density');
title(['Phase Offsets Between Coupled Neurons (\kappa = ', num2str(coupling_strength), ')']);
xlim([0,2*pi]);
lgd = legend(noise);
% lgd.Position = [0.7363 0.6782 0.1000 0.2226];
legend boxoff;