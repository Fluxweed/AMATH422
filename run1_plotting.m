% % plots the probability density functions and cumulative distribution
% functions of the phase offsets for the different noise models. the
% phase offsets are redefined to represent the how far the phase
% offets are from a state of asynchrony.

% coupling strength 0.1
load('run1_coupling01', 'delta_phi_su', 'delta_phi_cd', 'delta_phi_vc', 'delta_phi_mc');
ntrials = 100;
n_neurons = 2;
noise_intensity = 1;
coupling_strength = 0.1;
t = [0:0.01:300];
noise = {'Subunit', 'FoxLuSystemSize', 'VClamp', 'MarkovChain'};
noise_legend = {'Markov Chain', 'Subunit', 'V. Clamp', 'Syst. Size'};

% probability density function
figure(1);
hold on;
for i = 1:n_neurons-1
    dp_mc = dist_from_asynchrony(delta_phi_mc(:, i));
    [f_mc,kd_mc] = ksdensity(dp_mc, 'bandwidth', 0.1); 
    plot(kd_mc, f_mc, 'k');
    dp_su = dist_from_asynchrony(delta_phi_su(:, i));
    [f_su,kd_su] = ksdensity(dp_su, 'bandwidth', 0.1); 
    plot(kd_su, f_su, 'b');
    dp_vc = dist_from_asynchrony(delta_phi_vc(:, i));
    [f_vc,kd_vc] = ksdensity(dp_vc, 'bandwidth', 0.1); 
    plot(kd_vc, f_vc, 'color', [0 0.5 0.5]);
    dp_cd = dist_from_asynchrony(delta_phi_cd(:, i));
    [f_cd,kd_cd] = ksdensity(dp_cd, 'bandwidth', 0.1); 
    plot(kd_cd, f_cd, 'r');
end
xlabel('Phase Offset, \Delta\Phi', 'fontsize', 15);
ylabel('Probability Density', 'fontsize', 15);
title(['Phase Offsets Between Coupled Neurons (\kappa = ', num2str(coupling_strength), ')'], 'fontsize', 15, 'fontweight', 'normal');
xlim([0,pi]);
xticks([0 pi/4 pi/2 3*pi/4 pi]);
xticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'});
lgd = legend(noise_legend);
legend boxoff; 

% probability density function of phase offsets from 0 to 2pi
figure(2);
hold on;
for i = 1:n_neurons-1
    [f_mc,kd_mc] = ksdensity(delta_phi_mc(:, i), 'bandwidth', 0.1); 
    plot(kd_mc, f_mc, 'k');
    [f_su,kd_su] = ksdensity(delta_phi_su(:, i), 'bandwidth', 0.1); 
    plot(kd_su, f_su, 'b');
    [f_vc,kd_vc] = ksdensity(delta_phi_vc(:, i), 'bandwidth', 0.1); 
    plot(kd_vc, f_vc, 'color', [0 0.5 0.5]);
    [f_cd,kd_cd] = ksdensity(delta_phi_cd(:, i), 'bandwidth', 0.1); 
    plot(kd_cd, f_cd, 'r');
end
xlabel('Phase Offset, \Delta\Phi', 'fontsize', 15);
ylabel('Probability Density', 'fontsize', 15);
title(['Phase Offsets Between Coupled Neurons (\kappa = ', num2str(coupling_strength), ')'], 'fontsize', 15, 'fontweight', 'normal');
xlim([0,2*pi]);
lgd = legend(noise_legend);
legend boxoff;

% cumulative distribution function
figure(3);
hold on;
for i = 1:n_neurons-1
    dp_mc = dist_from_asynchrony(delta_phi_mc(:, i));
    [f_mc,kd_mc] = ksdensity(dp_mc, 'bandwidth', 0.1, 'function', 'cdf'); 
    plot(kd_mc, f_mc, 'k');
    dp_su = dist_from_asynchrony(delta_phi_su(:, i));
    [f_su,kd_su] = ksdensity(dp_su, 'bandwidth', 0.1, 'function', 'cdf'); 
    plot(kd_su, f_su, 'b');
    dp_vc = dist_from_asynchrony(delta_phi_vc(:, i));
    [f_vc,kd_vc] = ksdensity(dp_vc, 'bandwidth', 0.1, 'function', 'cdf'); 
    plot(kd_vc, f_vc, 'color', [0 0.5 0.5]);
    dp_cd = dist_from_asynchrony(delta_phi_cd(:, i));
    [f_cd,kd_cd] = ksdensity(dp_cd, 'bandwidth', 0.1, 'function', 'cdf'); 
    plot(kd_cd, f_cd, 'r');
end
xlabel('Phase Offset, \Delta\Phi', 'fontsize', 15);
ylabel('Cumulative Probability', 'fontsize', 15);
title(['Phase Offsets Between Coupled Neurons (\kappa = ', num2str(coupling_strength), ')'], 'fontsize', 15, 'fontweight', 'normal');
xlim([0,pi]);
xticks([0 pi/4 pi/2 3*pi/4 pi]);
xticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'});

% coupling strength 0.3
load('run1_coupling03', 'delta_phi_su', 'delta_phi_cd', 'delta_phi_vc', 'delta_phi_mc');
ntrials = 100;
n_neurons = 2;
noise_intensity = 1;
coupling_strength = 0.3;
t = [0:0.01:300];
noise = {'Subunit', 'FoxLuSystemSize', 'VClamp', 'MarkovChain'};

% probability density function
figure(1);
hold on;
for i = 1:n_neurons-1
    dp_mc = dist_from_asynchrony(delta_phi_mc(:, i));
    [f_mc,kd_mc] = ksdensity(dp_mc, 'bandwidth', 0.1); 
    plot(kd_mc, f_mc, 'k');
    dp_su = dist_from_asynchrony(delta_phi_su(:, i));
    [f_su,kd_su] = ksdensity(dp_su, 'bandwidth', 0.1); 
    plot(kd_su, f_su, 'b');
    dp_vc = dist_from_asynchrony(delta_phi_vc(:, i));
    [f_vc,kd_vc] = ksdensity(dp_vc, 'bandwidth', 0.1); 
    plot(kd_vc, f_vc, 'color', [0 0.5 0.5]);
    dp_cd = dist_from_asynchrony(delta_phi_cd(:, i));
    [f_cd,kd_cd] = ksdensity(dp_cd, 'bandwidth', 0.1); 
    plot(kd_cd, f_cd, 'r');
end
xlabel('Phase Offset, \Delta\Phi', 'fontsize', 15);
ylabel('Probability Density', 'fontsize', 15);
title(['Phase Offsets Between Coupled Neurons (\kappa = ', num2str(coupling_strength), ')'], 'fontsize', 15, 'fontweight', 'normal');
xlim([0,pi]);
xticks([0 pi/4 pi/2 3*pi/4 pi]);
xticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'});
lgd = legend(noise_legend);
legend boxoff; 

% probability density function of phase offsets from 0 to 2pi
figure(2);
hold on;
for i = 1:n_neurons-1
    [f_mc,kd_mc] = ksdensity(delta_phi_mc(1:lengths_mc(i), i), 'bandwidth', 0.1); 
    plot(kd_mc, f_mc, 'k');
    [f_su,kd_su] = ksdensity(delta_phi_su(1:lengths_su(i), i), 'bandwidth', 0.1); 
    plot(kd_su, f_su, 'b');
    [f_vc,kd_vc] = ksdensity(delta_phi_vc(1:lengths_vc(i), i), 'bandwidth', 0.1); 
    plot(kd_vc, f_vc, 'color', [0 0.5 0.5]);
    [f_cd,kd_cd] = ksdensity(delta_phi_cd(1:lengths_cd(i), i), 'bandwidth', 0.1); 
    plot(kd_cd, f_cd, 'r');
end
xlabel('Phase Offset, \Delta\Phi', 'fontsize', 15);
ylabel('Probability Density', 'fontsize', 15);
title(['Phase Offsets Between Coupled Neurons (\kappa = ', num2str(coupling_strength), ')'], 'fontsize', 15, 'fontweight', 'normal');
xlim([0,2*pi]);
lgd = legend(noise_legend);
legend boxoff;

% cumulative density function
figure(3);
hold on;
for i = 1:n_neurons-1
    dp_mc = dist_from_asynchrony(delta_phi_mc(:, i));
    [f_mc,kd_mc] = ksdensity(dp_mc, 'bandwidth', 0.1, 'function', 'cdf'); 
    plot(kd_mc, f_mc, 'k');
    dp_su = dist_from_asynchrony(delta_phi_su(:, i));
    [f_su,kd_su] = ksdensity(dp_su, 'bandwidth', 0.1, 'function', 'cdf'); 
    plot(kd_su, f_su, 'b');
    dp_vc = dist_from_asynchrony(delta_phi_vc(:, i));
    [f_vc,kd_vc] = ksdensity(dp_vc, 'bandwidth', 0.1, 'function', 'cdf'); 
    plot(kd_vc, f_vc, 'color', [0 0.5 0.5]);
    dp_cd = dist_from_asynchrony(delta_phi_cd(:, i));
    [f_cd,kd_cd] = ksdensity(dp_cd, 'bandwidth', 0.1, 'function', 'cdf'); 
    plot(kd_cd, f_cd, 'r');
end
xlabel('Phase Offset, \Delta\Phi', 'fontsize', 15);
ylabel('Cumulative Probability', 'fontsize', 15);
title(['Phase Offsets Between Coupled Neurons (\kappa = ', num2str(coupling_strength), ')'], 'fontsize', 15, 'fontweight', 'normal');
xlim([0,pi]);
xticks([0 pi/4 pi/2 3*pi/4 pi]);
xticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'});

% coupling strength 0.9
load('run1_coupling09', 'delta_phi_su', 'delta_phi_cd', 'delta_phi_vc', 'delta_phi_mc');
ntrials = 100;
n_neurons = 2;
noise_intensity = 1;
coupling_strength = 0.9;
t = [0:0.01:300];
noise = {'Subunit', 'FoxLuSystemSize', 'VClamp', 'MarkovChain'};

% probability density function
figure(1);
hold on;
for i = 1:n_neurons-1
    dp_mc = dist_from_asynchrony(delta_phi_mc(:, i));
    [f_mc,kd_mc] = ksdensity(dp_mc, 'bandwidth', 0.1); 
    plot(kd_mc, f_mc, 'k');
    dp_su = dist_from_asynchrony(delta_phi_su(:, i));
    [f_su,kd_su] = ksdensity(dp_su, 'bandwidth', 0.1); 
    plot(kd_su, f_su, 'b');
    dp_vc = dist_from_asynchrony(delta_phi_vc(:, i));
    [f_vc,kd_vc] = ksdensity(dp_vc, 'bandwidth', 0.1); 
    plot(kd_vc, f_vc, 'color', [0 0.5 0.5]);
    dp_cd = dist_from_asynchrony(delta_phi_cd(:, i));
    [f_cd,kd_cd] = ksdensity(dp_cd, 'bandwidth', 0.1); 
    plot(kd_cd, f_cd, 'r');
end
xlabel('Phase Offset, \Delta\Phi', 'fontsize', 15);
ylabel('Probability Density', 'fontsize', 15);
title(['Phase Offsets Between Coupled Neurons (\kappa = ', num2str(coupling_strength), ')'], 'fontsize', 15, 'fontweight', 'normal');
xlim([0,pi]);
xticks([0 pi/4 pi/2 3*pi/4 pi]);
xticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'});
lgd = legend(noise_legend);
legend boxoff; 

% pdf of phase offsets from 0 to 2pi
figure(2);
lineStyles = linspecer(n_neurons);
hold on;
for i = 1:n_neurons-1
    [f_mc,kd_mc] = ksdensity(delta_phi_mc(1:lengths_mc(i), i), 'bandwidth', 0.1); 
    plot(kd_mc, f_mc, 'k');
    [f_su,kd_su] = ksdensity(delta_phi_su(1:lengths_su(i), i), 'bandwidth', 0.1); 
    plot(kd_su, f_su, 'b')
    [f_vc,kd_vc] = ksdensity(delta_phi_vc(1:lengths_vc(i), i), 'bandwidth', 0.1); 
    plot(kd_vc, f_vc, 'color', [0 0.5 0.5]);
    [f_cd,kd_cd] = ksdensity(delta_phi_cd(1:lengths_cd(i), i), 'bandwidth', 0.1); 
    plot(kd_cd, f_cd, 'r');
end
xlabel('Phase Offset, \Delta\Phi', 'fontsize', 15);
ylabel('Probability Density', 'fontsize', 15);
title(['Phase Offsets Between Coupled Neurons (\kappa = ', num2str(coupling_strength), ')'], 'fontsize', 15, 'fontweight', 'normal');
xlim([0,2*pi]);
lgd = legend(noise_legend);
legend boxoff;

% cumulative distribution function
figure(3);
hold on;
for i = 1:n_neurons-1
    dp_mc = dist_from_asynchrony(delta_phi_mc(:, i));
    [f_mc,kd_mc] = ksdensity(dp_mc, 'bandwidth', 0.1, 'function', 'cdf'); 
    plot(kd_mc, f_mc, 'k');
    dp_su = dist_from_asynchrony(delta_phi_su(:, i));
    [f_su,kd_su] = ksdensity(dp_su, 'bandwidth', 0.1, 'function', 'cdf'); 
    plot(kd_su, f_su, 'b');
    dp_vc = dist_from_asynchrony(delta_phi_vc(:, i));
    [f_vc,kd_vc] = ksdensity(dp_vc, 'bandwidth', 0.1, 'function', 'cdf'); 
    plot(kd_vc, f_vc, 'color', [0 0.5 0.5]);
    dp_cd = dist_from_asynchrony(delta_phi_cd(:, i));
    [f_cd,kd_cd] = ksdensity(dp_cd, 'bandwidth', 0.1, 'function', 'cdf'); 
    plot(kd_cd, f_cd, 'r');
end
xlabel('Phase Offset, \Delta\Phi', 'fontsize', 15);
ylabel('Cumulative Probability', 'fontsize', 15);
title(['Phase Offsets Between Coupled Neurons (\kappa = ', num2str(coupling_strength), ')'], 'fontsize', 15, 'fontweight', 'normal');
xlim([0,pi]);
xticks([0 pi/4 pi/2 3*pi/4 pi]);
xticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'});