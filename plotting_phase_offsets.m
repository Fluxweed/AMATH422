clear all; close all; clc;

%% plot delta_phi as a function of coupling strength
n_neurons = 10;
noise_intensity = 1;
coupling_strength = 1;
t = [0:0.01:100];
noise = 'Subunit';

[delta_phi, lengths] = multiple_phase_offsets(t, n_neurons, coupling_strength, noise_intensity, noise);

%% probability density function estimation of phase offsets
legend_names = {'\Delta\Phi_{1,2}', '\Delta\Phi_{2,3}', '\Delta\Phi_{3,4}',
                '\Delta\Phi_{4,5}', '\Delta\Phi_{5,6}', '\Delta\Phi_{6,7}',
                '\Delta\Phi_{7,8}', '\Delta\Phi_{8,9}', '\Delta\Phi_{9,10}'};

figure(1);
lineStyles = linspecer(n_neurons);
hold on;
for i = 1:n_neurons-1
    dp = dist_from_asynchrony(delta_phi(:, i));
    [f_,kd_] = ksdensity(dp, 'bandwidth', 0.1, 'Function', 'cdf'); 
    plot(kd_, f_, 'color', lineStyles(i,:));
end
xlabel('Phase Offset, \Delta\Phi');
ylabel('Probability Density');
title(['Phase Offsets Between Coupled Neurons (\kappa = ', num2str(coupling_strength), ')']);
xlim([0,1]);
lgd = legend(legend_names());
% lgd.Position = [0.7363 0.6782 0.1000 0.2226];
legend boxoff; 

figure(2);
lineStyles = linspecer(n_neurons);
hold on;
for i = 1:n_neurons-1
    [f,kd] = ksdensity(delta_phi(1:lengths(i), i), 'bandwidth', 0.1, 'Function', 'cdf'); 
    plot(kd, f, 'color', lineStyles(i,:));
end
xlabel('Phase Offset, \Delta\Phi');
ylabel('Probability Density');
title(['Phase Offsets Between Coupled Neurons (\kappa = ', num2str(coupling_strength), ')']);
xlim([0,2*pi]);
lgd = legend(legend_names);
% lgd.Position = [0.7363 0.6782 0.1000 0.2226];
legend boxoff;

%% histograms of phase offsets         
figure(3); 
hold on; 
legend_names = {'\Delta\Phi_{1,2}', '\Delta\Phi_{2,3}', '\Delta\Phi_{3,4}', '\Delta\Phi_{4,5}'};
histogram(delta_phi_12, 'facealpha', 0.5);
histogram(delta_phi_23, 'facealpha', 0.5);
histogram(delta_phi_34, 'facealpha', 0.5);
histogram(delta_phi_45, 'facealpha', 0.5);
ylabel('Frequency');
xlabel('Phase Offset');
title('Histgram of Phase Offsets between Coupled Neurons')
legend(legend_names);
legend boxoff;
