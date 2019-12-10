% plots the mean phase offset for a range of coupling strengths
% for each noise model. 
%
% load data
load('run3_phase_offsets_withstd', 'mdp_su_means', 'mdp_cd_means', 'mdp_vc_means', 'mdp_ode_means', 'mdp_mc_means', ...
    'mdp_su_err', 'mdp_cd_err', 'mdp_vc_err', 'mdp_ode_err', 'mdp_mc_err');

ntrials = 25;
n_neurons = 2;
noise_intensity = 1;
coupling_strength = 0:0.1:1;
t = [0:0.01:300];
noise = {'MarkovChain', 'Subunit', 'VClamp', 'FoxLuSystemSize', 'ODE'};
noise_legend = {'Markov Chain', 'Subunit', 'V. Clamp', 'Syst. Size', 'ODE'};

figure(1);
hold on;
plot(coupling_strength, mdp_ode_means, 'k');
plot(coupling_strength, mdp_su_means, 'b');
plot(coupling_strength, mdp_vc_means, 'color', [0 0.5 0.5]);
plot(coupling_strength, mdp_cd_means, 'r');
plot(coupling_strength, mdp_mc_means, 'color', [0.5 0.1 0.8]);
errorbar(coupling_strength, mdp_ode_means,mdp_ode_err/5, 'k');
errorbar(coupling_strength, mdp_su_means, mdp_su_err/5, 'b');
errorbar(coupling_strength, mdp_vc_means,mdp_vc_err/5, 'color', [0 0.5 0.5]);
errorbar(coupling_strength, mdp_cd_means,mdp_cd_err/5, 'r');
errorbar(coupling_strength, mdp_mc_means, mdp_mc_err/5, 'color', [0.5 0.1 0.8]);
xlabel('Coupling Strength, \kappa', 'fontsize', 15);
ylabel('Mean Phase Offset, <\Delta\Phi>', 'fontsize', 15);
title('Mean Phase Offsets Between Coupled Neurons', 'fontsize', 15, 'fontweight', 'normal');
ylim([0, pi/2]);
yticks([0 pi/4 pi/2]);
yticklabels({'0', '\pi/4', '\pi/2'});
legend(noise_legend);
legend boxoff;
    