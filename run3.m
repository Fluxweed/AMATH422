% run 3

% mean phase offset for a range of coupling strengths for each  
% noise model 

n_neurons = 2;
ntrials = 25;
noise_intensity = 1;
coupling_strength = 0:0.1:1;
t = [0:0.01:300];
noise = {'Subunit', 'FoxLuSystemSize', 'VClamp', 'MarkovChain', 'ODE'};

mdp_su = zeros(1, length(ntrials));
mdp_cd = zeros(1, length(ntrials));
mdp_vc = zeros(1, length(ntrials));
mdp_ode = zeros(1, length(ntrials)); 
mdp_mc = zeros(1, length(ntrials));

mdp_su_means = zeros(1, length(coupling_strength));
mdp_cd_means = zeros(1, length(coupling_strength));
mdp_vc_means = zeros(1, length(coupling_strength));
mdp_ode_means = zeros(1, length(coupling_strength)); 
mdp_mc_means = zeros(1, length(coupling_strength));

for i = 1:length(coupling_strength)
    parfor j = 1:ntrials

        kappa = coupling_strength(i);   
        [delta_phi_su, lengths_su] = multiple_phase_offsets(t, n_neurons, kappa, noise_intensity, noise{1});
        [delta_phi_cd, lengths_cd] = multiple_phase_offsets(t, n_neurons, kappa, noise_intensity, noise{2});
        [delta_phi_vc, lengths_vc] = multiple_phase_offsets(t, n_neurons, kappa, noise_intensity, noise{3});
        [delta_phi_ode, lengths_ode] = multiple_phase_offsets(t, n_neurons, kappa, noise_intensity, noise{4});
        [delta_phi_mc, lengths_mc] = multiple_phase_offsets(t, n_neurons, kappa, noise_intensity, noise{5});
        
        mdp_su(j) = mean(dist_from_asynchrony(delta_phi_su(1:lengths_su, :)));
        mdp_cd(j) = mean(dist_from_asynchrony(delta_phi_cd(1:lengths_cd, :)));
        mdp_vc(j) = mean(dist_from_asynchrony(delta_phi_vc(1:lengths_vc, :)));
        mdp_ode(j) = mean(dist_from_asynchrony(delta_phi_ode(1:lengths_ode, :)));
        mdp_mc(j) = mean(dist_from_asynchrony(delta_phi_mc(1:lengths_mc, :)));

    end
    mdp_su_means(i) = mean(mdp_su);
    mdp_cd_means(i) = mean(mdp_cd);
    mdp_vc_means(i) = mean(mdp_vc);
    mdp_ode_means(i) = mean(mdp_ode);
    mdp_mc_means(i) = mean(mdp_mc);
end

save('run3_phase_offsets', 'mdp_su_means', 'mdp_cd_means', 'mdp_vc_means', 'mdp_ode_means', 'mdp_mc_means');
    