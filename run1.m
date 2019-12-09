% plots the probability density functions and cumulative distribution
% functions of the phase offsets for the different noise models

% coupling strength 0.1
ntrials = 100;
n_neurons = 2;
noise_intensity = 1;
coupling_strength = 0.1;
t = [0:0.01:300];
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

save('run1_coupling01', 'delta_phi_su', 'delta_phi_cd', 'delta_phi_vc', 'delta_phi_mc');

% coupling strength 0.3
ntrials = 100;
n_neurons = 2;
noise_intensity = 1;
coupling_strength = 0.3;
t = [0:0.01:300];
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

save('run1_coupling03', 'delta_phi_su', 'delta_phi_cd', 'delta_phi_vc', 'delta_phi_mc');

% coupling strength 0.9
ntrials = 100;
n_neurons = 2;
noise_intensity = 1;
coupling_strength = 0.9;
t = [0:0.01:300];
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

save('run1_coupling09', 'delta_phi_su', 'delta_phi_cd', 'delta_phi_vc', 'delta_phi_mc');


