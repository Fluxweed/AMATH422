% run 4
% mean interspike interval as a function of noise intensity 
% for a fixed coupling strength
ntrials = 10;
kappa = 0.3; 
t = [0:0.01:300];
noise = {'Subunit', 'FoxLuSystemSize', 'VClamp', 'MarkovChain'};
noise_intensity = 0:0.3:7;
n_neurons = 2;
isi_vec_subunit = zeros(1, length(noise_intensity));
isi_vec2_subunit = zeros(1, length(noise_intensity));
isi_vec_cond = zeros(1, length(noise_intensity));
isi_vec2_cond = zeros(1, length(noise_intensity));
isi_vec_vc = zeros(1, length(noise_intensity));
isi_vec2_vc = zeros(1, length(noise_intensity));
%isi_vec_mc = zeros(1, length(noise_intensity));
%isi_vec2_mc = zeros(1, length(noise_intensity));
parfor j = 1:length(noise_intensity)
    n_s = noise_intensity(j);
    [avg_mean_isi_su, avg_mean_isi2_su] = mean_isi_coupled(t, noise{1}, ntrials, n_neurons, kappa, n_s);
    isi_vec_subunit(1, j) = avg_mean_isi_su;
    isi_vec2_subunit(1, j) = avg_mean_isi2_su;
    [avg_mean_isi_cd, avg_mean_isi2_cd] = mean_isi_coupled(t, noise{2}, ntrials, n_neurons, kappa, n_s);
    isi_vec_cond(1, j) = avg_mean_isi_cd;
    isi_vec2_cond(1, j) = avg_mean_isi2_cd;
    [avg_mean_isi_vc, avg_mean_isi2_vc] = mean_isi_coupled(t, noise{3}, ntrials, n_neurons, kappa, n_s);
    isi_vec_vc(1, j) = avg_mean_isi_vc;
    isi_vec2_vc(1, j) = avg_mean_isi2_vc;
    %[avg_mean_isi_mc, avg_mean_isi2_mc] = mean_isi_coupled(t, noise{4}, ntrials, n_neurons, kappa, n_s);
    %isi_vec_mc(1, j) = avg_mean_isi_mc;
    %isi_vec2_mc(1, j) = avg_mean_isi2_mc;
end 

save('run4_isis_new', 'isi_vec_subunit', 'isi_vec_cond', 'isi_vec_vc', ...
    'isi_vec2_subunit', 'isi_vec2_cond', 'isi_vec2_vc');

