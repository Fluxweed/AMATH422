% run 4
% mean interspike interval as a function of noise intensity 
% for a fixed coupling strength
ntrials = 10;
kappa = 0.3; 
t = [0:0.01:300];
noise = {'Subunit', 'FoxLuSystemSize', 'VClamp'};
noise_intensity = 0:0.2:7;
n_neurons = 2;
isi_vec_su = zeros(1, length(noise_intensity));
isi_err_su = zeros(1, length(noise_intensity));
isi_vec_cd = zeros(1, length(noise_intensity));
isi_err_cd = zeros(1, length(noise_intensity));
isi_vec_vc = zeros(1, length(noise_intensity));
isi_err_vc = zeros(1, length(noise_intensity));
parfor j = 1:length(noise_intensity)
    n_s = noise_intensity(j);
    [avg_mean_isi_su, avg_mean_isi2_su, err_su, isi2_err_su] = mean_isi_coupled_err(t, noise{1}, ntrials, n_neurons, kappa, n_s);
    isi_vec_su(1, j) = avg_mean_isi_su;
    isi_err_su(1, j) = err_su;
    [avg_mean_isi_cd, avg_mean_isi2_cd, err_cd, isi2_err_cd] = mean_isi_coupled_err(t, noise{2}, ntrials, n_neurons, kappa, n_s);
    isi_vec_cd(1, j) = avg_mean_isi_cd;
    isi_err_cd(1, j) = err_cd;
    [avg_mean_isi_vc, avg_mean_isi2_vc, err_vc, isi2_err_vc] = mean_isi_coupled_err(t, noise{3}, ntrials, n_neurons, kappa, n_s);
    isi_vec_vc(1, j) = avg_mean_isi_vc;
    isi_err_vc(1, j) = err_vc;
end 

save('run4_isis_witherr', 'isi_vec_su', 'isi_vec_cd', 'isi_vec_vc', ...
    'isi_err_su', 'isi_err_cd', 'isi_err_vc');

