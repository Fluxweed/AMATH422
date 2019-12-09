function [isi_matrix, isi_matrix2] = isi_matrix(t, noise, ntrials, n_neurons, noise_intensity, coupling_strength)

    isi_matrix = zeros(length(coupling_strength), length(noise_intensity));
    isi_matrix2 = zeros(length(coupling_strength), length(noise_intensity));
    for i = 1:length(coupling_strength)
        for j = 1:length(noise_intensity)
            kappa = coupling_strength(i);
            n_s = noise_intensity(j);
            [avg_mean_isi, avg_mean_isi2] = mean_isi_coupled(t, noise, ntrials, n_neurons, kappa, n_s);
            isi_matrix(i, j) = avg_mean_isi;
            isi_matrix2(i, j) = avg_mean_isi2;
        end 
    end 
end