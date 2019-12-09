% computes the mean interspike intervals for two coupled neurons given 
% the time, type of noise, the number of trials, and the coupling strength
% and noise intensity
function [avg_mean_isi, avg_mean_isi2] = mean_isi_coupled(t, noise, ntrials, n_neurons, coupling_strength, noise_intensity)
    
    mean_isi_vec = zeros(1, ntrials);
    mean_isi_vec2 = zeros(1, ntrials);
    
    for j = 1:ntrials 
        
        % Y = coupled_HH_func(t, @(t) 10, 0.5, 100, noise, coupling_strength, noise_intensity);
        Y = noise_coupling(n_neurons, coupling_strength, noise_intensity, t, @(t) 10, 10, 100, noise);
        
        time = Y.t;
        voltage = Y.V(:,1);
        [mean_isi, spike_times] = isi_coupled(time, voltage);
        mean_isi_vec(j) = mean_isi;
                
        time2 = Y.t;
        voltage2 = Y.V(:,2);
        [mean_isi2, spike_times2] = isi_coupled(time2, voltage2);
        mean_isi_vec2(j) = mean_isi2;

    end
    avg_mean_isi = sum(mean_isi_vec) / length(mean_isi_vec);
    avg_mean_isi2 = sum(mean_isi_vec2) / length(mean_isi_vec2);
end