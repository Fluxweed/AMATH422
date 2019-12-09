% returns a matrix containing the phase offsets corresponding to n_neurons,
% and a given type of noise, coupling strength, and noise intensity.
function [delta_phi, lengths] = multiple_phase_offsets(t, n_neurons, coupling_strength, noise_intensity, noise)

    Y = noise_coupling(n_neurons, coupling_strength, noise_intensity, t, @(t) 10, 0.5, 100, noise);
    time = Y.t;

    M = zeros(length(time),n_neurons);
    x = zeros(1, n_neurons);
    for i = 1:n_neurons
        [phi, spike_times] = phases(time, Y.V(:,i));
        x(i) = length(phi);
        M(1:length(phi), i) = phi;
    end

    lengths = zeros(1, n_neurons - 1);
    for i = 1:length(x)-1
        lengths(i) = min(x(i), x(i+1));
    end

    delta_phi = zeros(length(time),n_neurons-1);
    for i = 1:n_neurons-1
        phi = M(1:lengths(i), i);
        phi2 = M(1:lengths(i), i+1);
        delta_phi(1:lengths(i), i) = mod((phi - phi2), 2*pi);
    end
end