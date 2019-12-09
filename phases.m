% returns the phase and spike times for a single neuron across the 
% given time interval
function [phi, spike_times] = phases(t, voltage)

    [mean_isi, spike_times] = isi_coupled(t, voltage);
    
    spike_times = [0 spike_times];
    nspikes = length(spike_times);
    phi = zeros(1, length(t(t < spike_times(end-1))));
    start = 0;
    total_time = 1;
    for i = 0:(nspikes-2)
        times = t(t < spike_times(i+1)); 
        times = times(times >= start);
        for j = 1:length(times)
            p = 2*pi*i + 2*pi*((times(j) - spike_times(i+1)) / (spike_times(i+2) - spike_times(i+1)));
            phi(total_time) = p;
            total_time = total_time + 1;
        end 
        start = spike_times(i+1);
    end 
    
    difs = diff(spike_times);
    for i = 1:length(difs)
        if difs(i) < 1
            spike_times(i) = 0;
        end
    end
    spike_times(spike_times == 0) = [];
end