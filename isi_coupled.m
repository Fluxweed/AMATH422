% returns the mean interspike interval for a single voltage trace given
% the time and the voltages
function [mean_isi, spike_times] = isi_coupled(time, voltage)        

    spiking_thresh = 40; 
    spike_voltages = zeros(length(time), 1);

    for k = 1:length(time)
        if voltage(k) > spiking_thresh
            spike_voltages(k) = voltage(k);
        end
        spike_voltages(k) = 0;
    end

    % find spike times
    [peaks, idxs] = findpeaks(voltage);
    new_peaks = peaks(peaks > spiking_thresh);
    idxs = idxs(peaks > spiking_thresh);
    spike_times = time(idxs);

    % interspike interval
    interspike_int = diff(spike_times);

    % mean interspike interval
    nspikes = length(interspike_int);
    mean_isi = sum(interspike_int) / nspikes;
end