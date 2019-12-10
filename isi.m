% computes the mean interspike interval and the coefficient of variance
% for a given amount of time, noise, input current, and number of trials
function [mean_cur_isi, std_cur_isi, mean_cv_isi, std_cv_isi] = isi(t, noise, DC_current, ntrials)

    mean_cur_isi = zeros(1, length(DC_current));
    std_cur_isi = zeros(1, length(DC_current));
    mean_cv_isi = zeros(1, length(DC_current));
    std_cv_isi = zeros(1, length(DC_current));

    parfor i = DC_current
        mean_isi_vec = zeros(1, ntrials);
        cv_isi_vec = zeros(1, ntrials);
        for j = 1:ntrials 
            Y = StochasticHH_func(t, @(t) DC_current(i + 1), 0, 100, noise);

            time = Y(:, 1);
            voltage = Y(:, 2);
            spiking_thresh = 40; 
            spike_voltages = zeros(length(time), 1);

            for k = 1:length(time)
                if voltage(k) > spiking_thresh
                    spike_voltages(k) = voltage(k);
                end
                spike_voltages(k) = 0;
            end

            [peaks, idxs] = findpeaks(voltage);
            new_peaks = peaks(peaks > spiking_thresh);
            idxs = idxs(peaks > spiking_thresh);
            spike_times = time(idxs);

            % interspike interval
            interspike_int = diff(spike_times);

            % mean interspike interval
            nspikes = length(interspike_int);
            mean_isi = sum(interspike_int) / nspikes;
            mean_isi_vec(j) = mean_isi;

            % cv of interspike interval 
            cv_isi = sqrt(var(interspike_int)) / mean_isi;
            cv_isi_vec(j) = cv_isi;

        end

        mean_cur_isi(i + 1) = mean(mean_isi_vec);
        std_cur_isi(i + 1) = std(mean_isi_vec);
        mean_cv_isi(i + 1) = mean(cv_isi_vec);
        std_cv_isi(i + 1) = std(cv_isi_vec);

    end 
end