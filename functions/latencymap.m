function latency_map = latencymap(sample)
f_temp_sample = sample;
latency_map = zeros(size(f_temp_sample,1),size(f_temp_sample,2));
for cc = 1:size(f_temp_sample,1);
    for rr = 1:size(f_temp_sample,2);
        temp_signal = squeeze(f_temp_sample(cc,rr,:));
        complex_signal = hilbert(temp_signal);
        Amp = abs(complex_signal);
        phase = angle(complex_signal);
        latency_cand = find(abs(diff(phase./abs(phase)))==2);
        if length(latency_cand) == 1
            latency_map(cc,rr) = latency_cand;
        end
    end
end
[~,start_pp] = max(latency_map(:));
[start_cc,start_rr] = ind2sub(size(latency_map),start_pp);

for cc = start_cc:size(f_temp_sample,1);
    for rr = start_rr:size(f_temp_sample,2);
        temp_signal = squeeze(f_temp_sample(cc,rr,:));
        complex_signal = hilbert(temp_signal);
        Amp = abs(complex_signal);
        phase = angle(complex_signal);
        latency_cand = find(abs(diff(phase./abs(phase)))==2);
        ref_latency_list = latency_map(max(1,cc-1):min(size(latency_map,1),cc+1),max(1,rr-1):min(size(latency_map,2),rr+1));
        ref_latency = nanmean(ref_latency_list(ref_latency_list~=0));
        [~,lat_ii] = min(abs(latency_cand-ref_latency));
        latency = latency_cand(lat_ii);
        if isempty(latency)
            latency = 0;
        end
        latency_map(cc,rr) = latency;
    end
    for rr = start_rr:-1:1;
        temp_signal = squeeze(f_temp_sample(cc,rr,:));
        complex_signal = hilbert(temp_signal);
        Amp = abs(complex_signal);
        phase = angle(complex_signal);
        latency_cand = find(abs(diff(phase./abs(phase)))==2);
        ref_latency_list = latency_map(max(1,cc-1):min(size(latency_map,1),cc+1),max(1,rr-1):min(size(latency_map,2),rr+1));
        ref_latency = nanmean(ref_latency_list(ref_latency_list~=0));
        [~,lat_ii] = min(abs(latency_cand-ref_latency));
        latency = latency_cand(lat_ii);
        if isempty(latency)
            latency = 0;
        end
        latency_map(cc,rr) = latency;
    end
end
for cc = start_cc:-1:1;
    for rr = start_rr:size(f_temp_sample,2);
        temp_signal = squeeze(f_temp_sample(cc,rr,:));
        complex_signal = hilbert(temp_signal);
        Amp = abs(complex_signal);
        phase = angle(complex_signal);
        latency_cand = find(abs(diff(phase./abs(phase)))==2);
        ref_latency_list = latency_map(max(1,cc-1):min(size(latency_map,1),cc+1),max(1,rr-1):min(size(latency_map,2),rr+1));
        ref_latency = nanmean(ref_latency_list(ref_latency_list~=0));
        [~,lat_ii] = min(abs(latency_cand-ref_latency));
        latency = latency_cand(lat_ii);
        if isempty(latency)
            latency = 0;
        end
        latency_map(cc,rr) = latency;
    end
    for rr = start_rr:-1:1;
        temp_signal = squeeze(f_temp_sample(cc,rr,:));
        complex_signal = hilbert(temp_signal);
        Amp = abs(complex_signal);
        phase = angle(complex_signal);
        latency_cand = find(abs(diff(phase./abs(phase)))==2);
        ref_latency_list = latency_map(max(1,cc-1):min(size(latency_map,1),cc+1),max(1,rr-1):min(size(latency_map,2),rr+1));
        ref_latency = nanmean(ref_latency_list(ref_latency_list~=0));
        [~,lat_ii] = min(abs(latency_cand-ref_latency));
        latency = latency_cand(lat_ii);
        if isempty(latency)
            latency = 0;
        end
        latency_map(cc,rr) = latency;
    end
end
