function phase_lat = latencymap2(sample)
phase_lat = zeros(101,101);
max_time = size(sample,3);
for cc = 1:101
    for rr = 1:101
        temp_signal = squeeze(sample(cc,rr,:));
        complex_signal = hilbert(temp_signal);
        Amp = abs(complex_signal);
        phase = angle(complex_signal);
        lat = find(phase == 0);
        lat(lat==1) = [];
        lat(lat==max_time) = [];
        if ~isempty(lat)
            if length(lat) ==1
                phase_lat(cc,rr) = lat;
            else
                phase_lat(cc,rr) = 0;
            end
        else
            sign_phase = sign(phase);
            lat = find(abs(diff(sign_phase)) ==2)+1;
            trans_point = find(abs(diff(phase))>2.5)+1;
            for ii = 1:length(trans_point)
                lat(lat == trans_point(ii)) =[];
            end
            lat(lat==1) = [];
            lat(lat==max_time) = [];
            if length(lat) ==1
                phase_lat(cc,rr) = lat;
            else
                phase_lat(cc,rr) = 0;
            end
        end
    end
end
