function time_table = auto_sample_mod2(ifr,f_std)

% find local maxima & minima
[max_peak,max_peak_loc] = findpeaks(ifr);
[min_peak,min_peak_loc] = findpeaks(-ifr);
min_peak = -min_peak;

%significant peak
max_peak_loc_sig = max_peak_loc(max_peak>0);

% define sample range
if isempty(max_peak_loc_sig)
    time_table = [];
    %plot
    figure
    subplot(2,1,1)
    plot(ifr)
    xlim([0 length(ifr)])
    subplot(2,1,2)
    plot([0 length(ifr)],[0 0]);
else
    
    delete_peak = [];
    for pp = 1:length(max_peak_loc_sig)
        % find right boundary
        R_find = 0;
        current_peak_loc = max_peak_loc_sig(pp);
        current_peak = ifr(current_peak_loc);
        iter = 1;
        while R_find == 0
            
            temp_R_bound_loc = min_peak_loc(find(min_peak_loc>current_peak_loc,1));
            if isempty(temp_R_bound_loc)
                R_bound_loc(pp) = length(ifr);
                R_find = 1;
                break;
            end
            
            temp_R_bound = ifr(temp_R_bound_loc);
            if temp_R_bound < 0
                R_bound_loc(pp) = current_peak_loc + find(ifr(current_peak_loc+1:end)<0,1);
                R_find = 1;
                break;
            end
            
            right_peak_loc = max_peak_loc(find(max_peak_loc>current_peak_loc,1));
            right_peak = ifr(right_peak_loc);
            
            %             if temp_R_bound/min(current_peak,right_peak) > 0.9
            if temp_R_bound/right_peak > 1
                current_peak = right_peak;
                current_peak_loc = right_peak_loc;
                if ~isempty(find(max_peak_loc_sig==right_peak_loc))
                    delete_peak = [delete_peak, pp+iter];
                    iter = iter+1;
                end
            else
                R_bound_loc(pp) = temp_R_bound_loc;
                R_find = 1;
                break;
            end
        end
        
        % find left boundary
        L_find = 0;
        current_peak_loc = max_peak_loc_sig(pp);
        current_peak = ifr(current_peak_loc);
        while L_find == 0
            
            temp_L_bound_loc = min_peak_loc(find(min_peak_loc<current_peak_loc,1,'last'));
            if isempty(temp_L_bound_loc)
                L_bound_loc(pp) = 1;
                L_find = 1;
                break;
            end
            
            temp_L_bound = ifr(temp_L_bound_loc);
            if temp_L_bound < 0
                L_bound_loc(pp) = find(ifr(1:current_peak_loc)<0,1,'last');
                L_find = 1;
                break;
            end
            
            left_peak_loc = max_peak_loc(find(max_peak_loc<current_peak_loc,1,'last'));
            left_peak = ifr(left_peak_loc);
            
            %             if temp_L_bound/min(current_peak,left_peak) > 0.9
            if temp_L_bound/left_peak > 1
                current_peak = left_peak;
                current_peak_loc = left_peak_loc;
            else
                L_bound_loc(pp) = temp_L_bound_loc;
                L_find = 1;
                break;
            end
        end
    end
    time_table = [L_bound_loc' R_bound_loc'];
    time_table(delete_peak,:) = [];
    max_std = zeros(size(time_table,1),1);
    for ii = 1:size(time_table,1)
        temp_std = f_std(time_table(ii,1):time_table(ii,2));
        max_std(ii) = max(temp_std);
    end
    time_table(max_std<0.2,:) = [];
%     
%     figure
%     subplot(2,1,1)
%     plot(ifr)
%     xlim([0 length(ifr)]);
%     ylim([-2 2]);
%     subplot(2,1,2)
%     hold on
%     for pp = 1:size(time_table,1);
%         H = [pp*(1/size(time_table,1)) 1 1];
%         M = hsv2rgb(H);
%         xx = time_table(pp,1):time_table(pp,2);
%         yy = ifr(xx);
%         area(xx,yy,'EdgeColor',[0 0 0],'FaceColor',M)
%     end
%     xlim([0 length(ifr)]);
%     ylim([-2 2]);
end


