FOI = {};

%%
folder = 'data/' %data folder;
for ii = 1:56
    if ii<10
        filename = ['cropped_Reg_Trans_0',num2str(ii),'_f.mat'];
    else
        filename = ['cropped_Reg_Trans_',num2str(ii),'_f.mat'];
    end
    FOI{size(FOI,1)+1,1} = filename;
end

%%
save filename.mat FOI

















