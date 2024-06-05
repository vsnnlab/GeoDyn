clear;
close all;
load('VSDI_profiles_oldWT4AD3_180409.mat');

geo_error = zeros(size(sample_profile,1));
dyn_error = zeros(size(sample_profile,1));
max_error = zeros(size(sample_profile,1));
lat_error = zeros(size(sample_profile,1));
for ii = 1:size(sample_profile,1)
    tic
    for jj = ii:size(sample_profile,1)
        geo_error(ii,jj) = nansum(nansum((sample_profile{ii,3}-sample_profile{jj,3}).^2));
        dyn_error(ii,jj) = nansum(nansum((sample_profile{ii,5}-sample_profile{jj,5}).^2));
        max_error(ii,jj) = nansum(nansum((sample_profile{ii,6}-sample_profile{jj,6}).^2));
        lat_error(ii,jj) = nansum(nansum((sample_profile{ii,7}-sample_profile{jj,7}).^2));
        
        geo_error(jj,ii) = geo_error(ii,jj);
        dyn_error(jj,ii) = dyn_error(ii,jj);
        max_error(jj,ii) = max_error(ii,jj);
        lat_error(jj,ii) = lat_error(ii,jj);
    end
    toc;
end
% geo_error = geo_error/max(geo_error(:));
% geo_error = 1-geo_error;
% dyn_error = dyn_error/max(dyn_error(:));
% dyn_error = 1-dyn_error;
% max_error = max_error/max(max_error(:));
% max_error = 1-max_error;
% lat_error = lat_error/max(lat_error(:));
% lat_error = 1-lat_error;

save VSDI_profile_error_oldWT4AD3_180409.mat geo_error dyn_error max_error lat_error