function mask = fermi_filter2(y_size,x_size,sigma,n_kT)

[xx,yy] = meshgrid(1:x_size,1:y_size);
xx = xx - x_size/2;
yy = yy - y_size/2;

r = sqrt(xx.^2 + yy.^2);
kT = sigma/n_kT;
mask = 1./(exp((r-sigma)/kT)+1);
end

