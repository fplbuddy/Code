%
clear all
%
n = 32;
kmax = n/3;
%
for i = 1:n/2
   k(i) = i-1;
   k(i+n/2) = i-n/2-1;
end
for i = 1:n/2
    ks(i) = k(i+n/2);
    ks(i+n/2) = k(i);
end
%
% Data
%
fid = fopen('/mi/share/scratch/dallas/Coding/HD/2D_box/HD2Df_v3/Spectra/spectrum2D_UU.000.dat','r');
fread(fid,1,'real*4');
b3 = fread(fid,inf,'real*8');
fclose(fid);
%
b3 = reshape(b3,n,n);
%
for j = 1:n
for i = 1:n/2
    bs(i,j) = b3(i+n/2,j);
    bs(i+n/2,j) = b3(i,j);
end
end
%
for j = 1:n/2
for i = 1:n
    bs2(i,j) = bs(i,j+n/2);
    bs2(i,j+n/2) = bs(i,j);
end
end
% %
% for j = 1:n/2-1
%     bs3(:,n/2-j) = flipud(bs2(:,j+n/2));
% end
% for j = 1:n/2
%     bs3(:,j+n/2) = bs2(:,j+n/2);
% end
stop
%%
% Plots
%
figure
% contour(ks,ks,log10(bs2));
pcolor(ks,ks,log10(bs2));
% pcolor(log10(bs3));
% shading flat
% shading interp
colorbar
% caxis([1e-6 1e-4])
% axis([-kmax kmax+1 -kmax kmax+1])
% axis([-16 16 -16 16])
xlabel('\fontsize{16}\bf k_x')
ylabel('\fontsize{16}\bf k_y')
stop