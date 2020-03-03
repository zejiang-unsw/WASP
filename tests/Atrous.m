%Toolbox to use
addpath(genpath("c:\Program Files\MATLAB\R2019a\toolbox\Wavelab850\"))

%CAPTION
fprintf('\n');
disp('Figure 4 JoH')
disp('The top curve is the original signal.')
disp('Its dyadic wavelet transform is shown below at scales 2^j for 0<j<6.')
disp('The bottom curve carries the lower frequencies corresponding ')
disp('to scales larger than $2^5$.')

  close all;
  global WLVERBOSE;
  WLVERBOSE='No';
  load('C:\Users\z5154277\OneDrive - UNSW\R_Package\WASP\data-raw\x.mat');
  L = 1;
  dwt = FWT_ATrou(x,L);
  figure(1);clf
  DisplayDWT(dwt,L-1,x);
  WLVERBOSE='Yes';
  
  save('C:\Users\z5154277\OneDrive - UNSW\R_Package\WASP\data-raw\x-atrous.mat', 'dwt')

% Written by Maureen Clerc and Jerome Kalifa, 1997
% clerc@cmapx.polytechnique.fr, kalifa@cmapx.polytechnique.fr
   D = size(dwt,2); 
   n = size(x,1);
   idwt = zeros(n,D);
   for d = 0:D-1
        dwt_i = FWT_ATrou(x,d);
%         p = IWT_ATrou(dwt_i,d)';
%         idwt(:,d+1) = p;
        idwt(:,d+1) = dwt_i(:,1);   
   end
   
   
% Plot observation and each component
    figure
    subplot(D+1,1,1);
    plot(x);
    ylabel('obs');
    axis([1 n min(x) max(x)])
    for is=1:D
        subplot(D+1,1,is+1);
        plot(idwt(:,D-is+1));
        ylabel(['c',int2str(is-1)]);
        %axis([1 n min(idwt(:,is)) max(idwt(:,is))])
        axis([1 n min(x) max(x)])
    end
   
    save('C:\Users\z5154277\OneDrive - UNSW\R_Package\WASP\data-raw\x-atrous-mra.mat', 'idwt')
    

%  Part of Wavelab Version 850
%  Built Tue Jan  3 13:20:39 EST 2006
%  This is Copyrighted Material
%  For Copying permissions see COPYING.m
%  Comments? e-mail wavelab@stat.stanford.edu 
