addpath(genpath('.'))

iter = 500; % maximum number of solves
tolerance = 1e-4; % residual tolerance
plot_stack = 0;
plot_convergence = 1;
plot_wavelet = 0;
precond = 0;


% % LSQR
smoothing = 0;
Wdomain = 0;
l1=0;
prefix = 'script_figs/LSQR';


mkdir('LSQR')
 IG = int_grad_inv_proj_test_for_ubc('dummypath','800','8','1','36','0', ...
     '0','1',plot_stack,plot_convergence,plot_wavelet, l1, Wdomain,...
                                iter, tolerance, smoothing, prefix,0);
movefile('digi_results/*', 'LSQR')
save 'LSQR/IG.mat', IG
% 
% 


% % LSQR tikonov
smoothing = 1;
prefix = 'script_figs/LSQR_tik';

 mkdir('LSQR_tik')
 IG = int_grad_inv_proj_test_for_ubc('dummypath','800','8','1','36','0', ...
     '0','1',plot_stack,plot_convergence,plot_wavelet, l1, Wdomain,...
                                iter, tolerance, smoothing, prefix,0);
 movefile('digi_results/*', 'LSQR_tik')
    save 'LSQR_tik/IG.mat', IG

% SPGL1
 mkdir('SPGL1')
 smoothing = 0;
 l1 = 1;
 prefix = 'script_figs/SPGL1';
 
 IG = int_grad_inv_proj_test_for_ubc('dummypath','800','8','1','36','0', ...
     '0','1',plot_stack,plot_convergence,plot_wavelet, l1, Wdomain,...
                                iter, tolerance, smoothing, prefix,0);
movefile('digi_results/*', 'SPGL1')
save 'SPGL1/IG.mat', IG

% SPGL1 WL
prefix =  'script_figs/SPGL1_WL';
Wdomain = 1;

mkdir('SPGL1_WL')
 IG = int_grad_inv_proj_test_for_ubc('dummypath','800','8','1','36','0', ...
     '0','1',plot_stack,plot_convergence,plot_wavelet, l1, Wdomain,...
                                iter, tolerance, smoothing, prefix,0);
movefile('digi_results/*', 'SPGL1_WL')
save 'SPGL1_WL/IG.mat', IG




