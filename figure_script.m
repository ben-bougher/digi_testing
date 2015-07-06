addpath(genpath('.'))

iter = 500; % maximum number of solves
tolerance = 1e-4; % residual tolerance


% % LSQRcd
mkdir('LSQR')
 int_grad_inv_proj_test_for_ubc('dummypath','800','8','1','36','0', ...
     '0','1',0,1,0,0,0,0,iter, tolerance, 0, 'script_figs/LSQR',0);
movefile('digi_results/*', 'LSQR')
% 
% 
% % LSQR tikonov
 mkdir('LSQR_tik')
 int_grad_inv_proj_test_for_ubc('dummypath','800','8','1','36','0', ...
     '0','1',1,1,0,0,0,0,iter, tolerance, 1, 'script_figs/LSQR_tik',0);
 movefile('digi_results/*', 'LSQR_tik')

% SPGL1
 mkdir('SPGL1')
 int_grad_inv_proj_test_for_ubc('dummypath','800','8','1','36','0', ...
     '0','1',0,1,0,1,0,0,iter, tolerance, 0, 'script_figs/SPGL1',0);
movefile('digi_results/*', 'SPGL1')


% SPGL1 WL
mkdir('SPGL1_WL')
int_grad_inv_proj_test_for_ubc('dummypath','800','8','1','36','0', ...
    '0','1',0,1,0,1,0,1,iter, tolerance, 0, 'script_figs/SPGL1_WL',0);
movefile('digi_results/*', 'SPGL1_WL')



% $$$ iter = 75; % maximum number of solves
% $$$ tolerance = 1e-4; % residual tolerance
% $$$ 
% $$$ 
% $$$ % SPGL1
% $$$ mkdir('SPGL1_early')
% $$$ int_grad_inv_proj_test_for_ubc('dummypath','800','8','1','36','0', ...
% $$$     '0','1',0,1,0,1,0,0,iter, tolerance, 0, 'script_figs/SPGL1_start_stop_early',0);
% $$$ movefile('digi_results/*', 'SPGL1_early')
% $$$ 
% $$$ 
% $$$ % SPGL1 WL
% $$$ mkdir('SPGL1_WL_early')
% $$$ int_grad_inv_proj_test_for_ubc('dummypath','800','8','1','36','0', ...
% $$$     '0','1',0,1,0,1,0,1,iter, tolerance, 0, 'script_figs/SPGL1_WL_early',0);
% $$$ movefile('digi_results/*', 'SPGL1_WL_early')




% $$$ iter = 500; % maximum number of solves
% $$$ tolerance = 1e-4; % residual tolerance

% $$$ 
% $$$ % SPGL1
% $$$ mkdir('SPGL1_pc')
% $$$ int_grad_inv_proj_test_for_ubc('dummypath','800','8','1','36','0', ...
% $$$     '0','1',0,1,0,1,0,0,iter, tolerance, 0, 'script_figs/SPGL1_pc',1);
% $$$ movefile('digi_results/*', 'SPGL1_pc')
% $$$ 
% $$$ 
% $$$ % SPGL1 WL
% $$$ mkdir('SPGL1_WL_pc')
% $$$ int_grad_inv_proj_test_for_ubc('dummypath','800','8','1','36','0', ...
% $$$     '0','1',0,1,0,1,0,1,iter, tolerance, 0, 'script_figs/SPGL1_WL_pc',1);
% $$$ movefile('digi_results/*', 'SPGL1_WL_pc')
% $$$ 
% $$$ iter = 75; % maximum number of solves
% $$$ tolerance = 1e-4; % residual tolerance
% $$$ 
% $$$ 
% $$$ % SPGL1
% $$$ mkdir('SPGL1_pc_early')
% $$$ int_grad_inv_proj_test_for_ubc('dummypath','800','8','1','36','0', ...
% $$$     '0','1',0,1,0,1,0,0,iter, tolerance, 0, 'script_figs/SPGL1_pc_early',1);
% $$$ movefile('digi_results/*', 'SPGL1_pc_early')
% $$$ 
% $$$ 
% $$$ % SPGL1 WL
% $$$ mkdir('SPGL1_WL_pc_early')
% $$$ int_grad_inv_proj_test_for_ubc('dummypath','800','8','1','36','0', ...
% $$$     '0','1',0,1,0,1,0,1,iter, tolerance, 0, 'script_figs/SPGL1_WL_pc_early',1);
% $$$ movefile('digi_results/*', 'SPGL1_WL_pc_early')

