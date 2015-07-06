function [] = int_grad_inv_proj_test_for_ubc(job_meta_path,i_block,...
    startvol,volinc,endvol,tottracerun,maxzout,wavevar, plot_stack, ...
    plot_convergence, plot_wavelet, l1, Wdomain, max_iterations, ...
    tolerance, smoothing, prefix, pc)

%% ------------------ FUNCTION DEFINITION ---------------------------------
% INT_GRAD_INV_PROJ: function to run Intercept Gradient Inversion using
% dynamic wavelet set.
%   Inputs:
%       job_mat_path = path of metadata .mat file.
%       i_block = current block to be processed.
%       startvol = number of the first angle trace/volume to read
%       volinc = angle trace/volume increment
%       endvol = number of the last angle trace/volume to read
%       tottracerun = number of traces to output, 0 = all, and il1234 would
%       only output inline 1234
%       maxzout = To be given in samples if this is initialized as 0 this will be
%       converted to the maximum number of samples in the volume
%       wavevar = Flag 1: for spatially varying wavelet, Flag 0: for
%       spatially stationary wavelet
%
%   Additional behaviour and plotting flags
%       plot_stack: Flag to plot the stacked section
%       plot_convergence: Flag to plot the convergence curve of the
%                         solver
%       plot_wavelet: Flag to plot the wavelet
%       l1: Flag to use SPGL1 as a solver
%       Wdomain: Flag to regularize in the Discrete Wavelet domain.
%       max_iterations: Maximimum number of iterations to use in
%                       the inversion.
%       tolerance: Tolerance to use for least squares solution.
%       smoothing: Flag specifying that smoothing via Tikonov
%                  regularization.
%       prefix: File prefix for plotting output.
%       pc: Flag to use a preconditioner matrix in the inversion.
%
%
%   Outputs:
%       digi_intercept SEGY files.
%       digi_gradient SEGY files.
%       digi_minimum_energy_eer_projection SEGY files.
%
% Authors: James Selvage, Jonathan Edgar and a bit of Charles Jones and Simon Wrigley
% -------------------------------------------------------------------------




%% Parameters
% would normaly convert all parameters to double, but keep i_block as string as being passed to
% other modules; it does happen at the bottom of this program for output
%i_block = str2double(i_block);
%
% angle trace data ranges to use, vol is an angle trace either as a
% seperate angle volume or as a angle trace in an angle gather

startvol = str2double(startvol);    % number of the first angle trace/volume to read
volinc = str2double(volinc);        % angle trace/volume increment 
%endvol = job_meta.nvols;
endvol = str2double(endvol);        % number of the last angle trace/volume to read

% number of traces to run, put to zero to make it run all traces in the
% block, this is the default, this is also used to pass an inline (pkey
% number to use in testing has to be ilnnnn format
useselectemode = 0;

if isempty(regexp(tottracerun,'il','once')) == 0
    useselectemode = 1;
    requiredinline =  str2double(regexprep(tottracerun,'il',''));
    %requiredinline =  str2double(strrep(tottracerun,'il',''));
    tottracerun = 0;
else
    tottracerun = str2double(tottracerun);
end
% tottracerun = 500;
%maxzout = 8000;
maxzout = str2double(maxzout);

output_std = 0;     % do you want to calculate and output standard line fitting intercept and gradient and eer
plot_on = 0;        % do you want interactive plots to pop up when running (just for debug)

% do you want to have a background starting model, default 0 is all zeros,
% setting this to value 1 uses standard line fitting to make a background
% model; lsq seems to reach the same answer with any model , did try random
% numbers and no difference in result just more iterations required
background = 0;
needconf = 0;                   % do you want a confidence volume calculating and outputing, 0 = no , 1 = yes
noofreports = 20;               % how many % progress reports do you want in the running of the job, 5 = 20% increment for reporting
chi_model_type = 'empirical';   % what is the chi model that you want for the eer...
                                % ...empirical is 19 degrees with 1 degree increase per km or second

% Tikhonov regularisation weight, 1 = little smoothing, 10 = moderate,
% 100 = smooth, 1000 = very smooth
%wsmooth = 7000; uruguay gathers
%wsmooth = 10; tza angle gathers
%wsmooth = str2double(wsmooth);

eer_weight = 0.1;   % Weight for EER constraint between 0 and 1, value of 1 forces it much  closer to the intercept values
iter = 300;         % maximum number of iterations in the lsq solver - normally 300
tol = 1e-3;         % convergence tolerance in the lsq solver can move to tol = 5e-4 or 1e-4.....but larger value speeds it up as less iterations to run 
%tol = 8e-4;  
padding = 50;       % water bottom pick padding
extrapad = 0;     % extra padding to top of dataset, this many samples will get zeroed below the wb 
%extrapad = (floor(extrapad/32))*32;
use_spatial_wavelets = wavevar; % FLag 1 for using spatial wavelets

warning off all;    % to reduce printout in compilied version turned all warning off

relfreq = 6;        % frequency to taper down below in hz, it is a sine taper so gentle ramp off, so can be a bit higher than expected, ie 6 cuts in at 4.5Hz
lowfreqtaperstart = 6; % Hz to start tapering the wavelet down from to 0 hz
wsmo_scal = 5; % scaler to make the weighting in the tikonov regularisation , normally set to 5 * the std dev of the amplitudes of the input data as measured
               % in the wavelet estimation code  

               
load('ubc_int_grad_test1.mat');               

% Plot traces and the stack
if plot_stack==1
    stack = squeeze(sum(vol_traces, 2));
    figure
    stack_std = std(stack(:));
    imagesc(stack);caxis([-stack_std stack_std]);
    colormap(gray);title('Stacked Seismic');
    xlabel('Trace #');ylabel('samples');
    saveas(gcf,[prefix,'stacked_seismic'],'epsc')
    
    figure
    
    plot(stack(:,2), 1:size(stack(:,5),1));title('Stacked trace (5)');
    axis ij
    xlabel('Amplitude');ylabel('samples');
    saveas(gcf,[prefix,'stacked_seismic_trace.eps'],'epsc')

    stack_spectra = 20.*log10(abs(fft(stack)));
    
    figure;
    plot(linspace(0,pi,ns/2),stack_spectra(1:end/2,5));title('Stack Amplitude Spectra (trace 5)');
    ylabel('Amplitude dB [ref arb]'); xlabel('w [radians]');
    saveas(gcf,[prefix, 'stacked_spectra.eps'],'epsc')
    
    % do a wigliorama
    wigglePlot(stack, [prefix, 'stacked_seismic']);
end
    
    

% end of parameters

%%
ns_wavelet = size(wavelets.all_wavelets_time{1},1)-1;                               % Number of samples in a wavelet?
hns_wavelet = floor(ns_wavelet/2);                                                  % Half of number of samples?

%=======================================
%make a low freq taper for the wavelet to avoid blowing up noise

taperlen = floor((lowfreqtaperstart/(500000/job_meta.s_rate))*hns_wavelet)+1;
taperst = (sin(linspace((-(pi*0.94)/2),((pi*0.8)/2),taperlen)')+1)/2;
%taperst(1:(end-2)) = taperst(1:(end-2)) ./ 10;
taperst= power(taperst,4);
taperend = flipud(taperst);
taperapply = [taperst;ones((ns_wavelet-(taperlen*2)),1);taperend];

%==========================================

vol_count = 1;
%Loop though the different angle volumes
for i_vol = startvol:volinc:endvol
    wavelet_z_grid = wavelets.all_wavelets_freq{i_vol}(1,:);
    wavelet_z_grid = wavelet_z_grid + padding;
    wavelet = wavelets.all_wavelets_freq{i_vol}(2:end,:);
    wavelet(isnan(wavelet)) = 0;
    % apply the low freq taper
    wavelet = bsxfun(@times,wavelet,taperapply);
    if plot_on == 1
        figure(2)
        subplot(1,totalvol,vol_count); imagesc(wavelet);
        figure(3)
        subplot(1,totalvol,vol_count); plot(wavelet(:,10));
    end
    start_interp = min(wavelet_z_grid);%-hns_wavelet;
    end_interp = max(wavelet_z_grid);%+hns_wavelet;
    wavelet_interp{vol_count} = interp1(wavelet_z_grid,wavelet',start_interp:1:end_interp,'linear');
    %wavelet_interp{vol_count} = interpolate_wavelets(wavelet_z_grid,wavelet,start_interp,end_interp);
    wavelet_interp{vol_count} = circshift(ifft(wavelet_interp{vol_count}','symmetric'),floor(job_meta.ns_win/2));
    %wavelet_interp{i_vol} = wavelet_interp{i_vol};
    
    if start_interp > 1;
        % Pad with zeros or first wavelet
        pad_size = start_interp-1;
        wavelet_interp{vol_count} = [repmat(wavelet_interp{vol_count}(:,1),1,pad_size),...
            wavelet_interp{vol_count}];
    end
    if end_interp < ns;
        % Pad with zeros or last wavelet
        pad_size = ns-end_interp;
        wavelet_interp{vol_count} = [wavelet_interp{vol_count},...
            repmat(wavelet_interp{vol_count}(:,end),1,pad_size)];
    end
    if floor(ns_wavelet/2) == ns_wavelet/2
        wavelet_interp{vol_count}(end+1,:) = wavelet_interp{vol_count}(end,:);
    end
    if plot_on == 1
        figure(4)
        subplot(1,totalvol,vol_count); imagesc(wavelet_interp{vol_count});
    end
    wavelet_interp{vol_count} = wavelet_interp{vol_count}(:,1:ns);
    vol_count = vol_count + 1;    
end
interp_wavelet_z_grid = 1:ns;
ns_wavelet = size(wavelet_interp{1},1);

%==========================================================================================================================
% normalise the wavelets
if job_meta.is_gather == 0
    [~, wbliveidx] =  min(abs((startvol:volinc:endvol)-pick_wb_ind));
    %wavelet_norm = wavelet_rms_norm(totalvol,wavelet_interp,interp_wavelet_z_grid,ceil(totalvol*0.6667));
    wavelet_norm = wavelet_rms_norm(totalvol,wavelet_interp,wbliveidx);
else
    % need to store the water bottom live offset in the job_meta
    %job_meta.livewb = 12;
    %[~, wbliveidx] = min(abs(input_angles-job_meta.livewb));
    [~, wbliveidx] =  min(abs((startvol:volinc:endvol)-pick_wb_ind));
    wavelet_norm = wavelet_rms_norm(totalvol,wavelet_interp,wbliveidx);
end

clear wavelet_interp wavelets wavelet_z_grid interp_wavelet_z_grid wavelet ;

for i_vol = 1:1:totalvol
    wavelet_norm{i_vol}(isnan(wavelet_norm{i_vol})) = 0;
end

if plot_wavelet == 1
    
    figure
    subplot(211);
    plot(mean(waveletck,2));xlabel('samples');ylabel('amplitude');title('Mean Wavelet');
    subplot(212) 
    spectra = 20*log10(abs(fft(mean(waveletck,2))));
    plot(linspace(0,pi,size(spectra,1)/2),spectra(1:end/2));
        xlabel('w [rad]');ylabel('Amplitude [dB ref arb]');
        
    saveas(gcf,prefix+'estimated_wavelet.eps','epsc')
    
%     imagesc(waveletck(:,1000:ns:end));
%     plot(sum(abs(waveletck(:,1600:ns:end)),1));



end


%%
% start building the inversion operators
% Chi model
switch chi_model_type
    case 'empirical' 
        chi = (job_meta.s_rate/1e6)*(0:1:ns-1)'.*-2 + 19;

    case 'raw'
        
    case 'bootstrap'

    otherwise
        warning('Unexpected plot type. No plot created.');
end

% Build operator matrix
% Build blanking matrix used to ensure the convolution operator matrix is correct
IGblank = spdiags(ones(ns,2*hns_wavelet+1),(-hns_wavelet:hns_wavelet),ns,ns);
IGblank = repmat(IGblank,1+totalvol,2);

% Tikhonov regularisation matrix % check tik order
%##smooth = spdiags([-wsmooth*ones(2*ns,1) wsmooth*ones(2*ns,1)],[-1 0],2*ns,2*ns);
%smooth = spdiags([-wsmooth*ones(2*ns,1) wsmooth*ones(2*ns,1)],[-1 0],ns,ns);

smooth = spdiags([wsmooth*ones(2*ns,1) -2*wsmooth*ones(2*ns,1) wsmooth*ones(2*ns,1)],[-1 0 1],2*ns,2*ns);
if l1 ==1 | smoothing==0
    smooth = smooth.*0;
end
%smooth = smooth.*0;
%cj edit
%smooth = spdiags([wsmooth*ones(2*ns,1) -2*wsmooth*ones(2*ns,1) wsmooth*ones(2*ns,1)],[-1 0 1],ns,ns);
%smooth2 = spdiags([wsmooth*10*ones(2*ns,1) -2*wsmooth*10*ones(2*ns,1) wsmooth*10*ones(2*ns,1)],[-1 0 1],ns,ns);
% smooth = [smooth, smooth];

% Extract the angle of the stack closest to the mean chi angle for this trace
% This is for the IG crossplot correlation constraint
theta = asind(sqrt(tand(mean(chi))));
[~, angle_index] = min(abs(input_angles-theta));

IGmatrix = build_operator(totalvol,input_angles,ns_wavelet,wavelet_norm,ns,hns_wavelet,angle_index,chi,eer_weight,IGblank,smooth);

clear wavelet_norm IGblank smooth ;


%% Inversion loop

% Set first and last traces to loop over
first_iter = 1;
if tottracerun ~= 0;
    ntraces = tottracerun;
end
last_iter = ntraces;

% Begin inversion loop
tic
ava = zeros(2*ns,((last_iter - first_iter)+1),'single');
%residcount = zeros(iter,(last_iter - first_iter)+1);
%cjmodel = rand(4100,1);
%
% make a taper to apply before any filterinf
tflen = 0.05;
%zeropad = 5;
up = ((sin(linspace((-(pi)/2),((pi)/2),round(tflen*ns))')+1)/2);
%down = flipud(up);
%outtaper = [up; ones(ns-(length(up)*2)-zeropad,1); down; zeros(zeropad,1); up; ones(ns-(length(up)*2)-zeropad,1); down; zeros(zeropad,1) ];
down = ((sin(linspace(((pi)/2),(-(pi)/2),round((tflen*2)*ns))')+1)/2);
outtaper = [up; ones(ns-(length(up)+length(down)),1); down; up; ones(ns-(length(up)+length(down)),1); down ];


% now loop round all the traces in the input and run the inversion we have
% added a parfor loop for testing but will not compile
%for kk = 2 % first_iter:last_iter 
for kk = first_iter:last_iter
%     for ii = 1:totalvol
%         data(:,ii) = traces{ii}(:,kk); % Read the angle stack data for the inversion of this trace
%     end

if useselectemode == 0
    requiredinline = ilxl_read{vol_index_wb}(kk,1); 
end

if ilxl_read{vol_index_wb}(kk,1) == requiredinline;

    if job_meta.is_gather == 1
        data = vol_traces(:,:,kk);
    else
        data = vol_traces(:,:,kk)';
    end
    data(isnan(data)) = 0;  % Set NaNs to zero
      
    fold = sum(data ~= 0,2); % Get angle fold  
    
    
    
    
    fmask = data == 0;
    %fmask = low_amp_mute(trim_data);

    % filter the low freq out of the dataset
        %[dataqc scaleused] = time_balence(data);
        %figure(57); imagesc(data); colormap(gray); caxis([-10000 10000]); 
        %amp_spec_plot(data,job_meta.s_rate)
    data = bandpass_filter(data,(job_meta.s_rate/1000000),0,relfreq,top3dpt,topfreq);
        %amp_spec_plot(data,job_meta.s_rate)
        %[dataqc scaleused] = time_balence(dataqc);
        %figure(58); imagesc(data); colormap(gray); caxis([-10000 10000]); 
    data = data .* (1.-fmask);
    %figure(59); imagesc(data); colormap(gray);
    data_tmp = data(:);  % Make temporary column vector from data
    
    %data_zeros = abs(data_tmp) < abs(mean(data_tmp)/range(data_tmp)); % Find zones where data is zero (due to mute angle mute functions)
    data_zeros = data_tmp == 0;
    data_zeros2 = data(:,wbliveidx) == 0;
    % or close to zero
    %data(reshape(data_zeros,[],totalvol)) = 0;
    %data_zeros = logical([data_zeros;zeros(3*ns,1)]); 
    %data_zeros = logical([data_zeros;data_zeros2;zeros(2*ns,1)]);
    data_zeros = logical([data_zeros;data_zeros2;data_zeros2;data_zeros2]); 
    
    if plot_on == 1;
        figure(6)
        subplot(1,2,1); spy(IGmatrix)
    end
    IGiter = IGmatrix;
    IGiter(data_zeros,:) = 0; % Set operator rows to zero if there are zeros in the data vector
    
    if plot_on == 1;
        figure(6)
        subplot(1,2,2); spy(IGiter)
    end
    
    % Model intercept = [near], model gradient = [far-near]/[far_angle - near_angle]
    if background == 1 || output_std == 1;
        model_tmp = zeros(2,ns);
        for ii = 1:ns
            model_op = [ones(totalvol,1),sind(input_angles').*sind(input_angles')];
            model_zeros = data(ii,:) == 0;
            model_op(model_zeros,:) = 0;
            model_tmp(:,ii) = model_op\data(ii,:)';
        end
        model_tmp = model_tmp';
        %[traces{1}(:,1) traces{2}(:,1) traces{3}(:,1) model_tmp]; 
        if output_std == 1;
            Imodel(:,kk) = model_tmp(:,1)/norm(model_tmp(:,1)); %data(:,1)/norm(data(:,1));
            Gmodel(:,kk) = model_tmp(:,2)/norm(model_tmp(:,1)); %-Imodel./tand(chi);
            model = [Imodel(:,kk);Gmodel(:,kk)];
            model(isnan(model)) = 0;
        else
            Imodel = model_tmp(:,1)/norm(model_tmp(:,1)); %data(:,1)/norm(data(:,1));
            Gmodel = model_tmp(:,2)/norm(model_tmp(:,1)); %-Imodel./tand(chi);
            %model = [Imodel;Gmodel];
            model(isnan(model)) = 0;
        end
    end
    % Set NaNs to zero
     
    % Make the data a column vector and add zeros on the end for the EER constraint and the Tikhonov regularisation
    data = double([data(:);zeros(3*ns,1,'single')]);
    %data = [data(:);zeros(2*ns,1)];
    % Do the inversion
    if background == 1;
        [ava_tmp,lsqflag,~,fiternotmp] = lsqr(IGiter,data,tol,iter,[],[],model);
        ava(:,kk) = single(ava_tmp);
        % try a different solving method?
        ava(:,kk) = bsxfun(@times,cjava,taperrev);
    else     
        %[ava(:,kk),~] = lsqr(IGiter,data,tol,iter,[],[],cjmodel);
        %[ava_tmp,lsqflag,relres,fiternotmp,residvec] = lsqr(IGiter,data,tol,iter,[],[]);
        if l1 == 1 
            
            opts = spgSetParms('optTol',tolerance, 'iterations', max_iterations);
            %opts = spgSetParms('verbosity',0);
            
            if(Wdomain ==1)
                padded_len = 2^(nextpow2(ns)); % Hard coded
                dim_out = [padded_len, 1];
            
                % Set up the spot operators
                P = opKron(eye(2), opPad2D([ns,1], dim_out, 0));
                IGSpot = opMatrix(IGiter);
                W = opKron(eye(2), opSplineWavelet(padded_len, 1, padded_len, 3, 5, '*bspline'));
            
       
                dirac = zeros(4096,1);
                dirac(1)=1;
                wavelet_test = W'*dirac;
                plot(wavelet_test);xlim([1028,3028]);
                title('Symmetric Spline Wavelet');
                set(gca,'xtick',[],'ytick',[])
                
                %preconditioning: Scale the data and operator
                if pc
                input = zeros(ns,1);
                input(50:50:end) = 1;
                t_in = [50:50:size(input,1)];
            
                response = IGSpot'*IGSpot*[input;input];
                
                % peaks
                mag_resp = response(t_in);

                scale1 = interp1(t_in, mag_resp, [1:ns]);
                scale1(isnan(scale1)) = scale1(51);
                scale1 = scale1';
                scale1(900:end) = scale1(900);

                scale = repmat(scale1,32,1) + 1;

                precond = spdiags(1 ./ scale,0, size(scale,1), ...
                    size(scale,1));
                end
                
                A = IGSpot * P'*W';
                
                % Adjoint test
                x = randn(size(W',2),1);
                y = randn(size(IGSpot, 1),1);
                assert((y'*(A*x))/((A'*y)'*x) ~=1)
            else
                A = IGiter;
                 %preconditioning: Scale the data and operator
                 if pc
                input = zeros(ns,1);
                input(50:50:end) = 1;
                t_in = [50:50:size(input,1)];

                response = A'*A*[input;input];
                % peaks
                mag_resp = response(t_in);

                scale1 = interp1(t_in, mag_resp, [1:ns]);
                scale1(isnan(scale1)) = scale1(51);
                scale1 = scale1';
                scale1(900:end) = scale1(900);

                scale = repmat(scale1,32,1) + 1;

                precond = spdiags(1 ./ scale,0, size(scale,1), ...
                    size(scale,1));
                 end
            end
            
           
            if pc
            op_precond = opMatrix(precond);
            end
            
            clear scale response mag_resp
           
            
            
        
                if pc
                [ava_tmp,r,g,info] = spgl1(op_precond*A, op_precond*data, 0, 1e-3, [], opts );
                else
                [ava_tmp,r,g,info] = spgl1(A, data, 0, 1e-3, [], opts );
                end
         
            if plot_convergence==1 & kk==2
                iterations = [1:info.iter];
                figure
                subplot(211)
                [ax,p1,p2] = plotyy(iterations,...
                    info.rNorm2,iterations,info.xNorm1);
                title('SPGL1 Summary');
                ylabel(ax(1), 'Residual');ylabel(ax(2), '|x|');
                xlabel('interation')
                subplot(212)
                plot(info.xNorm1, info.rNorm2);xlabel('|x|');
                ylabel('Residual');title('Pareto Curve');
                
                saveas(gcf,strcat(prefix,'l1_convergence.eps'),'epsc')
            end
                
            if Wdomain ==1
                ava_tmp = P'*W'* ava_tmp;
            end
            
            
            fiternotmp = 0;
            %[ava_tmp, r, info] = lb(IGiter,data,[],[],[],[]);
            if kk==5
                figure;
                subplot(121)
                plot(ava_tmp(1:end/2),[1:size(ava_tmp)/2]);title('Intercept');
                xlabel('Amplitude'); ylabel('sample');
                axis ij
                ylim([1,ns]);xlim([-4,4]);
                subplot(122)
                plot(ava_tmp(end/2 +1:end),[1:size(ava_tmp)/2]);title('Gradient');
                xlabel('Amplitude'); ylabel('sample');
                axis ij
                ylim([1,ns]);xlim([-4,4]);
                
                saveas(gcf,strcat(prefix,'intercept_gradient.eps'),'epsc')
                
                figure;
                intercept_spec = 20*log10(abs(fft(ava_tmp(1:end/2))));
                gradient_spec = 20*log10(abs(fft(ava_tmp(end/2+1:end))));
                
                subplot(211)
                plot(linspace(0,pi,ns/2),intercept_spec(1:end/2));
                title('Intercept');ylabel('Amplitude [db ref arb]');
                xlabel('w [rad]');
                subplot(212)
                plot(linspace(0,pi,ns/2),gradient_spec(1:end/2));
                title('Gradient');ylabel('Amplitude [db ref arb]');
                xlabel('w [rad]');
                
                saveas(gcf,strcat(prefix,'intercept_gradient_spec.eps'),'epsc')
                
            end
                
                
        else
            [ava_tmp,lsqflag,~,fiternotmp,res] = lsqr(IGiter,data,tol,iter,[],[]);
                        if kk==5
                figure;
                subplot(121)
                plot(ava_tmp(1:end/2),[1:size(ava_tmp)/2]);title('Intercept');
                xlabel('Amplitude'); ylabel('sample');
                axis ij
                ylim([1,ns]);xlim([-4,4])
                subplot(122)
                plot(ava_tmp(end/2 +1:end),[1:size(ava_tmp)/2]);title('Gradient');
                xlabel('Amplitude'); ylabel('sample');
                axis ij
                ylim([1,ns]);xlim([-4,4])
                
                saveas(gcf,strcat(prefix,'intercept_gradient.eps'),'epsc')
                
                figure;
                intercept_spec = 20*log10(abs(fft(ava_tmp(1:end/2))));
                gradient_spec = 20*log10(abs(fft(ava_tmp(end/2+1:end))));
                
                subplot(211)
                plot(linspace(0,pi,ns/2),intercept_spec(1:end/2));
                title('Intercept');ylabel('Amplitude [db ref arb]');
                xlabel('w [rad]');
                subplot(212)
                plot(linspace(0,pi,ns/2),gradient_spec(1:end/2));
                title('Gradient');ylabel('Amplitude [db ref arb]');
                xlabel('w [rad]');
                
                saveas(gcf,strcat(prefix,'intercept_gradient_spec.eps'),'epsc')
                
                figure;
                plot(res);
                title('LSQR Summary');xlabel('iteration');
                ylabel('misfit');
                saveas(gcf,strcat(prefix,'convergence.eps'),'epsc')
                
                
            end
        end
        
        if(kk==5)
        d_real = reshape(data, ns, 32);
    
        d_sim  = IGiter*ava_tmp;
        d_sim = reshape(d_sim, ns, 32);
        
        
        wigglePlot(d_real(300:500,:), strcat(prefix,'real_data'));
        wigglePlot(d_sim(300:500,:), strcat(prefix,'modelled_data'));
        wigglePlot(d_real(300:500,:) - d_sim(300:500,:), ...
                   strcat(prefix,'misfit'));
        
        wigglePlot(d_real(1300:1600,:), strcat(prefix,'real_data_end'));
        wigglePlot(d_sim(1300:1600,:), strcat(prefix,'modelled_data_end'));
        wigglePlot(d_real(1300:1600,:) - d_sim(1300:1600,:), strcat(prefix,'misfit'));
     
        
        end
        %residcount(1:length(residvec),kk) = residvec;
        %apply band pass filter to the output
        ava_zeros = ava_tmp == 0;
        % need to apply tapering to the data before the bandpass filter at
        % the start of each cube
        ava_tmp = ava_tmp.*outtaper;
        ava_tmp = bandpass_filter(ava_tmp,(job_meta.s_rate/1000000),(relfreq/3),(relfreq+1),top3dpt,topfreq);
        ava(:,kk) = ava_tmp .* (1.-ava_zeros);
        ava(1,kk) = single(fiternotmp);
        % left out the taper reverse to see if that just kept high
        % amplitudes down
        %ava(:,kk) = bsxfun(@times,ava_tmp,taperrev);
    end
    %[ava(:,kk),~] = lsqr(IGiter,data,1e-2,100,[],[],model);
    
    % mute out any hard zeros from the data
    %ava([fold;fold]==0,kk)=0;
    
%     % apply band pass filter to the output 
%     ava_zeros = ava(:,kk) == 0;
%     ava(:,kk) = bandpass_filter(ava(:,kk),(job_meta.s_rate/1000000),0,relfreq,top3dpt,topfreq);
%     ava(:,kk) = ava(:,kk) .* (1.-ava_zeros);
    
    % Estimate the R^2 confidence in the result. This is the variance ratio:
    % 1-[(sum(data-Gm)^2)/(sum(data-mean)^2)], where the inversion was
    % solving data = Gm.
    if needconf == 1;
        data = reshape(data(1:ns*totalvol,:),[],totalvol);
        digi_confidence(:,kk) = 1-(sum((data-reshape(IGiter(1:totalvol*ns,:)*ava(:,kk),[],totalvol)).^2,2)./sum(bsxfun(@minus,data,sum(data,2)./fold).^2,2));
    end
    % Give a status report, add this to a log file
    reportinc = round((last_iter-first_iter+1)/noofreports);
    curtr = (kk-first_iter+1)/reportinc;
    curtrresid = curtr - round(curtr); 
    if kk == first_iter;
        fprintf('Completed 0 percent: trace %d of %d lsqflag was: %d it took: %d iterations\n',kk-first_iter+1,last_iter-first_iter+1,spgl,fiternotmp)
    elseif (curtrresid == 0) 
        fprintf('Completed %5.2f percent: trace %d of %d lsqflag was: %d it took: %d iterations\n',((100/noofreports)*curtr),kk-first_iter+1,last_iter-first_iter+1,spgl,fiternotmp)
    end
end    
end
toc

%% images
figure
imagesc(ava(1:end/2,:));title('Intercept');ylabel('samples');xlabel('traces');
caxis([-2,2]);colormap(gray);
saveas(gcf,strcat(prefix,'intercept_section.eps'),'epsc')
figure
imagesc(ava(end/2+1:end,:));title('Gradient');ylabel('samples');xlabel('traces');
caxis([-2,2]);colormap(gray);
saveas(gcf,strcat(prefix,'gradient_section.eps'),'epsc')
%% wiggliorams
wigglePlot(ava(2:end/2,:), [prefix, 'intercept_wiggler']);
wigglePlot(ava(end/2+1:end,:), [prefix, 'gradient_wiggler']);
% just the interesting part
wigglePlot(ava(300:500,:), [prefix, 'intercept_wiggler_sub']);
wigglePlot(ava(end/2+300:end/2 + 500,:), [prefix, 'gradient_wiggler_sub']);



%fprintf('Completed %8d traces: lsqflag was: %d last iteration no was: %d \n',(kk,lsqflag,fiternotmp);
clear IGiter IGmatrix;      % remove the preconditioning matrices

%%


if needconf == 1;
    digi_confidence(digi_confidence<0)=0;
    digi_confidence(digi_confidence>1)=1;
    digi_confidence(isnan(digi_confidence))=0;
    digi_confidence = [digi_confidence;zeros(job_meta.n_samples{vol_index_wb}-ns,ntraces)];
end
if output_std == 1;
    std_intercept = [Imodel;zeros(job_meta.n_samples{vol_index_wb}-ns,ntraces)];
    std_gradient = [Gmodel;zeros(job_meta.n_samples{vol_index_wb}-ns,ntraces)];
    %calculate the mimum energy based on chi angle provided
    std_minimum_energy_eer_projection = [bsxfun(@times,Imodel,cosd(chi))+bsxfun(@times,Gmodel,sind(chi));zeros(job_meta.n_samples{vol_index_wb}-ns,ntraces)];
end



resultno = 1;
% Save outputs into correct structure to be written to SEGY.
results_out{resultno,1} = 'Meta data for output files';
results_out{resultno,2}{1,1} = ilxl_read{vol_index_wb}(1:ntraces,:);
%results_out{resultno,2}{2,1} = uint32(zeros(size(traces{vol_index_wb},2),1));
results_out{resultno,2}{2,1} = uint32(zeros(ntraces,1));

ebcstrtowrite = sprintf('%-3200.3200s',[results_out{resultno,1} '  ' ebdichdr '  ' tmpebc]);
results_out{resultno,1} = ebcstrtowrite;

resultno = resultno + 1;
clear ilxl_read;

results_out{1,3} = 'is_gather'; % 1 is yes, 0 is no

results_out{resultno,1} = strcat('digi_intercept',testdiscpt);
%results_out{2,2} = digi_intercept;
results_out{resultno,2} = zeros(job_meta.n_samples{vol_index_wb},ntraces);
% Unflatten data using the window index
%results_out{resultno,2}(win_ind(:,1:ntraces)) = 1000.*ava(1:ns,:);
results_out{resultno,2} = 1000.*ava(1:ns,:);
results_out{resultno,3} = 0;
resultno = resultno + 1;

results_out{resultno,1} = strcat('digi_gradient',testdiscpt);
%results_out{3,2} = digi_gradient;
results_out{resultno,2} = zeros(job_meta.n_samples{vol_index_wb},ntraces);
% Unflatten data using the window index
%results_out{resultno,2}(win_ind(:,1:ntraces)) = 1000.*ava(1+ns:end,:);
results_out{resultno,2} = 1000.*ava(1+ns:end,:);
results_out{resultno,3} = 0;
resultno = resultno + 1;

results_out{resultno,1} = strcat('digi_minimum_energy_eer_projection',testdiscpt); 
%results_out{4,2} = digi_minimum_energy_eer_projection;
digi_minimum_energy_eer_projection = [bsxfun(@times,ava(1:ns,:),cosd(chi))+bsxfun(@times,ava(1+ns:end,:),sind(chi));zeros(job_meta.n_samples{vol_index_wb}-ns,ntraces)];
results_out{resultno,2} = zeros(job_meta.n_samples{vol_index_wb},ntraces);
% Unflatten data using the window index
%results_out{resultno,2}(win_ind(:,1:ntraces)) = 1000.*digi_minimum_energy_eer_projection(1:ns,:);
results_out{resultno,2} = 1000.*digi_minimum_energy_eer_projection(1:ns,:);
results_out{resultno,3} = 0;
resultno = resultno + 1;

if needconf == 1;
    results_out{resultno,1} = strcat('digi_confidence',testdiscpt);
    %results_out{5,2} = digi_confidence;
    digi_confidence = [digi_confidence;zeros(job_meta.n_samples{vol_index_wb}-ns,ntraces)];
    results_out{resultno,2} = zeros(job_meta.n_samples{vol_index_wb},ntraces);
    % Unflatten data using the window index
    %results_out{resultno,2}(win_ind(:,1:ntraces)) = 1000.*digi_confidence(1:ns,:);
    results_out{resultno,2} = 1000.*digi_confidence(1:ns,:);
    results_out{resultno,3} = 0;
    resultno = resultno + 1;
end

% output standard intercept and gradient result
if output_std == 1;
    results_out{resultno,1} = strcat('std_intercept',testdiscpt);
    results_out{resultno,2} = 1000.*std_intercept;
    results_out{resultno,3} = 0;
    resultno = resultno + 1;
    
    results_out{resultno,1} = strcat('std_gradient',testdiscpt);
    results_out{resultno,2} = 1000.*std_gradient;
    results_out{resultno,3} = 0;
    resultno = resultno + 1;

    results_out{resultno,1} = strcat('std_minimum_energy_eer_projection',testdiscpt);
    results_out{resultno,2} = 1000.*std_minimum_energy_eer_projection;
    results_out{resultno,3} = 0;
    resultno = resultno + 1;
end
%%
% segy write function
job_meta.output_dir = './';
if exist(strcat(job_meta.output_dir,'digi_results/'),'dir') == 0
    output_dir = strcat(job_meta.output_dir,'digi_results/');
    mkdir(output_dir);    
else
    output_dir = strcat(job_meta.output_dir,'digi_results/');
end

i_block = str2double(i_block);
node_segy_write(results_out,i_block,job_meta.s_rate/1000,output_dir);

end

%%
function [IGmatrix] = build_operator(totalvol,input_angles,ns_wavelet,wavelet_tmp,ns,hns_wavelet,angle_index,chi,alpha,IGblank,smooth)
    % Build operator for inversion
    for ii = 1:totalvol
        Iwavelet_interp(:,1+(ii-1)*ns_wavelet:ii*ns_wavelet) = wavelet_tmp{ii}';
        Gwavelet_interp(:,1+(ii-1)*ns_wavelet:ii*ns_wavelet) = wavelet_tmp{ii}'*(sind(input_angles(ii)).*sind(input_angles(ii)));
    end

    IGdiagnals = sort(reshape([(-hns_wavelet:hns_wavelet)',bsxfun(@plus,(-hns_wavelet:hns_wavelet)',(-ns:-ns:-ns*(totalvol-1)))],1,[]),'descend');

    Imatrix = spdiags(Iwavelet_interp,IGdiagnals,ns*totalvol,ns);
    Gmatrix = spdiags(Gwavelet_interp,IGdiagnals,ns*totalvol,ns);

    EERmatrix = alpha*[bsxfun(@times,Imatrix(1+ns*(angle_index-1):+ns*angle_index,:),cosd(chi)),bsxfun(@times,Imatrix(1+ns*(angle_index-1):+ns*angle_index,:),sind(chi))];

    IGmatrix = [[Imatrix,Gmatrix;EERmatrix].*IGblank;smooth];
end

%%
function wavelet_norm = wavelet_rms_norm(totalvol,wavelet_interp,vol_index)
    % Normalise the wavelets to have constant energy w.r.t. angle. The energy
    % is set to that of the nearest angle wavelets. Wavelet energy still varies
    % w.r.t. time.
    norm_to = sum(abs(wavelet_interp{vol_index}));
    for ii=1:totalvol
        curwav = sum(abs(wavelet_interp{ii}));
        ratio = norm_to./curwav;
        wavelet_norm{ii} = bsxfun(@times,wavelet_interp{ii},ratio);
    end
end
