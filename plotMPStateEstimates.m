function plotMPStateEstimates(ss_est,ss_mdl,f_onoff,x_true,alpha_true,ss_gen,f_add_vol,f_transform,ofile_pref)
% plotMPStateEstimates plots estimated states with uncertainties. 
% Parameters for each peak are plotted in individual figures, with the
% parameters as separate subplots. 
%
% Parameters for each peak are plotted using the internal function 
% plotStateEstimates. Additional internal functions getPeakVolume, 
% combo2PickParams, combo2PickPeaks, and getCI are used to 
% compute/extract necessary data.
%
% INPUTS:
%   ss_est      -- struct containing filter estimates. Required.
%   ss_mdl      -- StateSpaceMultiPeak object with parameters used to 
%                  obtain filter estimates in ss_est. Required.
%   f_onoff     -- binary flag indicating whether estimates and uncertainties 
%                  are plotted for all times or only their peak is on. 
%                  0 - estimates displayed at all times.
%                  1 - estimates are not shown when peak is off.
%   x_true      -- dim_x x NT+1 matrix of true state values
%   alpha_true  -- num_combos x NT+1 matrix of true On/Off-peak
%                  combinations 
%   ss_gen      -- depracated. No longer used.
%   f_add_vol   -- binary flag indicating whether to compute and plot the 
%                  approximate volume of the peak in an additional subplot.
%                  0 - do not compute and plot volume.
%                  1 - compute and plot volume. 
%                  Default is 0.
%   f_transform -- binary flag indicating whether estimates and uncertainties 
%                  are plotted as raw state values or as values transformed 
%                  by the link/bound functions. 
%                  0 - raw values.
%                  1 - transformed values.
%                  Default is 1.
%   ofile_pref  -- path and filename prefix for saving figures. 
%                  Figures not saved if empty. Default is empty ''.
%
% OUTPUTS:
%   Figures generated. Figures saved if ofile_pref provided. 
%
% Modified:
%   20200210 -- Adds option to save by providing an outile prefix
%   20200203 -- This version adds On-peak indicator subplot at top.
%               To do so, it has to add a distinction b/w on/off of the indicator 
%               plot of the peak and on/off of the state plots.
%               This has to be handled in this outer function and the inner
%               plotStateEstimates function.
% Created: Patrick Stokes            
%
%

%*************************
% Handle variable inputs *
%*************************
if nargin < 9
    ofile_pref = '';
end
if nargin < 8
    f_transform = [];
end
if nargin < 7
    f_add_vol = [];
end
if nargin < 6
    ss_gen = [];
end
if nargin < 5
    alpha_true = [];
end
if nargin < 4
    x_true = [];
end
if nargin < 3
    f_onoff = [];
end
if nargin < 2
    error('StateSpaceMultiPeak object, ss_mdl, used to run filter required.');
end
if nargin < 1
    error('Estimate structure, ss_est, containing filter estimates required.');
end
    
%*****************************
% Form some needed variables *
%*****************************
    
% Matrix of color values
% Current initialized to empty. Set in plotStateEstimates as get(gca,'colororder'). 
clr_vals = []; 

% Matrix of On/Off indicators for each parameter at each time
% Rows for parameters from same peak are identical
pick_var_on = combo2PickParams(ss_est.alpha,ss_mdl);
% Set to empty num_params x 0 vector, if true On/Off-combinations not provided 
pick_on_true = combo2PickParams(alpha_true,ss_mdl);

% Vector of peak volumes
if f_add_vol
    [peak_vol,~,~] = getPeakVolume(ss_est.comps_hat_f,ss_mdl.omega);
else
    peak_vol = [];
end

%********************************
% Plot parameters for each peak *
%********************************
for ii = 1:ss_mdl.numPeaks
    % Get indices of parameters for current peak
    idx_params = find(ss_mdl.getPeakIdxs(ii));
    
    % Subselect relevant state and covariance estimates 
    xf_hat = ss_est.xf_hat(idx_params,:);
    Pf = ss_est.Pf(idx_params,idx_params,:);
    if ~isempty(x_true)
        sub_x_true = x_true(idx_params,:);
    else
        sub_x_true = [];
    end
    sub_on_true = pick_on_true(idx_params,:);
    
    % Plot parameters for current peak
    if ~isempty(peak_vol)
        plotStateEstimates(xf_hat,Pf,pick_var_on(idx_params,:),ss_est.t,sub_x_true,sub_on_true,ss_mdl.peakModels(ii),peak_vol(ii,:),[],clr_vals,f_onoff,f_transform,{ss_mdl.peakModels(ii).modelParams.linkFunParams});
    else
        plotStateEstimates(xf_hat,Pf,pick_var_on(idx_params,:),ss_est.t,sub_x_true,sub_on_true,ss_mdl.peakModels(ii),[],[],clr_vals,f_onoff,f_transform,{ss_mdl.peakModels(ii).modelParams.linkFunParams});
    end
    
    % Save file
    if ~isempty(ofile_pref)
        print(gcf,[ofile_pref '_estState_peak' num2str(ii-1)],'-depsc');
    end
end
end

function pick_on_params = combo2PickParams(pick_combo,m)
% Converts a matrix of On/Off-combination indicators to a matrix 
% of On/Off-parameter indicators. Rows (parameters) from the same peak are 
% identical.
%
% INPUTS:
%   pick_combo -- num_combos x num_time matrix of binary On/Off indicators.
%                 Default is empty [].
%   m          -- StateSpaceMultiPeak objects with corresponding combos.
%                 Required if pick_combo provided.
% OUTPUTS:
%   pick_on_params -- num_params x num_time matrix of binary On/Off
%                     indicators. **NOTE: Returns a num_params x 0 empty vector
%                     if pick_combos is empty.**
%

if nargin<2
    m = [];
end
if nargin<1
    pick_combo = [];
end
if isempty(pick_combo)
    pick_on_params = zeros(m.totalPeakParams,0);
else
    if isempty(m)
       error('StateSpaceMultiPeak object, m, is required.'); 
    end
    % Number of time points
    N = size(pick_combo,2);
    pick_on_params = zeros(m.totalPeakParams,N);
    % For each time point, determine which peaks are On, and set their 
    % corresponding parameters to On
    for nn = 1:N
        tmp_combo = pick_combo(:,nn)>0;
        if sum(tmp_combo)>1
            % Use all peaks with positive probability, if more than one
            pick_peaks = sum(m.peakCombos(tmp_combo,:),1)>0;
        else
            pick_peaks = m.peakCombos(tmp_combo,:);
        end
        pick_on_params(:,nn) = m.getPeakIdxs(logical(pick_peaks));
    end
end
end

function [vol,CI_lo,CI_hi] = getPeakVolume(comps,omega,p_lims)
% Computes the approximate volumes of the estimated peak components,
% as well frequency limits of a central interval containing a 
% specified percentile of the peak volume.
%
% INPUTS:
%   comps  -- num_peaks x dim_y x N matrix of state estimates. Required.
%   omega  -- dim_y x 1 vector of frequency bin values. Required.
%   p_lims -- 1 x 2 vector of lower and upper percentile limits of peak
%             volumes. Default is [0.05 0.95].
% OUTPUTS:
%   vol   -- num_peaks x N matrix of estimated peak volumes 
%   CI_lo -- num_peaks x N matrix of lower interval values
%   CI_hi -- num_peaks x N matrix of upper interval values
%

% Variable inputs and defaults
if nargin < 3
    p_lims = [];
end
if isempty(p_lims)
    p_lims = [0.05 0.95];
end

num_peaks = size(comps,1);
N = size(comps,3); % number of time points
omega = omega(:);

% Finds frequency bin widths, replicates and reshapes to match 
% dimensions of comps for element-wise operation
dOMEGA = permute(repmat(omega(2:end)-omega(1:(end-1)),[1 num_peaks N]),[2 1 3]);
% Average bin height
vol = (comps(:,2:end,:)+comps(:,1:(end-1),:))/2;
% Height x width
vol = squeeze(sum(vol.*dOMEGA,2));
        
% For each peak at each time point, compute limits of frequency interval 
% containing specified percentile of volume 
CI_hi = zeros(num_peaks,N);
CI_lo = zeros(num_peaks,N);
for ii = 1:num_peaks
    for jj = 1:N
        % Current cdf of peak heights
        peak_CDF = cumsum(comps(ii,:,jj));
        peak_CDF = peak_CDF/max(peak_CDF);
        % Find current lower limit
        low_idx = find(peak_CDF<=p_lims(1),1,'last');
        if isempty(low_idx)
            low_idx = 1;
        end
        % Find current upper limit
        high_idx = find(peak_CDF>=p_lims(2),1,'first');
        if isempty(high_idx)
            high_idx = length(omega);
        end
        % Store
        CI_hi(ii,jj) = omega(high_idx);
        CI_lo(ii,jj) = omega(low_idx);
    end
end
end

function [CI_lo, CI_hi] = getCI(xhat,P)
% Compute approximate 95%-confidence intervals for state estimates as
% 2 times the square root of the variance from the diagonal of the 
% estimated error covariance.
%
% INPUTS:
%   xhat -- dim_x x NT+1 matrix of state estimates
%   P    -- dim_x x dim_x x NT+1 array of covariance matrices
% OUTPUTS:
%   CI_lo -- dim_x x NT+1 matrix of lower interval values
%   CI_hi -- dim_x x NT+1 matrix of upper interval values
%

CI_lo = zeros(size(xhat));
CI_hi = zeros(size(xhat));

for ii = 1:size(xhat,2)
    CI_lo(:,ii) = xhat(:,ii) - 2*sqrt(diag(P(:,:,ii)));
    CI_hi(:,ii) = xhat(:,ii) + 2*sqrt(diag(P(:,:,ii)));
end

end

function plotStateEstimates(xhat,P,pick_var_on,t,x_true,pick_on_true,peak_mdl,peak_vol,ax,clr_vals,f_onoff,f_transform,x_bounds)
% plotMPStateEstimates plots estimated states with uncertainties. 
% Called in plotMPStateEstimes for each peak separately.  
% Each parameter is ploted in its own subplot. An indicator of when the
% peak is estimated to be On is also plotted.
% True state values and On-indicator values can be added if available. 
% Peak volume can also be plotted is provided.
%
% INPUTS:
%   xhat         -- dim_x x dim_t matrix of state estimates. Required.
%   P            -- dim_x x dim_x x dim_t arrary of state error covariance
%                   matrices. Required.
%   pick_var_on  -- dim_x x dim_t matrix of indicators of when paramters
%                   are estimated to be on. For a single peak, all rows are identical, but 
%                   the function was originally applied to parameters for 
%                   multiple peaks. Default is a matrix of all ones.
%   t            -- 1 x dim_t vector of times. Default t = 1:size(xhat,2).
%   x_true       -- dim_x x dim_t matrix of true state values. Optional.
%                   Default is empty [].
%   pick_on_true -- dim_x x dim_t matrix of indicators of when paramters
%                   are truly on. For a single peak, all rows are identical, but 
%                   the function was originally applied to parameters for 
%                   multiple peaks. Default is is empty [].
%   peak_mdl     -- PeakModel object containing its parameter names and link
%                   functions. Default is empty []. Required if f_transform==1.
%   peak_vol     -- 1 x dim_t vector of peak volume values. Optional.
%                   Default is empty [].
%   ax           -- vector of figure axis handles if figure and axes were
%                   formed outside of the function. 
%                   Default is empty [], in which case a new figure is formed. 
%   clr_vals     -- dim_x x 3 matrix of color values for plotting each
%                   each parameter. Default is 7 x 3 matrix returned by
%                   get(gca,'colororder'), which will need to be augmented 
%                   if a peak is made that has more than 7 parameters.
%   f_onoff      -- binary flag indicating whether estimates and uncertainties 
%                   are plotted for all times or only their peak is on. 
%                   0 - estimates displayed at all times.
%                   1 - estimates are not shown when peak is off.
%                   Default is 0.
%   f_transform  -- binary flag indicating whether estimates and uncertainties 
%                   are plotted as raw state values or as values transformed 
%                   by the link/bound functions. 
%                   0 - raw values.
%                   1 - transformed values.
%                   Default is 1.
%   x_bounds     -- 1 x dim_x cell array of 1 x 2 vectors of bounds for each
%                   parameter used to set ylims in plots. Default is cell
%                   array of empty vectors.
%
% OUTPUTS:
%   Figure generated if ax not provided. Otherwise plotted in axes of ax.  
%
% Modified:
%   20200210 -- Adds option to save by providing an outile prefix
%   20200203 -- This version adds On-peak indicator subplot at top.
%               To do so, it has to add a distinction b/w on/off of the indicator 
%               plot of the peak and on/off of the state plots.
%               This has to be handled in this outer function and the inner
%               plotStateEstimates function.
% Created: -- Patrick Stokes              
%
%

%*************************
% Handle variable inputs *
%*************************
if nargin<13
    x_bounds = [];
end
if nargin<12
    f_transform = [];
end
if nargin<11
    f_onoff = [];
end
if nargin<10
    clr_vals = [];
end
if nargin<9
    ax = [];
end
if nargin<8
    peak_vol = [];
end
if nargin<7
    peak_mdl = [];
end
if nargin<6
    pick_on_true = [];
end
if nargin<5
    x_true = [];
end
if nargin<4
    t = [];
end
if nargin<3
    pick_var_on = [];
end

%***************
% Set defaults *
%***************
dim_x = size(xhat,1);

if isempty(x_bounds)
   x_bounds = cell(1,dim_x); 
end

if isempty(t)
    dim_t = size(xhat,2);
    t = 1:dim_t;
else
    dim_t = length(t);
end

if isempty(pick_var_on)
    pick_var_on = ones(dim_x,dim_t);
end

if isempty(pick_on_true)
    % This is set to all ones, but only used if x_true is provided,
    % i.e., if true states are given but true On/Off-combos unknown,
    %       then the true state is always plotted.
    pick_on_true = ones(dim_x,dim_t);
end

if isempty(peak_vol)
    f_plot_vol = 0;
else
    f_plot_vol = 1;
end

% Generate figure if existing axes not provided
if isempty(ax)
    figure;
    num_axes = max([dim_x 4])+f_plot_vol;
    merges = cell(1,num_axes);
    for pp = 1:(num_axes)
        merges{pp} = (1:3) + 3*(pp-1) + 1;
    end
    ax = figdesign(3*num_axes+1,1,'orient','landscape','type','usletter','merge',merges,'margins',[.05 .05 .05 .05 .05 .05]);
end

if isempty(clr_vals)
    clr_vals = get(ax(1),'colororder');
end

%****************************************************
% Plot indicator of when peak is estimated to be On *
%****************************************************
axes(ax(1));
tmp_pick_var_on = pick_var_on;
tmp_pick_var_on(tmp_pick_var_on==0) = nan;
plot(t,tmp_pick_var_on(1,:),'Color',[.5 .5 .5],'Linewidth',10);

% Add indicator of when peak is truly On if available
if ~isempty(x_true)
    hold on;
    tmp_pick_on_true = pick_on_true;
    tmp_pick_on_true(tmp_pick_on_true==0) = nan;
    plot(t,tmp_pick_on_true(1,:),'k','Linewidth',2);
end

% Remove tick labels and add title
set(ax(1),'xticklabels','','yticklabels','');
if ~isempty(peak_mdl)
    title(ax(1), [peak_mdl.peakName ' - On']);
else
    title(ax(1), ['Peak On']);
end

%****************************************************
% Plot parameter estimates and confidence intervals *
%****************************************************

% Compute confidence intervals
[x_CI_lo,x_CI_hi] = getCI(xhat,P);

% If displaying estimates at all times, set on indicators to all ones.
if ~f_onoff
    pick_var_on = ones(size(pick_var_on));
end

% Plot each parameter in a separate subplot
for pp = 1:dim_x
    % Find chunks of time when peaks are estimated to be On
    % Will be a single chunk if always on or f_onoff==0
    [~, idx_chunks] = consecutive(pick_var_on(pp,:));
    
    axes(ax(pp+1));
    hold on;
    % Plot estimate and shaded confidence interval for each chunk
    for cc = 1:length(idx_chunks)
        idx_chunk = idx_chunks{cc};
        % Transform values with link function if needed 
        if f_transform
            mid = peak_mdl.modelParams(pp).linkFun(xhat(pp,idx_chunk));
            hi = peak_mdl.modelParams(pp).linkFun(x_CI_hi(pp,idx_chunk));
            lo = peak_mdl.modelParams(pp).linkFun(x_CI_lo(pp,idx_chunk));
        else
            mid = xhat(pp,idx_chunk);
            hi = x_CI_hi(pp,idx_chunk);
            lo = x_CI_lo(pp,idx_chunk);
        end
        lh = shadebounds(t(idx_chunk),mid,hi,lo,clr_vals(pp,:),min(clr_vals(pp,:)+.3,1),'none');
        set(lh,'linewidth',1');
    end
    
    % Plot true state values if provided
    if ~isempty(x_true)
        % Plot true values in solid line when peak is on
        [~, idx_chunks] = consecutive(pick_on_true(pp,:));
        for cc = 1:length(idx_chunks)
            idx_chunk = idx_chunks{cc};
            if f_transform
                plot(t(idx_chunk),peak_mdl.modelParams(pp).linkFun(x_true(pp,idx_chunk)),'k','linewidth',2);
            else
                plot(t(idx_chunk),x_true(pp,idx_chunk),'k','linewidth',2);
            end
        end
        % Plot true values in dashed line when peak is off
        [~, idx_not_chunks] = consecutive(abs(pick_on_true(pp,:)-1));
        for cc = 1:length(idx_not_chunks)
            idx_chunk = idx_not_chunks{cc};
            % Extend chunks by +/-1 to connect plotting segments
            if idx_chunk(1)>1
                idx_chunk = [idx_chunk(1)-1 idx_chunk];
            end
            if idx_chunk(end)<size(x_true(pp,:),2)
                idx_chunk = [idx_chunk idx_chunk(end)+1];
            end
            if f_transform
                plot(t(idx_chunk),peak_mdl.modelParams(pp).linkFun(x_true(pp,idx_chunk)),'k--','linewidth',2);
            else
                plot(t(idx_chunk),x_true(pp,idx_chunk),'k--','linewidth',2);
            end
        end
    end
    
    % Set y-axis bounds
    if ~isempty(x_bounds{pp})
        ylim(x_bounds{pp});
    end
    % Clear x-axis ticklabels for all but bottom plot
    if pp < dim_x
        set(ax(pp+1),'xticklabels','');
    end
    % Add title
    if ~isempty(peak_mdl)
        title(ax(pp+1), [peak_mdl.peakName ' - ' peak_mdl.paramNames{pp}]);
    else
        title(ax(pp+1), ['x' num2str(pp)]);
    end    
end

%********************************
% Plot peak volume if available *
%********************************
if f_plot_vol
    axes(ax(dim_x+1));
    plot(t, peak_vol,'b','linewidth',1);
    title('Peak Volume');
    axis tight
end

%******************************
% Set x-axis limits and label *
%******************************
linkaxes(ax,'x');
xlim([min(t) max(t)]);
xlabel('Time (sec)');

end