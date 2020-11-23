function plotMPSpectWithSlices(stimes,sfreqs,spect_true,true_alpha,ss_mdl,ss_est,f_plot_resid,ofile_pref)
% plotMPSpectWithSlices is a function to plot an observed spectrogram and 
% its corresponding estimated spectrogram. It optionally plots the 
% estimated residuals or forms an interactive subplot showing the spectra, 
% observation estimate, and individual peak components at time slices. 
% Estimated On-peaks (with true On-peaks, if available) can also be plotted.
%
% Example call for one of the filter estimates in compareEKFs_jobs_wRandSSMs:
%   plotMPSpectWithSlices(stimes,sfreqs,sim_spect,alpha,s,est,2);
% Example call from ex_real_data: 
%   plotMPSpectWithSlices(stimes,sfreqs_trunc,nanpow2db(spect_trunc),[],ss,ss_est,2);
%
% INPUTS:
%   stimes       -- 1 x NT vector of times. Required.
%   sfreqs       -- 1 x NF vector of frequency bins. Required.
%   spect_true   -- NF x NT spectrogram matrix. Required.
%   true_alpha   -- num_combos x NT+1 matrix indicating true On/Off-peak
%                   combinations. Optional, only possible for simulated data.
%   ss_mdl       -- StateSpaceMultiPeak object with parameters used to 
%                   obtain filter estimates in ss_est. The model is needed 
%                   to recompute the filter observation estimate from the 
%                   state estimates in the case when filter what estimated
%                   using a subsampling of frequency bins. Only required 
%                   if estimated frequency bins differ from sfreqs, or if
%                   plotting the estimated On-peaks is desired.
%   ss_est       -- struct containing filter estimates. Required.
%   f_plot_resid -- integer flag indicating arrangement of output figures.
%                   0 - Plots main figure with observed spectrogram in
%                       first axis; estimated spectrogram in second axis;  
%                       an interactive time slice of spectrum, estimate,
%                       and individual peak components in third axis;
%                       and, if ss_mdl.peak_combos is provided, estimated 
%                       On-peaks in fourth axis. 
%                       Does not plot estimated residuals.
%                   1 - Plots estimated residuals as a separate figure.
%                   2 - Plots the estimated residuals in the third axis of 
%                       the main figure in place of the interactive slice
%                       plot.
%                   Default is 1.
%   ofile_pref   -- path and filename prefix for saving figures. 
%                   Figures not saved if empty. Default is empty ''.
% 
% OUTPUTS:
%   Figures generated. Figures saved if ofile_pref provided.
%   
%   20200224 -- (v3) changed to input model (peak_combos and names no
%               longer input separately) so that the observation estimate 
%               can be resampled at different frequency bins. Inputs and
%               their order consequently were changed
%   20200131 -- (v2) added true combos.
%               added option to plot residuals in place of slices. 
% Created by Patrick Stokes
%

%******************************************
% Handle variable inputs and set defaults *
%******************************************
if nargin<8
    ofile_pref = '';
end
if nargin<7
    f_plot_resid = [];
end
if nargin<5
    ss_mdl = [];
end
if nargin<4
    true_alpha = [];
end
if isempty(f_plot_resid)
    f_plot_resid = 1;
end

%**********************************************************************
% Extract necessary variables from model class and estimate structure * 
%**********************************************************************
% If model is provided, determine On-peaks and recompute estimated 
% spectrogram and its components if model frequency bins do not match 
% the frequency bins of the observation 
if ~isempty(ss_mdl)
    peak_combos = ss_mdl.peakCombos;
    peak_names = {ss_mdl.peakModels.peakName};
    if length(sfreqs)~=length(ss_mdl.omega) || any(sfreqs~=ss_mdl.omega)
        yf_hat = zeros(size(spect_true));
        comps_hat_f = zeros([size(ss_est.comps_hat_f,1) size(spect_true)]);
        for tt = 1:size(ss_est.xf_hat,2)
            pick_peaks_on = 1*any(peak_combos(logical(ceil(ss_est.alpha(:,tt))),:),1);    
            [yf_hat(:,tt), comps_hat_f(:,:,tt)] = ss_mdl.getH(ss_est.xf_hat(:,tt),sfreqs, pick_peaks_on);
        end
    else
        yf_hat = [];
        comps_hat_f = [];        
    end
else
    peak_combos = [];
    peak_names = [];  
end
% Extract estimated spectrogram and its components, if not already done
if isempty(yf_hat)
    yf_hat = ss_est.yf_hat;
end
if isempty(comps_hat_f)
    comps_hat_f = ss_est.comps_hat_f;
end

%*********************
% Set up main figure *
%*********************
f1 = figure;
% If peak_combos are available, On-peaks are plotted in fourth axis
if ~isempty(peak_combos)
    ax = figdesign(4,1,'orient','landscape','type','usletter','margins',[.05 .05 .05 .05 .05 .05]);
else
    ax = figdesign(3,1,'orient','landscape','type','usletter','margins',[.05 .05 .05 .05 .05 .05]);
end
if f_plot_resid==2 
    % Residuals plotted in third axis
    linkaxes(ax(1:end),'x');
    linkaxes(ax(1:3),'y');
else
    % Third axis is interactive time slices
    if ~isempty(peak_combos)
        linkaxes(ax([1:2 4]),'x');
        linkaxes(ax(1:2),'y');
    else
        linkaxes(ax(1:2),'xy');
    end
end

%******************************************
% Plot observed spectrogram in first axis *
%******************************************
axes(ax(1));
imagesc(stimes, sfreqs, spect_true);
axis xy;
colormap(jet(1024));
title('Spectrogram');
climscale
clim = caxis(ax(1));
ylabel('Frequency (Hz)');

%********************************************
% Plot estimated spectrogram in second axis *
%********************************************
axes(ax(2));
imagesc(stimes, sfreqs, yf_hat(:,2:end));
axis xy;
colormap(jet(1024));
title('Filter Estimate');
caxis(clim);
ylabel('Frequency (Hz)');
ylim([0 40]);
axis tight

%************************************************************
% Plot true estimated On-peaks in fourth axis, if available *
%************************************************************
if ~isempty(peak_combos)
    num_peaks = size(comps_hat_f,1);
    peaks_on = zeros(num_peaks,size(ss_est.alpha,2));
    for ii = 1:size(peaks_on,2)
        % If have probabilities over combos (e.g. initial state),
        % take all peaks in any positive probability combo to be on.
        peaks_on(:,ii) = 1*any(peak_combos(logical(ceil(ss_est.alpha(:,ii))),:),1);
    end
    % peaks_on = peaks_on .* repmat((num_peaks:-1:1)',[1 size(peaks_on,2)]);
    peaks_on = peaks_on .* repmat((1:num_peaks)',[1 size(peaks_on,2)]);
    peaks_on(peaks_on == 0) = nan;
    axes(ax(4));
    plot(stimes,peaks_on(:,2:end),'Linewidth',10);
    hold on;
    if ~isempty(true_alpha)
        true_peaks = zeros(num_peaks,size(true_alpha,2));
        for ii = 1:size(true_peaks,2)
            % If have probabilities over combos (e.g. initial state),
            % take all peaks in any positive probability combo to be on.
            true_peaks(:,ii) = 1*any(peak_combos(logical(ceil(true_alpha(:,ii))),:),1);
        end
        % true_peaks = true_peaks .* repmat((num_peaks:-1:1)',[1 size(true_peaks,2)]);
        true_peaks = true_peaks .* repmat((1:num_peaks)',[1 size(true_peaks,2)]);
        true_peaks(true_peaks == 0) = nan;
        plot(stimes,true_peaks(:,2:end),'k','Linewidth',2);
        hold on;        
    end
    ylim([0.5 num_peaks+0.5]);
    xlabel('Time (sec)');
    yticks(1:num_peaks);
    if isempty(peak_names)
        ylabel('Peak Number'); 
    else
        % yticklabels(fliplr(peak_names));
        yticklabels(peak_names);
    end
    title('On-Peaks');
end
    
%****************************************
% Plot interactive slices in third axis *
%****************************************
if f_plot_resid~=2
    N = length(stimes);
    num_freqs = length(sfreqs);
    num_peaks = size(comps_hat_f,1);
    aug_comps = zeros(num_freqs,N,num_peaks+2);
    aug_comps(:,:,1) = spect_true;
    aug_comps(:,:,2) = yf_hat(:,2:end);
    for ii = 1:num_peaks
        aug_comps(:,:,ii+2) = squeeze(comps_hat_f(ii,:,2:end));
    end
    
    drawnow;
    axis(ax(3));
    set(ax(3),'colororder',[0 0 0; 1 0 0;get(ax(3),'colororder')]);
    
    slicepopup(gcf, ax, stimes, sfreqs, aug_comps, 'Time','Frequency (Hz)','Power (db)', 'y', 1, ax(3));
end

%**************************************************************
% Plot estimated residuals either in third axis or new figure *
%**************************************************************
if f_plot_resid>0
    if f_plot_resid==1
        f2 = figure;
        ca = gca;
    else
        axes(ax(3));
        ca = gca; % ax(3);
    end
    imagesc(ca,stimes, sfreqs, spect_true-yf_hat(:,2:end));
    axes(ca);
    axis xy;
    colormap(jet(1024));
    climscale;
    ylabel('Frequency (Hz)');
    title('Filter Residual');
end

%*******************
% Save main figure *
%*******************
if ~isempty(ofile_pref)
   print(f1,[ofile_pref '_estSpect'],'-depsc');
end
 
end

