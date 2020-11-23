function peak_combos = makePeakCombos(isdynamic)
% makePeakCombos forms the indicator matrix of maximum available peak 
% combos based on an indicator vector of which peaks are dynamic.
%
% INPUTS:
%   isdynamic -- an indicator vector (1 x num_peaks) of which peaks are
%                dynamic. Required.
%
% OUTPUTS:
%   peak_combos -- indicator matrix (num_combos x num_peaks) of which peaks
%                  are on in each combo.
%
% Modified: 20200707 -- Commented - Patrick Stokes
% Created by Michael Prerau
%

if nargin < 1
    error('Vector indicating dynamic peaks isdynamic is required.');
end

num_peaks = length(isdynamic);
% 2^num_peaks is maximum possible number of combos if all are dynamic
D = (0:2^num_peaks-1)';

% pow2(-(num_peaks-1):0) gives zero to one values for each peak
% D*pow2 makes a (num_combos x num_peaks) matrix of values from 0 to 2^(num_combos-1)
% floor reduces values to integers from 0 to 2^(num_combos-1)
% rem reduces values to binary
% fliplr undoes the order reversal from pow2
peak_combos = fliplr(rem(floor(D*pow2(-(num_peaks-1):0)),2));

% Set non-dynamic peaks to be always on
peak_combos(:,~isdynamic) = 1;
% Remove redundant combos
peak_combos = unique(peak_combos,'rows');
