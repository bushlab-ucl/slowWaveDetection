
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IEDs detector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Detection of IEDs based on criteria defined in Gelinas et al. (2016)
% Nature Medicine.
% If used, thank you for citing:
% Sheybani et al. (2023) Nat Comm (in revision)
% Given its substantial theorical inputs, it would be fair to cite Gelinas
% et al. (2016) Nature Medicine as well

% The code is shared without any warranty

% Laurent Sheybani, laboratory of Prof. Matthew C Walker, University
% College London (UCL), London, UK

% The first step detects potential IEDs based on Gelinas et al. (2016). The
% second step asks the reviewer to refine the detection. 


% Data inputs
MAT                             = []     % MAT should be a matrix with rows as electrodes and columns as timeframes
fc_low                          = 20;    % Gelinas et al.: 60
fc_high                         = 80;    % same as in Gelinas et al.
refractory_ms                   = 1000;  % after each variable, all variables following in the 'refractory_ms' period will be deleted - not described in Gelinas et al. In practice, the visual confirmation overwrites this.
window_peak_ms                  = 200;   % window in ms around the marker, that identifies the IED peak
volume_conduction_window_ms     = 300;   % not described in Gelinas et al.
threshold_filtered_unit         = 3.5;   % 3 in Gelinas et al.
threshold_unfiltered_unit       = 4;     % 3 in Gelinas et al.
sf                              = ;      % sampling frequency of the file
filter_order                    = ;      % order of the filter
name_file                       = '';    % name of the file being process

clear idx idx2 IED_position

% Actual script starts here
[b_low a_low] = butter(filter_order, 2*fc_low/sf,'high'); % coefficients of the high pass-filter
[b_high a_high] = butter(filter_order, 2*fc_high/sf,'low'); % coefficients of the low pass-filter

% Filter data
thedata = filtfilt(b_low, a_low, MAT'); % filtering (high-pass)
thedata = filtfilt(b_high, a_high, thedata)'; % filtering (low-pass)
refractory = round(refractory_ms / 1000 * sf);
window_peak = [-round((window_peak_ms/2) / 1000 * sf) : round((window_peak_ms/2) / 1000 * sf)];
volume_conduction_window = ones(1,round(volume_conduction_window_ms / 1000 * sf));

% Select event with filtered envelope > 3 times the baseline
envl = abs(hilbert(thedata'))';
baseline_envelop = mean(envl,2);
threshold = baseline_envelop .* threshold_filtered_unit;
for j = 1 : size(envl,1)
    idx{j} = find(envl(j,:) > threshold(j));
end

% Remove events with unfiltered envelope < 3 times the baseline
unfiltered_envl = abs(hilbert(MAT'))';
baseline_unfiltered_envl = mean(unfiltered_envl,2);
threshold_unfiltered = baseline_unfiltered_envl .* threshold_unfiltered_unit;
for j = 1 : length(idx)
    idx{j}(unfiltered_envl(j,idx{j}) < threshold_unfiltered(j)) = [];
end

% Remove events if they occur during the refractory period
idx2 = idx;
for k = 1 : length(idx2)
    j = 1;
    while j < length(idx2{k})
        a = [zeros(j,1) ; ismember(idx2{k}(j+1 : end),idx2{k}(j)+1 : idx2{k}(j) + refractory)'];
        idx2{k}(a==1) = [];
        j = j + 1;
    end
end

% Select the peak of the IED and tag it whether it is positive or negative
for j = 1 : size(MAT,1)
    m = 1;
    o = 1;
    for k = 1 : length(idx2{j})
        if idx2{j}(k) + window_peak(1) > 0 && idx2{j}(k) + window_peak(end) < size(MAT,2)
            [maxi peak_IED_temp] = max(unfiltered_envl(j,idx2{j}(k) + window_peak));

            % Save IED depending on their polarity
            if MAT(j,peak_IED_temp + idx2{j}(k) + window_peak(1) - 1) > 0
                tf_of_POS_peak_IED{j}(m) = peak_IED_temp + idx2{j}(k) + window_peak(1) - 1;
                m = m + 1;
            else
                tf_of_NEG_peak_IED{j}(o) = peak_IED_temp + idx2{j}(k) + window_peak(1) - 1;
                o = o + 1;
            end
        end
    end
end

% If > 1 IED is detected within 'volume_conduction_window' ms, keep the
% spike with the highest amplitude
IED_position = zeros(size(MAT,1), size(MAT,2));
for j = 1 : size(MAT,1)
    if exist('tf_of_NEG_peak_IED')
        IED_position(j,tf_of_NEG_peak_IED{j}) = 1;
    end
    if exist('tf_of_POS_peak_IED')
        IED_position(j,tf_of_POS_peak_IED{j}) = 1;
    end
end
overlapping_IED = conv(sum(IED_position,1),volume_conduction_window,'same');
idx_overlapping_IED = zeros(1,length(overlapping_IED));
idx_overlapping_IED(overlapping_IED > 1) = 1;

% Each time a '1' appears in start_overlapping_IED, then > 1 electrode
% record an IED
start_overlapping_IED = diff(idx_overlapping_IED);
idx_vol_cond_start = find(start_overlapping_IED == 1);
idx_vol_cond_end = find(start_overlapping_IED == -1);

% The argument '(0.2*sf)' is made to include the first IED in the burst of
% IEDs recorded across several electrodes, as this IED will be counted as
% one
if length(idx_vol_cond_start) < length(idx_vol_cond_end)
    end_loop = length(idx_vol_cond_start);
else
    end_loop  = length(idx_vol_cond_end);
end
for j = 1 : end_loop
    if round(idx_vol_cond_start(j) - (0.2*sf)) < 1
        [max_IED maxIED_idx] = max(unfiltered_envl(:,[1 : idx_vol_cond_end(j)])');
    else
        [max_IED maxIED_idx] = max(unfiltered_envl(:,[round(idx_vol_cond_start(j) - (0.2*sf)) : idx_vol_cond_end(j)])');
    end
    [val max_IED_electrode] = max(max_IED);
    for k = 1 : size(MAT,1)
        if k ~= max_IED_electrode
            if ~isempty(tf_of_POS_peak_IED{k})
                tf_of_POS_peak_IED{k}(ismember(tf_of_POS_peak_IED{k}, round(idx_vol_cond_start(j) - (0.2*sf)) : idx_vol_cond_end(j))) = [];
            end
            if ~isempty(tf_of_NEG_peak_IED{k})
                tf_of_NEG_peak_IED{k}(ismember(tf_of_NEG_peak_IED{k}, round(idx_vol_cond_start(j) - (0.2*sf)) : idx_vol_cond_end(j))) = [];
            end
        end
    end
end

if exist('tf_of_POS_peak_IED')
    data.POS_IED = tf_of_POS_peak_IED;
else
    data.POS_IED = [];
end
if exist('tf_of_NEG_peak_IED')
    data.NEG_IED = tf_of_NEG_peak_IED;
else
    data.NEG_IED = [];
end

data.eeg = eegData.eeg;
data.time = eegData.time;
data.sf = sf;
data.name = name_file;

clear tf_of_NEG_peak_IED tf_of_POS_peak_IED idx idx2 idx_vol_cond_end idx_vol_cond_start start_overlapping_IED idx_overlapping_IED IED_position

data.fc_low = fc_low;
data.fc_high = fc_high;
data.refractory_ms = refractory_ms;
data.window_peak_ms = window_peak_ms;
data.volume_conduction_window_ms = volume_conduction_window_ms;

cd('Your directory for supplementary scripts\');
if isfield(data,'POS_IED')
    if isfield(data,'NEG_IED')
        mrks = [cell2mat(data.POS_IED)' ; cell2mat(data.NEG_IED)'];
    else
        mrks = cell2mat(data.POS_IED);
    end
else
    mrks = cell2mat(data.NEG_IED);
end

% From here on, use data.POS_IED and data.NEG_IED for the visual
% confirmation.

