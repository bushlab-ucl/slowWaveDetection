
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slow wave detector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Detection of slow waves based on criteria defined in Frauscher et al.
% (2015) Brain and Riedner et al. (2007) Sleep.
% If used, thank you for citing:
% Sheybani et al. (2023) Nat Comm (in revision)
% Sheybani et al. (2023) Brain Communications
% Given their substantial theorical inputs, it would be fair to cite
% Frauscher et al. (2015) Brain and Riedner et al. (2007) Sleep as well

% The code is shared without any warranty

% Laurent Sheybani, laboratory of Prof. Matthew C Walker, University
% College London (UCL), London, UK

rawdata                             = ;          % lines are electrodes, columns are timeframes
n                                   = ;          % filter order
sf                                  = ;          % sampling frequency
fc_low                              = 0.5;       % high-pass for SW detection
fc_high                             = 4;         % low-pass for SW detection
min_duration_ZeroCrossing_sec       = 0.25;      % minimal duration of a half-wave, in s
max_duration_ZeroCrossing_sec       = 1;         % maximal duration of a half-wave, in s
keep_above_prctle                   = [];        % threshold above which candidate waves are saved

% For simplicity, we remove the condition "no_IED_before", and assume that
% the post-processing pruning controls for that

%% Filter the data for SW detection
[b_low a_low] = butter(n, 2*fc_low/sf,'high'); % coefficients of the high pass-filter
[b_high a_high] = butter(n, 2*fc_high/sf,'low'); % coefficients of the low pass-filter

min_duration_ZeroCrossing = min_duration_ZeroCrossing_sec * sf;
max_duration_ZeroCrossing = max_duration_ZeroCrossing_sec * sf;

thedata = filtfilt(b_low, a_low, rawdata'); % filtering (high-pass)
thedata = filtfilt(b_high, a_high, thedata)'; % filtering (low-pass)

%% Look for SW of default polarity
idx_ZeroCrossing = iFindZeroCrossing(thedata);
clear PosToNeg_ZeroCrossing
for k = 1 : length(idx_ZeroCrossing)
    if thedata(k,idx_ZeroCrossing{k}(1) + 1) > thedata(k,idx_ZeroCrossing{k}(1))
        PosToNeg_ZeroCrossing{k} = idx_ZeroCrossing{k}(2 : 2 : end);
    else
        PosToNeg_ZeroCrossing{k} = idx_ZeroCrossing{k}(1 : 2 : end);
    end
end

% Initiate variables
onset_SWA = cell(1,length(PosToNeg_ZeroCrossing));
offset_SWA = cell(1,length(PosToNeg_ZeroCrossing));
SWA = cell(1, length(PosToNeg_ZeroCrossing));
SWA_onset = cell(1, length(PosToNeg_ZeroCrossing));
SWA_offset = cell(1, length(PosToNeg_ZeroCrossing));
SWA_middle = cell(1, length(PosToNeg_ZeroCrossing));

% Across electrodes
disp('Negative waves...')
clear diff_ZeroCrossing
for k = 1 : length(PosToNeg_ZeroCrossing)
    if isempty(PosToNeg_ZeroCrossing{k})
        SWA_onset{k} = [];
        SWA_offset{k} = [];
        SWA_middle{k} = [];
        SWA{k} = [];
    else
        diff_ZeroCrossing = diff(idx_ZeroCrossing{k});
        if thedata(k,idx_ZeroCrossing{k}(1) + 1) > thedata(k,idx_ZeroCrossing{k}(1))
            PosToNeg_diff_ZeroCrossing = diff_ZeroCrossing(2 : 2 : end);
        else
            PosToNeg_diff_ZeroCrossing = diff_ZeroCrossing(1 : 2 : end);
        end
        onset_SWA{k} = PosToNeg_ZeroCrossing{k}(PosToNeg_diff_ZeroCrossing > min_duration_ZeroCrossing...
            & PosToNeg_diff_ZeroCrossing < max_duration_ZeroCrossing);
        for m = 1 : length(onset_SWA{k})
            offset_SWA{k} = [offset_SWA{k} idx_ZeroCrossing{k}(find(onset_SWA{k}(m) == idx_ZeroCrossing{k}) + 1)];
        end

        for m = 1 : length(onset_SWA{k})
            amp_temp(m) = max(abs(thedata(k,onset_SWA{k}(m) : offset_SWA{k}(m))));
        end
        threshold_amp = prctile(amp_temp, keep_above_prctle);
        clear amp_temp

        m_trace = 1;
        for m = 1 : length(onset_SWA{k})

            if max(abs(thedata(k, onset_SWA{k}(m) : offset_SWA{k}(m)))) >= max(abs(threshold_amp))
                [maxi idx] = max(abs(thedata(k, onset_SWA{k}(m) : offset_SWA{k}(m))));
                if onset_SWA{k}(m) + idx - (3*sf) > 0 && onset_SWA{k}(m) + idx + (3*sf) < length(thedata(k,:))
                    SWA{k}(m_trace,:) = rawdata(k, onset_SWA{k}(m) + idx + [(-3*sf) : (3*sf)]);
                    SWA_onset{k}(m_trace) = onset_SWA{k}(m);
                    SWA_offset{k}(m_trace) = offset_SWA{k}(m);
                    SWA_middle{k}(m_trace) = onset_SWA{k}(m) + idx;
                    m_trace = m_trace + 1;
                end
            end
        end
        clear threshold_amp
    end
end

%% Look for SW of opposite polarity
clear rawdata_inv
for k = 1 : size(rawdata,1)
    rawdata_inv(k,:) = rawdata(k,:);
end
thedata_inv = filtfilt(b_low, a_low, rawdata_inv'); % filtering (high-pass)
thedata_inv = filtfilt(b_high, a_high, thedata_inv)'; % filtering (low-pass)

idx_ZeroCrossing_inv = iFindZeroCrossing(thedata_inv);
clear NegToPos_ZeroCrossing_inv
for k = 1 : length(idx_ZeroCrossing_inv)
    if thedata_inv(k,idx_ZeroCrossing_inv{k}(1) + 1) < thedata_inv(k,idx_ZeroCrossing_inv{k}(1))
        NegToPos_ZeroCrossing_inv{k} = idx_ZeroCrossing_inv{k}(2 : 2 : end);
    else
        NegToPos_ZeroCrossing_inv{k} = idx_ZeroCrossing_inv{k}(1 : 2 : end);
    end
end

% Initiate variables
onset_SWA_inv = cell(1,length(NegToPos_ZeroCrossing_inv));
offset_SWA_inv = cell(1,length(NegToPos_ZeroCrossing_inv));
SWA_inv = cell(1, length(NegToPos_ZeroCrossing_inv));
SWA_inv_onset = cell(1, length(NegToPos_ZeroCrossing_inv));
SWA_inv_offset = cell(1, length(NegToPos_ZeroCrossing_inv));
SWA_inv_middle = cell(1, length(NegToPos_ZeroCrossing_inv));

% Across electrodes
disp('Positive waves...')
clear diff_ZeroCrossing_inv
for k = 1 : length(NegToPos_ZeroCrossing_inv)
    if isempty(NegToPos_ZeroCrossing_inv{k})
        SWA_inv_onset{k} = [];
        SWA_inv_offset{k} = [];
        SWA_inv_middle{k} = [];
        SWA_inv{k} = [];
    else
        diff_ZeroCrossing_inv = diff(idx_ZeroCrossing_inv{k});
        if thedata_inv(k,idx_ZeroCrossing_inv{k}(1) + 1) < thedata_inv(k,idx_ZeroCrossing_inv{k}(1))
            NegToPos_diff_ZeroCrossing_inv = diff_ZeroCrossing_inv(2 : 2 : end);
        else
            NegToPos_diff_ZeroCrossing_inv = diff_ZeroCrossing_inv(1 : 2 : end);
        end
        onset_SWA_inv{k} = NegToPos_ZeroCrossing_inv{k}(NegToPos_diff_ZeroCrossing_inv > min_duration_ZeroCrossing...
            & NegToPos_diff_ZeroCrossing_inv < max_duration_ZeroCrossing);
        for m = 1 : length(onset_SWA_inv{k})
            offset_SWA_inv{k} = [offset_SWA_inv{k} idx_ZeroCrossing_inv{k}(find(onset_SWA_inv{k}(m) == idx_ZeroCrossing_inv{k}) + 1)];
        end

        for m = 1 : length(onset_SWA_inv{k})
            amp_temp_inv(m) = max(abs(thedata_inv(k,onset_SWA_inv{k}(m) : offset_SWA_inv{k}(m))));
        end
        threshold_amp_inv = prctile(amp_temp_inv, current_threshold);
        clear amp_temp_inv

        m_trace = 1;
        for m = 1 : length(onset_SWA_inv{k})

            if max(abs(thedata_inv(k, onset_SWA_inv{k}(m) : offset_SWA_inv{k}(m)))) >= max(abs(threshold_amp_inv))
                [maxi idx] = max(abs(thedata_inv(k, onset_SWA_inv{k}(m) : offset_SWA_inv{k}(m))));
                if onset_SWA_inv{k}(m) + idx - (3*sf) > 0 && onset_SWA_inv{k}(m) + idx + (3*sf) < length(thedata_inv(k,:))
                    SWA_inv{k}(m_trace,:) = rawdata_inv(k, onset_SWA_inv{k}(m) + idx + [(-3*sf) : (3*sf)]);
                    SWA_inv_onset{k}(m_trace) = onset_SWA_inv{k}(m);
                    SWA_inv_offset{k}(m_trace) = offset_SWA_inv{k}(m);
                    SWA_inv_middle{k}(m_trace) = onset_SWA_inv{k}(m) + idx;
                    m_trace = m_trace + 1;
                end
            end
        end
        clear threshold_amp_inv
    end
end




