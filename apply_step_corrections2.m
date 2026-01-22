function results = apply_step_corrections2(results, participant_id)
%APPLY_STEP_CORRECTIONS Apply manual step-detection corrections for one participant.
%
% This function applies participant-specific corrections (add/remove) to the
% detected step indices stored in:
%   results.(participant_id).steps.(foot).(trial)
%
% Inputs
%   results        - Struct containing detected steps per participant/foot/trial.
%   participant_id - Participant identifier string (e.g., 'P7').
%
% Output
%   results        - Same struct, with corrected step indices written back.

%% -------------------- Step 1: Define correction rules --------------------
% Corrections are defined as sample indices to add/remove for specific
% participant/foot/trial combinations.
corrections = struct();

% LFoot - Ascent1
corrections.P37.LFoot.Ascent1.add    = [221];
corrections.P30.LFoot.Ascent1.remove = [866];
corrections.P29.LFoot.Ascent1.add    = [240];
corrections.P23.LFoot.Ascent1.add    = [227];
corrections.P10.LFoot.Ascent1.add    = [249];

% LFoot - Ascent2
corrections.P33.LFoot.Ascent2.add    = [238];
corrections.P29.LFoot.Ascent2.add    = [242];
corrections.P10.LFoot.Ascent2.add    = [234, 474, 549, 627];
corrections.P10.LFoot.Ascent2.remove = [485, 448, 559, 521, 595, 631];
corrections.P3.LFoot.Ascent2.add     = [234];

% LFoot - Descent1
corrections.P36.LFoot.Descent1.add    = [565];
corrections.P36.LFoot.Descent1.remove = [591];
corrections.P35.LFoot.Descent1.add    = [356, 776];
corrections.P31.LFoot.Descent1.add    = [238];
corrections.P29.LFoot.Descent1.add    = [236];
corrections.P22.LFoot.Descent1.add    = [244];
corrections.P19.LFoot.Descent1.add    = [1163, 855];
corrections.P7.LFoot.Descent1.add     = [233];

% LFoot - Descent2
corrections.P36.LFoot.Descent2.add    = [913];
corrections.P36.LFoot.Descent2.remove = [943];
corrections.P19.LFoot.Descent2.add    = [1048, 725];

% RFoot - Ascent1
corrections.P39.RFoot.Ascent1.add = [252];
corrections.P36.RFoot.Ascent1.add = [235];
corrections.P24.RFoot.Ascent1.add = [260];
corrections.P24.RFoot.Ascent1.remove = [286];
corrections.P12.RFoot.Ascent1.add = [233];
corrections.P4.RFoot.Ascent1.add  = [251];
corrections.P2.RFoot.Ascent1.add  = [238];
corrections.P1.RFoot.Ascent1.add  = [258];

% RFoot - Ascent2
corrections.P34.RFoot.Ascent2.add = [231];
corrections.P30.RFoot.Ascent2.add = [245];
corrections.P30.RFoot.Ascent2.remove = [271];
corrections.P27.RFoot.Ascent2.add = [240];
corrections.P27.RFoot.Ascent2.remove = [264];
corrections.P24.RFoot.Ascent2.add = [239];
corrections.P24.RFoot.Ascent2.remove = [264];
corrections.P22.RFoot.Ascent2.add = [234];
corrections.P22.RFoot.Ascent2.remove = [258];
corrections.P19.RFoot.Ascent2.remove = [254, 667];
corrections.P15.RFoot.Ascent2.add = [201];
corrections.P12.RFoot.Ascent2.add = [235];
corrections.P12.RFoot.Ascent2.remove = [259];

% RFoot - Descent1
corrections.P37.RFoot.Descent1.remove = [1153];
corrections.P35.RFoot.Descent1.add    = [956];
corrections.P34.RFoot.Descent1.remove = [273, 904];
corrections.P26.RFoot.Descent1.add    = [220];
corrections.P24.RFoot.Descent1.add    = [247];
corrections.P17.RFoot.Descent1.add    = [875, 929];
corrections.P9.RFoot.Descent1.add     = [220];
corrections.P6.RFoot.Descent1.add     = [377, 471, 550];
corrections.P2.RFoot.Descent1.add     = [419];
corrections.P2.RFoot.Descent1.remove  = [441];

% RFoot - Descent2
corrections.P36.RFoot.Descent2.add    = [243, 313, 379];
corrections.P26.RFoot.Descent2.add    = [697];
corrections.P26.RFoot.Descent2.remove = [718];
corrections.P24.RFoot.Descent2.add    = [244];
corrections.P17.RFoot.Descent2.add    = [917];
corrections.P6.RFoot.Descent2.add     = [435];
corrections.P5.RFoot.Descent2.add     = [249];
corrections.P3.RFoot.Descent2.add     = [842];

%% -------------------- Step 2: Apply corrections --------------------
% Only apply corrections if this participant exists in the corrections struct.
if isfield(corrections, participant_id)

    feet_corr = fieldnames(corrections.(participant_id));
    for f = 1:length(feet_corr)
        foot = feet_corr{f};

        trials_corr = fieldnames(corrections.(participant_id).(foot));
        for t = 1:length(trials_corr)
            trial = trials_corr{t};

            % Verify that original step data exist before attempting edits.
            if isfield(results.(participant_id).steps.(foot), trial)
                steps = results.(participant_id).steps.(foot).(trial);

                % Add missing steps (union + unique)
                if isfield(corrections.(participant_id).(foot).(trial), 'add')
                    steps = unique([steps, corrections.(participant_id).(foot).(trial).add]);
                end

                % Remove incorrect steps, with a ±1 sample tolerance check.
                if isfield(corrections.(participant_id).(foot).(trial), 'remove')
                    to_remove = corrections.(participant_id).(foot).(trial).remove;

                    for i = 1:length(to_remove)
                        val = to_remove(i);

                        if any(steps == val)
                            steps = steps(steps ~= val);
                            fprintf('✔ Removed exact %d in %s - %s - %s\n', ...
                                    val, participant_id, foot, trial);

                        elseif any(steps == val + 1)
                            steps = steps(steps ~= val + 1);
                            warning('✔ %d not found, but %d was present → removed in %s - %s - %s', ...
                                    val, val + 1, participant_id, foot, trial);

                        elseif any(steps == val - 1)
                            steps = steps(steps ~= val - 1);
                            warning('✔ %d not found, but %d was present → removed in %s - %s - %s', ...
                                    val, val - 1, participant_id, foot, trial);

                        else
                            warning('✘ %d (and ±1) not found in %s - %s - %s', ...
                                    val, participant_id, foot, trial);
                        end
                    end
                end

                % Write corrected steps back to the results struct.
                results.(participant_id).steps.(foot).(trial) = steps;

            else
                warning('No original step data found for %s - %s - %s', ...
                        participant_id, foot, trial);
            end
        end
    end
end
end
