%% 02_segment_pulses_and_soc.m
% Break pulses using discharge pulse starts (from current)
% Segment voltage between discharge pulse starts
% Compute SOC from current with SOC(0)=100% and SOC(end)=0% (normalized integration)


% Make sure Step 1 ran (or load it here)
if ~exist("time_s","var") || ~exist("current_a","var") || ~exist("voltage_v","var")
    run("scripts/01_load_synthetic_data.m");
end

% ---- SETTINGS (you control these) ----
dischargeIsPositive = true;   % set false if discharge pulses are negative in your current signal
pulseThresholdFrac  = 0.10;   % 10% of max discharge current to detect pulse start
minGap_s            = 0.5;    % debounce: ignore pulse starts closer than this
assignSocPoint      = "mid";  % "start" | "mid" | "end"

% Put in a table
T = table(time_s(:), current_a(:), voltage_v(:), 'VariableNames', {'time_s','current_a','voltage_v'});

%% 1) Compute SOC inline (SOC start=100, end=0)
t = double(T.time_s);
i = double(T.current_a);

% Ensure time monotonic (safe)
[~, idx] = sort(t);
t = t(idx); i = i(idx);
T = T(idx,:);

% Convert to discharge-positive for SOC integration
if dischargeIsPositive
    i_dis = i;
else
    i_dis = -i;
end

% Only count discharge (ignore any negative/charge parts)
i_dis(i_dis < 0) = 0;

% Integrate discharged charge (Coulombs)
Q_c = cumtrapz(t, i_dis);     % A*s

Q_total = Q_c(end);
if Q_total <= 0
    error("Total discharged charge is zero. Flip dischargeIsPositive or check current signal.");
end

% Normalize to SOC 100->0 across dataset
soc_pct = 100 * (1 - (Q_c ./ Q_total));
soc_pct = max(0, min(100, soc_pct));  % clamp for numerical safety

T.soc_pct = soc_pct;

%% 2) Find discharge pulse starts (rising edge above threshold)
% Discharge-positive current for detection
if dischargeIsPositive
    i_det = T.current_a;
else
    i_det = -T.current_a;
end

thr = pulseThresholdFrac * max(abs(i_det));
isPulse = i_det > thr;

d = diff([false; isPulse]);
startIdx = find(d == 1);

% Debounce close starts
if ~isempty(startIdx)
    keep = true(size(startIdx));
    for k = 2:numel(startIdx)
        if (T.time_s(startIdx(k)) - T.time_s(startIdx(k-1))) < minGap_s
            keep(k) = false;
        end
    end
    startIdx = startIdx(keep);
end

%% 3) Segment from each start to just before the next start
segments = struct('idxStart',{},'idxEnd',{}, ...
                  'time_s',{},'current_a',{},'voltage_v',{}, ...
                  'soc_start',{},'soc_end',{},'soc_mid',{},'soc_point',{});

for k = 1:numel(startIdx)
    idxStart = startIdx(k);
    if k < numel(startIdx)
        idxEnd = startIdx(k+1) - 1;
    else
        idxEnd = height(T);
    end

    segments(k).idxStart  = idxStart;
    segments(k).idxEnd    = idxEnd;
    segments(k).time_s    = T.time_s(idxStart:idxEnd);
    segments(k).current_a = T.current_a(idxStart:idxEnd);
    segments(k).voltage_v = T.voltage_v(idxStart:idxEnd);

    segments(k).soc_start = T.soc_pct(idxStart);
    segments(k).soc_end   = T.soc_pct(idxEnd);
    segments(k).soc_mid   = T.soc_pct(round((idxStart + idxEnd)/2));

    switch lower(string(assignSocPoint))
        case "start"
            segments(k).soc_point = segments(k).soc_start;
        case "end"
            segments(k).soc_point = segments(k).soc_end;
        otherwise
            segments(k).soc_point = segments(k).soc_mid;
    end
end

%% 4) Push outputs to workspace
assignin("base","T",T);
assignin("base","segments",segments);
assignin("base","pulse_start_idx",startIdx);

fprintf("Detected %d discharge pulse starts -> %d segments.\n", numel(startIdx), numel(segments));

%% 5) Plots for sanity
figure; plot(T.time_s, T.current_a); grid on; hold on;
xline(T.time_s(startIdx));
xlabel("Time (s)"); ylabel("Current (A)");
title("Detected discharge pulse starts");

figure; plot(T.time_s, T.voltage_v); grid on; hold on;
xline(T.time_s(startIdx));
xlabel("Time (s)"); ylabel("Voltage (V)");
title("Voltage with segment boundaries");

figure; plot(T.time_s, T.soc_pct); grid on;
xlabel("Time (s)"); ylabel("SOC (%)");
title("SOC normalized from 100% (start) to 0% (end)");

