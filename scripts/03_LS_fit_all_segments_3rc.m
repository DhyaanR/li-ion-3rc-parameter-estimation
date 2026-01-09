%% 03_LS_fit_all_segments_3rc.m
% Fit 3RC ECM for ALL segments using nonlinear least squares.
% Plot each R and C vs SOC individually.
% Output final table: SOC and corresponding parameter values.

addpath(genpath(pwd));

% ---- Fit settings ----
cfg = struct();
cfg.minN = 60;                 % minimum samples per segment
cfg.maxFunEvals = 25000;
cfg.maxIters = 1200;
cfg.display = "off";           % set "iter" to watch convergence
cfg.useWarmStart = true;       % reuse last fit to speed convergence
cfg.skipIfNoExcitation = true; % skip near-zero current segments
cfg.minCurrentA = 0.02;        % A

opts = optimoptions("lsqnonlin", ...
    "Display", cfg.display, ...
    "MaxFunctionEvaluations", cfg.maxFunEvals, ...
    "MaxIterations", cfg.maxIters);

% ---- Initial guess (R/C) ----
R0_0 = 0.02;
R1_0 = 0.01;  C1_0 = 2000;     % tau ~ 20 s
R2_0 = 0.02;  C2_0 = 5000;     % tau ~ 100 s
R3_0 = 0.05;  C3_0 = 20000;    % tau ~ 1000 s

% Bounds (physical)
% x = [V0 R0 R1 C1 R2 C2 R3 C3]
lb_base = [-inf, 1e-6, 1e-6, 1,    1e-6, 1,    1e-6, 1];
ub_base = [ inf,  1,    1,   1e7,  1,    1e7,  1,    1e8];

% Storage matrix columns:
% soc, R0, R1, R2, R3, C1, C2, C3, V0, rmse, maxErr, exitflag
rows = zeros(0, 12);
skipped = 0;

x_prev = []; % warm start cache

for k = 1:numel(segments)
    seg = segments(k);
    t = double(seg.time_s(:));
    i = double(seg.current_a(:));
    v = double(seg.voltage_v(:));

    if numel(v) < cfg.minN
        skipped = skipped + 1;
        continue;
    end

    if cfg.skipIfNoExcitation && max(abs(i)) < cfg.minCurrentA
        skipped = skipped + 1;
        continue;
    end

    % Normalize time
    t = t - t(1);

    % Segment-specific V0 initial guess (offset)
    V0_0 = v(end);

    x0 = [V0_0, R0_0, R1_0, C1_0, R2_0, C2_0, R3_0, C3_0];

    if cfg.useWarmStart && ~isempty(x_prev)
        x0(2:end) = x_prev(2:end);
        x0(1) = V0_0;
    end

    lb = lb_base;
    ub = ub_base;
    lb(1) = min(v) - 2;
    ub(1) = max(v) + 2;

    try
        fun = @(x) (v - simulate_3rc(x, t, i));  % residual vector
        [xhat, ~, ~, exitflag] = lsqnonlin(fun, x0, lb, ub, opts);

        vhat = simulate_3rc(xhat, t, i);

        rmse = sqrt(mean((v - vhat).^2));
        maxerr = max(abs(v - vhat));

        % Physical sanity: skip non-positive R/C
        if any(xhat([2 3 5 7]) <= 0) || any(xhat([4 6 8]) <= 0)
            skipped = skipped + 1;
            continue;
        end

        soc = seg.soc_point;

        % Unpack
        V0 = xhat(1);
        R0 = xhat(2);
        R1 = xhat(3); C1 = xhat(4);
        R2 = xhat(5); C2 = xhat(6);
        R3 = xhat(7); C3 = xhat(8);

        rows(end+1,:) = [soc, R0, R1, R2, R3, C1, C2, C3, V0, rmse, maxerr, exitflag]; %#ok<AGROW>

        x_prev = xhat;

    catch
        skipped = skipped + 1;
        continue;
    end
end

assert(~isempty(rows), "No segments were fitted. Check segmentation or excitation threshold.");

% Build final table (as you requested)
paramsTbl = array2table(rows, 'VariableNames', ...
    {'soc_pct','R0_ohm','R1_ohm','R2_ohm','R3_ohm','C1_F','C2_F','C3_F','V0_V','rmse_V','maxAbsErr_V','exitflag'});

% Sort by SOC descending
paramsTbl = sortrows(paramsTbl, "soc_pct", "descend");

% Push to workspace
assignin("base","paramsTbl",paramsTbl);
assignin("base","cfg_fit",cfg);

fprintf("Done: fitted %d segments, skipped %d.\n", height(paramsTbl), skipped);

%% ---- Plots (each parameter vs SOC individually) ----
soc = paramsTbl.soc_pct;

% R plots (individual)
figure; plot(soc, paramsTbl.R0_ohm, "-o"); grid on;
xlabel("SOC (%)"); ylabel("R0 (ohm)"); title("R0 vs SOC");

figure; plot(soc, paramsTbl.R1_ohm, "-o"); grid on;
xlabel("SOC (%)"); ylabel("R1 (ohm)"); title("R1 vs SOC");

figure; plot(soc, paramsTbl.R2_ohm, "-o"); grid on;
xlabel("SOC (%)"); ylabel("R2 (ohm)"); title("R2 vs SOC");

figure; plot(soc, paramsTbl.R3_ohm, "-o"); grid on;
xlabel("SOC (%)"); ylabel("R3 (ohm)"); title("R3 vs SOC");

% C plots (individual)
figure; plot(soc, paramsTbl.C1_F, "-o"); grid on;
xlabel("SOC (%)"); ylabel("C1 (F)"); title("C1 vs SOC");

figure; plot(soc, paramsTbl.C2_F, "-o"); grid on;
xlabel("SOC (%)"); ylabel("C2 (F)"); title("C2 vs SOC");

figure; plot(soc, paramsTbl.C3_F, "-o"); grid on;
xlabel("SOC (%)"); ylabel("C3 (F)"); title("C3 vs SOC");

% Optional: fit quality plots
figure; plot(soc, paramsTbl.rmse_V, "-o"); grid on;
xlabel("SOC (%)"); ylabel("RMSE (V)"); title("Fit RMSE vs SOC");

figure; plot(soc, paramsTbl.maxAbsErr_V, "-o"); grid on;
xlabel("SOC (%)"); ylabel("Max Abs Error (V)"); title("Max Abs Error vs SOC");

%% ---- Local function: 3RC simulation ----
function v = simulate_3rc(x, t, i)
    % x = [V0 R0 R1 C1 R2 C2 R3 C3]
    V0 = x(1); R0 = x(2);
    R1 = x(3); C1 = x(4);
    R2 = x(5); C2 = x(6);
    R3 = x(7); C3 = x(8);

    n = numel(t);
    v1 = 0; v2 = 0; v3 = 0;
    v  = zeros(n,1);

    for k = 1:n
        Ik = i(k);
        v(k) = V0 - Ik*R0 - v1 - v2 - v3;

        if k == n, break; end
        dt = t(k+1) - t(k);

        a1 = exp(-dt/(R1*C1));  v1 = a1*v1 + R1*(1-a1)*Ik;
        a2 = exp(-dt/(R2*C2));  v2 = a2*v2 + R2*(1-a2)*Ik;
        a3 = exp(-dt/(R3*C3));  v3 = a3*v3 + R3*(1-a3)*Ik;
    end
end
