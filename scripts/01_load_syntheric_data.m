%% 01_load_synthetic_data.m
% Load Synthetic_LiPo_PulseDischarge.mat from MATLAB path
% Extract time, voltage, current into workspace

clc; clear; close all;

% Load MAT file (already on MATLAB path)
S = load("Synthetic_LiPo_PulseDischarge.mat");

% Inspect contents (for safety during first run)
disp("Variables in MAT file:");
disp(fieldnames(S));

% ---- ASSUMPTION (adjust after first run if needed) ----
% Most MathWorks synthetic battery datasets use:
%   S.time
%   S.current
%   S.voltage
% If names differ, you will tell me and we’ll lock it.

time_s    = S.time(:);
current_a = S.current(:);
voltage_v = S.voltage(:);

% Basic sanity checks
assert(numel(time_s) == numel(current_a) && numel(time_s) == numel(voltage_v), ...
    "Length mismatch between time, current, voltage");

% Push to base workspace (explicit, beginner-safe)
assignin("base","time_s",time_s);
assignin("base","current_a",current_a);
assignin("base","voltage_v",voltage_v);

% Also create a table (useful later)
T = table(time_s, current_a, voltage_v);
assignin("base","T",T);

fprintf("Loaded Synthetic LiPo dataset: %d samples\n", height(T));

% Quick verification plots
figure;
subplot(2,1,1);
plot(time_s, voltage_v); grid on;
ylabel("Voltage (V)");
title("Synthetic LiPo Pulse Discharge");

subplot(2,1,2);
plot(time_s, current_a); grid on;
xlabel("Time (s)");
ylabel("Current (A)");
