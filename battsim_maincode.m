clc;clear;close all;
% Charging Mode Selection
% Parameters
r0 = 0.1; % Reference internal resistance (Ohm)
alpha_r = 0.004; % Resistance temp coefficient (1/°C)
eta = 0.98; % Coulombic efficiency
k_SEI = 2e-3;
Ea = 25000;
R_gas = 8.314;
T_ref = 298; % Reference temperature (K)

% Get battery type first to set appropriate ranges
batterytype = input('Enter 1 for Li-ion; 2 for Ni-metal hydride: ');
switch batterytype
    case 1
        rho = 2500; Cp = 900; Vc = 5; Vb = 3.7; C_nom = 5000; % mAh
        min_PCharge = 5; max_PCharge = 25;
        min_Rload = 0.27; max_Rload = 7.3;
    case 2
        rho = 2000; Cp = 750; Vc = 2.2; Vb = 1.6; C_nom = 2500;
        min_PCharge = 1.5; max_PCharge = 5.5;
        min_Rload = 0.22; max_Rload = 6.3;
    otherwise
        error('Invalid battery type selected.');
end

% Get and validate charging power
valid_input = false;
while ~valid_input
    PCharge = input(['Charging power(W) [' num2str(min_PCharge) '-' num2str(max_PCharge) ']: ']);
    if PCharge >= min_PCharge && PCharge <= max_PCharge
        valid_input = true;
    else
        fprintf('Charging power must be between %.1f and %.1f W\n', min_PCharge, max_PCharge);
    end
end

% Get and validate load resistance
valid_input = false;
while ~valid_input
    Rload = input(['Load resistance (Ohm) [' num2str(min_Rload) '-' num2str(max_Rload) ']: ']);
    if Rload >= min_Rload && Rload <= max_Rload
        valid_input = true;
    else
        fprintf('Load resistance must be between %.2f and %.2f Ohm\n', min_Rload, max_Rload);
    end
end

% Get and validate cycle number
valid_input = false;
while ~valid_input
    n = input('Requested cycle number [1-700]: ');
    if n >= 1 && n <= 700
        valid_input = true;
    else
        fprintf('Cycle number must be between 1 and 700\n');
    end
end

% Perfect Charger Sweep (Dynamic, Based on User Input)
if max_PCharge - min_PCharge < 5
    PCharge_sweep = linspace(min_PCharge, max_PCharge, 10);
else
    PCharge_sweep = min_PCharge:5:max_PCharge;
end
total_discharge_time_sweep = zeros(length(PCharge_sweep), 1);

for i = 1:length(PCharge_sweep)
    Ic_sweep = PCharge_sweep(i) / Vc;
    Pwaste_sweep = Ic_sweep^2 * r0 + (Vc - Vb) * Ic_sweep;
    Prcv_sweep = eta * (PCharge_sweep(i) - Pwaste_sweep);

    I_discharge_sweep = Vb / (Rload + r0);
    Pused_sweep = Vb * I_discharge_sweep;

    num_cycles_sweep = 5000;
    BH_sweep = 100 * ones(num_cycles_sweep, 1);
    t2_sweep = zeros(num_cycles_sweep, 1);
    t3_sweep = zeros(num_cycles_sweep, 1);
    TendC_sweep = zeros(num_cycles_sweep, 1);
    TendD_sweep = zeros(num_cycles_sweep, 1);
    Constant_sweep = zeros(num_cycles_sweep, 1);
    CMax_sweep = C_nom * ones(num_cycles_sweep, 1);
    QMax_sweep = C_nom * 3.6 * Vb * ones(num_cycles_sweep, 1);

    PTC_sweep = Pwaste_sweep;
    PTD_sweep = Vb^2 * r0 / Rload^2;

    k_sweep = 2;
    while BH_sweep(k_sweep-1) >= 50 && k_sweep < num_cycles_sweep
        t2_sweep(k_sweep-1) = CMax_sweep(k_sweep-1) * Vb * 3.6 * 0.8 / max(Prcv_sweep, 1e-10);
        t3_sweep(k_sweep-1) = CMax_sweep(k_sweep-1) * Vb * 3.6 * 0.8 / max(Pused_sweep, 1e-10);

        TendC_sweep(k_sweep-1) = PTC_sweep * t2_sweep(k_sweep-1)^2 / (2 * rho * Cp) + T_ref;
        TendD_sweep(k_sweep-1) = PTD_sweep * t3_sweep(k_sweep-1)^2 / (2 * rho * Cp) + T_ref;

        Constant_sweep(k_sweep-1) = exp(-Ea * (1/TendC_sweep(k_sweep-1) - 1/T_ref) / R_gas) * sqrt(t2_sweep(k_sweep-1)) + ...
            exp(-Ea * (1/TendD_sweep(k_sweep-1) - 1/T_ref) / R_gas) * sqrt(t3_sweep(k_sweep-1));

        kT_sweep = k_SEI * Constant_sweep(k_sweep-1);
        dSEI_sweep = kT_sweep;

        CMax_sweep(k_sweep) = max(CMax_sweep(k_sweep-1) - 4 * dSEI_sweep, 0);
        BH_sweep(k_sweep) = CMax_sweep(k_sweep) / (C_nom / 100);
        QMax_sweep(k_sweep) = CMax_sweep(k_sweep) * 3.6 * Vb;

        k_sweep = k_sweep + 1;
    end

    total_discharge_time_sweep(i) = sum(t3_sweep(1:k_sweep-1)) / 3600; % in hours
end

% Plot Total Discharge Time vs Charging Power
figure;
plot(PCharge_sweep, total_discharge_time_sweep, 'b-o', 'LineWidth', 2);
hold on;
[~, idx_perfect] = max(total_discharge_time_sweep);
PerfectPCharge = PCharge_sweep(idx_perfect);
PerfectDischargeTime = total_discharge_time_sweep(idx_perfect);
plot(PerfectPCharge, PerfectDischargeTime, 'rp', 'MarkerSize', 14, 'MarkerFaceColor', 'r');
text(PerfectPCharge, PerfectDischargeTime, ...
    sprintf(' Perfect: %.2f W', PerfectPCharge), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 12, 'Color', 'r');
xlabel('Charging Power (W)');
ylabel('Total Discharge Time (h)');
title('Discharge Duration vs Charging Power');
grid on;
hold off;

fprintf('Perfect charging power for your settings: %.2f W\n', PerfectPCharge);

% Ask user if they want to use the perfect charging power 
use_perfect = input(['Do you want to use the perfect charging power (' num2str(PerfectPCharge) ' W)? (1=Yes, 0=No): ']);
if use_perfect
    PCharge = PerfectPCharge;
    fprintf('Charging power set to %.2f W (perfect charger)\n', PCharge);
end

% Main Simulation

Ic = PCharge / Vc;
Vgap = Vc - Vb;

% Initial Calculations 
Pwaste = Ic^2 * r0 + Vgap * Ic;
Prcv = eta * (PCharge - Pwaste); % With Coulombic efficiency
if Prcv <= 0
    error('Charging power too low or resistance too high, resulting in no usable power.');
end

% Preallocate
num_cycles = 15000; 
BH = 100 * ones(num_cycles, 1); % Battery Health
t2 = zeros(num_cycles, 1); % Charging time per cycle
t3 = zeros(num_cycles, 1); % Discharging time per cycle
TendC = zeros(num_cycles, 1); % Final Temp charging per cycle
TendD = zeros(num_cycles, 1); % Final Temp discharging per cycle
Constant = zeros(num_cycles, 1); % (mAh) a number to decay battery capacity per cycle
CMax = C_nom * ones(num_cycles, 1); % (mAh) Battery maximum capacity
QMax = C_nom * 3.6 * Vb * ones(num_cycles, 1); % (J) Batttery Maximum Energy
k = 2;

while BH(k-1) >= 50 && k < num_cycles
    if k > 2
        Tavg = 25 + ((TendC(k-2) + TendD(k-2))/2 - T_ref);
    else
        Tavg = 25;
    end
    rT = r0 * (1 + alpha_r * (Tavg - 25));
    Ic = PCharge / Vc;
    Pwaste = Ic^2 * rT + Vgap * Ic;
    Prcv = eta * (PCharge - Pwaste);
    I_discharge = Vb / (Rload + rT);
    Pused = Vb * I_discharge;

    % Ensure minimum power values to prevent division by zero
    Prcv_safe = max(Prcv, 1e-10);
    Pused_safe = max(Pused, 1e-10);

    % Calculate times with safety checks
    t2(k-1) = CMax(k-1) * Vb * 3.6 * 0.8 / Prcv_safe;
    t3(k-1) = CMax(k-1) * Vb * 3.6 * 0.8 / Pused_safe;

    % Enforce minimum time of 10 seconds
    min_time = 10;
    if t2(k-1) <= min_time
        t2(k-1) = min_time;
    end
    if t3(k-1) <= min_time
        t3(k-1) = min_time;
    end

    PTC = Pwaste; % PTC - Power of Temperature increase for charging 
    PTD = Vb^2 * rT / Rload^2; % PTD - Power of Temperature increase for discharging 
    TendC(k-1) = PTC * t2(k-1)^2 / (2 * rho * Cp) + T_ref;
    TendD(k-1) = PTD * t3(k-1)^2 / (2 * rho * Cp) + T_ref;
    Constant(k-1) = exp(-Ea * (1/TendC(k-1) - 1/T_ref) / R_gas) * sqrt(t2(k-1)) + ...
        exp(-Ea * (1/TendD(k-1) - 1/T_ref) / R_gas) * sqrt(t3(k-1));
    kT = k_SEI * Constant(k-1);
    dSEI = kT; % For which growth of SEI later, battery gets damaged
    CMax(k) = max(CMax(k-1) - 4 * dSEI, 0); % SEI layer contributes 25% 
    BH(k) = CMax(k) / (C_nom / 100);
    QMax(k) = CMax(k) * 3.6 * Vb;
    k = k + 1;
end

if n >= k || n <= 1 % n - cycle count and k is that iteration no for which BH = 50%
    n = k - 1;
end

% Calculate reaction orders with safety checks
ord1 = Constant(1) / (t2(1) + t3(1));
ordn = Constant(n) / (t2(n) + t3(n));

% Display results
disp(['Order of reaction at first cycle: ', num2str(ord1)]);
disp(['Order of reaction at cycle ', num2str(n), ': ', num2str(ordn)]);

%  Plotting 
min_time = 10;
safe_QMax = max(QMax(n), 1);
safe_Prcv = max(Prcv, 1);
safe_Pused = max(Pused, 1);
safe_PTC = max(Pwaste, 1);
safe_PTD = max(PTD, 1);

% Ensure minimum time for plotting
if round(t2(n)) <= 1
    t2(n) = min_time;
end
if round(t3(n)) <= 1
    t3(n) = min_time;
end

% Charging curve
figure;
tcharging = 1:round(t2(n));
Qtc = safe_QMax * 0.1 + safe_Prcv * tcharging;
Qtc(Qtc > safe_QMax) = safe_QMax;
plot(tcharging / 60, Qtc * 100 / safe_QMax, 'b', 'LineWidth', 2);
xlabel('Time (min)'); ylabel('State of Charge (%)');
title('Charging Curve'); grid on;

% Discharging curve
figure;
tdischarging = 1:round(t3(n));
Qtd = safe_QMax * 0.9 - safe_Pused * tdischarging;
Qtd(Qtd < 0) = 0;
plot(tdischarging / 60, Qtd * 100 / safe_QMax, 'r', 'LineWidth', 2);
xlabel('Time (min)'); ylabel('State of Charge (%)');
title('Discharging Curve'); grid on;

% Temperature rise during charging
figure;
t = 0:1:round(t2(n));
TCC = safe_PTC * t.^2 / (2 * rho * Cp) + 25;
plot(t / 60, TCC, 'm', 'LineWidth', 2);
xlabel('Time (min)'); ylabel('Temperature (°C)');
title('Temperature Rise during Charging'); grid on;

% Temperature rise during discharging
figure;
t = 0:1:round(t3(n));
TDD = safe_PTD * t.^2 / (2 * rho * Cp) + 25;
plot(t / 60, TDD, 'k', 'LineWidth', 2);
xlabel('Time (min)'); ylabel('Temperature (°C)');
title('Temperature Rise during Discharging'); grid on;

% Battery health over cycles
figure;
cycles = 1:k-1;
colors = zeros(length(cycles), 3);
for i = 1:length(cycles)
    if BH(i) > 80
        colors(i,:) = [0 0.6 0];
    elseif BH(i) > 50
        colors(i,:) = [1 0.7 0];
    else
        colors(i,:) = [1 0 0];
    end
end
hold on;
for i = 1:length(cycles)-1
    plot(cycles(i:i+1), BH(i:i+1), 'Color', colors(i,:), 'LineWidth', 2);
end
xlabel('Cycle Number'); ylabel('Battery Health (%)');
title('Battery Health over Cycles (Color Coded)'); grid on;
hold off;

% Discharge time per cycle
figure;
plot(cycles, t3(1:k-1) / 60, 'c', 'LineWidth', 2);
xlabel('Cycle Number'); ylabel('Discharge Time (min)');
title('Discharge Time per Cycle'); grid on;

% Total discharge time
total_discharge_time = sum(t3(1:k-1)) / 3600;
disp(['Total battery discharge time (hour): ', num2str(total_discharge_time)]);

% 3D Plots: Cycle Number vs. Temperature vs. Capacity 

valid_idx = (CMax(1:k-1) > 0);

cycles_plot = cycles(valid_idx);
TendC_plot = TendC(valid_idx);
TendD_plot = TendD(valid_idx);
CMax_plot = CMax(valid_idx);

% Charging temperature
figure;
plot3(cycles_plot, TendC_plot, CMax_plot, 'bo-', 'LineWidth', 2);
xlabel('Cycle Number');
ylabel('End Charging Temp (K)');
zlabel('Capacity (mAh)');
title('3D Plot: Cycle vs. Charging Temp vs. Capacity');
grid on;
view(135, 30);

% Discharging temperature
figure;
plot3(cycles_plot, TendD_plot, CMax_plot, 'ro-', 'LineWidth', 2);
xlabel('Cycle Number');
ylabel('End Discharging Temp (K)');
zlabel('Capacity (mAh)');
title('3D Plot: Cycle vs. Discharging Temp vs. Capacity');
grid on;
view(135, 30);
