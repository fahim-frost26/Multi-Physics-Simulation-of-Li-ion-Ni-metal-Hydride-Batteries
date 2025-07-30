# Multi-Physics-Simulation-of-Li-ion-Ni-metal-Hydride-Batteries

Objectives: 
The objective of this project is to develop a MATLAB-based simulation toolkit to evaluate the condition and 
performance of lithium-ion battery and Ni-metal hydride battery using user-provided input data. The toolkit will 
integrate electrochemical reaction mechanisms, chemical kinetics, and thermal effects to simulate battery 
behavior during charging and  discharging cycles. 

The key goals of the simulation model include: 
- Capacity Estimation: Accurately simulate battery capacity during charging and discharging, accounting for 
heat generation and thermal losses.
- Degradation Analysis: Model the degradation of battery life under various load conditions and evaluate 
how temperature impacts internal resistance and long-term performance.
- Cycle Health Monitoring: Assess battery condition at any given cycle and estimate the maximum usable 
battery and remaining lifetime for that cycle.
- Thermal Effects of Fast Charging: Simulate temperature rise and its effects on battery health during fast 
charging scenarios. 

Input Parameters: 
This model will take the following input parameters from the user: 
- Load  resistances
- Initial voltage of the battery
- Total number of charging/discharging cycles to simulate
- Power ratings of the chargers used during charging cycles
  
Methodology: 
To simulate battery degradation, we will develop a MATLAB-based multi-physics model that couples 
electrochemical reaction kinetics, heat transfer, and degradation mechanisms. The simulation will take user 
inputs—including load resistances, battery voltage, initial internal resistance, ambient temperature, and 
charger power ratings—to initialize the system. 
First, we will model the electrochemical charge/discharge dynamics using coupled reaction kinetics. 
Simultaneously, the thermal module will compute heat generation (Joule heating, entropic effects) and 
dissipation (convection, conduction) to track temperature rise during cycling. Temperature at anytime t will 
depend on load resistances and we will get the temperature by solving the pde of Heat Equation in 1D.
Next, we will incorporate degradation modeling, where capacity fade and internal resistance growth are 
tracked over cycles, considering temperature-dependent aging effects such as SEI layer formation. 
The simulation will run iteratively over multiple cycles, updating degradation parameters at each step. Finally, 
we will generate performance plots—including capacity vs. cycles, temperature profiles, and lifetime 
comparisons. 

Simulation Outputs and Plots: 
The model will generate the following key plots and outputs: 
1. Battery life vs. time during discharge under load resistance per cycle 
2. Battery life vs. time during charging under load resistance per cycle 
3. Temperature vs. time during charging phases. 
4. Temperature vs. time discharging phases. 
5. Discharging time vs Optimal power  
6. Battery Health vs. number of cycles. 
7. Reaction rate

