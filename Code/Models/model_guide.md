## OCV (Open Circuit Voltage)
- Most basic model(not connected to any external circuit)
- Parameters: SOC, voltage
    - Qn: negative electrode capacity (As)
    - nu: negative/positive electrode capacity ratio
    - miu: cyclable lithium/positive electrode capacity ratio
  
    - Q: effective cell capacity(=Qn)
- Usage: basic SOC estimation

## OCVT (Open Circuit Voltage with Temperature)
- OCV + Temperature
- Parameters: SOC, voltage, temperature
    - Cp: heat capacity of the core
    - Cps: heat capacity of the surface
    - tauT: internal heat transfer timescale
    - tauA: external heat transfer timescale
    - CE: coulombic efficiency
    - Q: effective negative electrode capacity(=Qn/CE)
- Usage: when temperature effects non-negligible, not suitable for dynamic simulation

## ROCV (Resistance Open Circuit Voltage)
- OCV + Series Resistance
- Parameters: SOC, voltage, series resistance
    - Rs: series resistance
- Usage: when taking into account internal resistance

## RORC (Resistance, Open Circuit Voltage, Resistance, Capacitance)
- OCV + Series Resistance, RC pair
- Parameters: SOC, voltage, series resistance, RC pair parameters(resistance, capacitance)
    - tau1: time constant of RC pair(=R1*C1)
    - C1: capacitance of RC pair
- Usage: transient behaviour(more complex)
  
---
## EHM (Equivalent Hydraulic Model)
- Parameters:
    - tau_ref: diffusion time constant
    - b: negative electrode surface/particle volume ratio(higher b=larger reaction surface)
    - Ip_ref: reference exchange currents for positive and negative electrodes
    - Rf: film resistance
- Constants:
    - alph: charge transfer coefficients
    - Faraday: Faraday's constant
    - Rg: gas constant
- Activation Energies
    - E_Dsn, E_kn, E_kp
- Usage: electrochemical, diffusion effects
    
## EHMT (Equivalent Hydraulic Model with Temperature)
- EHM + Thermal dynamics(heat during charge/discharge)
- Usage: thermal management