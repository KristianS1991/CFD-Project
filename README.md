# Computational-Fluid-Dynamics-Project

## Timeframe
5 days

## Technologies Used
* MatLab
* Simulink

## Installation
1. Clone or download the repository.
2. Open each `.m` file MatLab.
3. Click the `Run` button in the navigation menu.

## Project Scope
The purpose of this project was to model airflow through a rotor in an airduct. Under given conditions, we were to plot various properties of the air at multiple positions along the axial direction of the duct. We were also to plot the shape of the rotor blade based on the given inlet and outlet conditions of the flow. We were to model the airflow using the Radial Equilibrium Equation (REE), obtaining the results for the compressible flow and incompressible flow cases. 

## Radial Equilibrium Equation
The REE equation is based on the premise that the radial force experienced by the air contributes to the radial movement of the flow. These forces must be balanced by the static forces exerted by the pressure gradient existing in the flow. Therefore, at any moment of time the system is assumed to be in radial equilibrium. 

#### Radial Equilibrium Equation
![REE](https://i.imgur.com/8jVC85G.png)

The majority of the calculations take place within a for-loop, iterating over each position in the air duct, calculating the constants required to satisfy the REE equation, then using these constants to calculate various properties of the air such as temperature and pressure. This took place for both the compressible and incompressible cases.

#### Flow Conditions
![Flow Conditions](https://i.imgur.com/3MRQVMk.png)

#### REE Property Calculations
![REE Property Calculations](https://i.imgur.com/CQYKvU6.png)

## Results
### Incompressible Flow Case
The incompressible flow case assumes that the density of air remains constant throughout the duct, even under changes in temperature, pressure, and velocity.

#### Incompressible Flow Results
![Incompressible Flow Results](https://i.imgur.com/NNM9BLG.png)

### Compressible Flow Case
The compressible flow case accounts for changes in the density of air under changes in pressure, temperature and velocity. 

#### Incompressible Flow Results
![Compressible Flow Results](https://i.imgur.com/8e0cGWz.png)

## Challenges


