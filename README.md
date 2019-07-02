# Computational Fluid Dynamics Project - Free Vortex Design

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

The majority of the calculations take place within a `for` loop, iterating over each position in the air duct, calculating the constants required to satisfy the REE equation, then using these constants to calculate various properties of the air such as temperature and pressure. This took place for both the compressible and incompressible cases.

#### Flow Conditions
![Flow Conditions](https://i.imgur.com/3MRQVMk.png)

#### REE Property Calculations
![REE Property Calculations](https://i.imgur.com/CQYKvU6.png)

## Results
### Incompressible Flow Case
The incompressible flow case assumes that the density of air remains constant throughout the duct, even under changes in temperature, pressure, and velocity. The first diagram below shows the relationship between the density of the air and the radial distance from the center of the duct. As this case is incompressible, the density remains constant at any radial distance from the center. The second diagram shows the relationship between beta, the inlet angle of the airflow on the rotor blade, and the radial distance from the center of the rotor. As the distance from the center of the rotor increases, the inlet angle beta decreases. The third diagram is a 3-D plot of the rotor blade, based on the beta values of the air in a 3-D spatial geometry.

#### Incompressible Flow Results
![Incompressible Flow Results](https://i.imgur.com/NNM9BLG.png)

### Compressible Flow Case
The compressible flow case accounts for changes in the density of air under changes in pressure, temperature and velocity. The first diagram below shows the relationship between the density of the air and the radial distance from the center of the duct. In the compressible case, the density of air varies depending on the radial distance from the center. The density of air increases as the radial distance increases. This is because air at the tips of the rotor blades has a lower velocity and is under higher pressure than air at the center of the rotor blades. The second diagram shows the relationship between beta, the inlet angle of the airflow on the rotor blade, and the radial distance from the center of the rotor. As the distance from the center of the rotor increases, the inlet angle beta decreases. The third diagram is a 3-D plot of the rotor blade, based on the beta values of the air in a 3-D spatial geometry.

#### Incompressible Flow Results
![Compressible Flow Results](https://i.imgur.com/8e0cGWz.png)

## Challenges
Calculating properties of air at a singular point is simple if you are given the necessary conditions at that point. Modelling these properties across a whole 3-D spatial geometry (the air duct) proves to be a lot more difficult and requires a series of nested `for` loops to account for each dimension. Ensuring that each loop is nested properly, and boundary conditions are applied and met can be very challenging when so many properties have to be calculated. There is a lot of room for error when the whole model depends on the syntax of multiple complex formulas. Completing this project took an iterative trial and error process, involving running the code, logging results at certain points, and checking that the graphs lined up with the theory.

