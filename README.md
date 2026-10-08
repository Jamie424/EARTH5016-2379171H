# Helmsdale groundwater flow and geothermal model

A collaborative MSc Computational Geoscience project at the University of Glasgow, using MATLAB to investigate geothermal potential near Helmsdale in the Scottish Highlands.

## Overview

This 2D numerical model combines groundwater circulation with heat transport through a geological cross-section containing granite, sediments and a normal fault.

The project investigates how granite heat production and fault flow properties affect underground temperatures and potential geothermal drilling locations.

## Methods demonstrated

* Groundwater flow modelling using Darcy’s law.

* Coupled heat advection and diffusion.

* Assigning geological properties and boundary conditions.

* Parameter sensitivity analysis.

* Comparing modelled temperatures with borehole observations using RMSE.

* Numerical verification against an analytical solution and convergence tests.

## Example result

The higher granite heat-production scenario achieved a temperature RMSE of 1.86°C, compared with 2.18°C for the reference scenario.



Limited borehole observations leave fault-zone properties uncertain, so the results explore plausible scenarios.

## Running the code

Requires MATLAB with Image Processing Toolbox.

Set MATLAB’s Current Folder to usr and run Thermal_ref.m. Figures are saved to out/thermal_ref/.

Other Thermal scripts explore sensitivity scenarios. Run run_convtest_dt.m and run_convtest_dx.m for convergence tests.
