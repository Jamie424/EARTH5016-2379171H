# Helmsdale groundwater flow and geothermal model

A MSc Computational Geoscience project at the University of Glasgow, using MATLAB to investigate geothermal potential near Helmsdale in the Scottish Highlands.

## Overview

This 2D numerical model combines groundwater circulation with heat transport through a geological cross section containing granite, sediments and a normal fault. The project investigates how granite heat production and fault flow properties affect underground temperatures and potential geothermal drilling locations.

## Model setup

### Geological cross section

![Helmsdale geological cross-section](helmsdale-cross-section.png)

Geological vertical cross section showing the Helmsdale granite phases, surrounding geological units and fault zone represented in the model. D denotes the location of the supplied borehole data at depth, while point C denotes the proposed drill site.

### Reference model starting conditions

![Reference model starting conditions](out/thermal_ref/thermal_ref_initialfig_0.png)

Initial temperature distribution, Darcy mobility, radiogenic heat production and thermal conductivity.

## Methods demonstrated

* Groundwater flow modelling using Darcy’s law.

* Coupled heat advection and diffusion.

* Assigning geological properties and boundary conditions.

* Parameter sensitivity analysis.

* Comparing modelled temperatures with borehole observations using RMSE.

* Numerical verification against an analytical solution and convergence tests.

## Example results

### 1. Comparing predictions with borehole observations

Modelled temperature profiles were compared with supplied borehole observations. Higher granite heat production improved the fit in the scenarios tested: RMSE was 3.47°C for lower heat production, 2.18°C for the reference model and 1.86°C for higher heat production.

| Lower radiogenic heat production | Higher radiogenic heat production |
| :---: | :---: |
| ![Low heating: borehole comparison](out/thermal_lowQr/thermal_lowQr_plotdrill_16.png) | ![High heating: borehole comparison](out/thermal_uppQr/thermal_uppQr_plotdrill_16.png) |

Blue curves show modelled temperature profiles, while red points show borehole observations.

### 2. Implications for the proposed drilling site

At the proposed site, higher radiogenic heat production brought target temperatures to shallower depths. The predicted depth to 100°C decreased from 2,981 m in the lower heat production scenario 
to 2,792 m in the higher scenario.

| Lower radiogenic heat production | Higher radiogenic heat production |
| :---: | :---: |
| ![Low heating: proposed drilling site](out/thermal_lowQr/thermal_lowQr_isothermplot_16.png) | ![High heating: proposed drilling site](out/thermal_uppQr/thermal_uppQr_isothermplot_16.png) |

These are predicted temperature depth profiles at the proposed site, rather than measured borehole temperatures.

### 3. Sensitivity to Darcy mobility

Changing the fault Darcy mobility altered groundwater circulation and temperature distribution. In the scenarios tested, higher mobility increased downward transport of cooler water and slightly increased the predicted depth to 100°C at the proposed site.

| Lower Darcy mobility | Higher Darcy mobility |
| :---: | :---: |
| ![Lower Darcy mobility](out/thermal_lowKD/thermal_lowKD_isotherm2D_16.png) | ![Higher Darcy mobility](out/thermal_uppKD/thermal_uppKD_isotherm2D_21.png) |

The predicted depth to 100°C increased from 2,856 m to 2,921 m between the lower and higher mobility scenarios.

Fault zone properties remain poorly constrained as borehole observations were limited and located away from the fault. Different Darcy mobility scenarios produced similar RMSE values at the borehole, suggesting these observations provide limited information about fault zone flow. Therefore, these results explore plausible scenarios rather than establish a confirmed drilling recommendation.

## Running the code

Requires MATLAB with Image Processing Toolbox.

Set MATLAB’s Current Folder to usr and run Thermal_ref.m. Figures are saved to out/thermal_ref/.

Other Thermal scripts explore sensitivity scenarios. Run run_convtest_dt.m and run_convtest_dx.m for convergence tests.

## Limitations

Geometry and boundaries: The model represents a 2D cross section, with closed fluid flow boundaries, a fixed surface temperature and uniform basal heat flux. These assumptions simplify the regional groundwater and thermal system.

Material properties: Properties are assigned by geological unit, with radiogenic heating included only in granite. Low mobility in granite and some sediment units limits the representation of fracture flow and marine sediment circulation.

Numerical approach: Diffusion-only spin up uses Forward Euler whereas RK2 is used in the analytical verification case. Temperature updates are relaxed after Darcy flow starts. The simulations aim to approach near equilibrium conditions, so post-Darcy elapsed times should be interpreted cautiously.

## Code Origins and acknowledgements

This project began with a 1D heat transport template and supporting code snippets supplied by Dr Tobias Keller during the Numerical Geodynamics course at the University of Glasgow.

I completed and extended the implementation into a 2D heat transport and Darcy-flow model for the Helmsdale application, tracking development using Git and GitHub. The project includes parameter sensitivity studies, borehole temperature comparisons and additional visualisations.

Thanks to Tobias Keller for his teaching, guidance and feedback.
