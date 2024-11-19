
# Energy Balance and Crop Water Stress Index (CWSI) with Actual Evapotranspiration (ETa) Calculations
This repository combines resources and codes for modeling Energy Balance (SEBAL and SAT_EB)-based Actual Evapotranspiration (ETa) and calculating Crop Water Stress Index (CWSI)-based ETa for the course `CIVE 519: Irrigation Water Management`. These simplified exercises serve as a foundation for remote sensing-based water stress and ETa modeling and can be extended for more advanced applications.

# Repository Contents
## 1. EB_ETa (Energy Balance and ETa)
This section contains simplified implementations of Surface Energy Balance Algorithm for Land (SEBAL) and Surface Aerodynamic Temperature-Energy Balance (SAT-EB). The exercise is aimed at demonstrating key concepts of energy balance approaches to estimate ETa. While basic, the framework can be expanded for functional modeling and research purposes.

<img src="EB_ETa/figs/Ts_dT_plot_SEBAL.png" height="350"/>

Fig: Ts-dT linear relationship formulation for SEBAL in EB_ETa scripts.

## 2. CWSI_ETa (Crop Water Stress Index (CWSI) and ETa)
This section provides tools to calculate CWSI using infrared thermometer (IRT) sensor data and estimate ETa based on meteorological inputs. The included scripts integrate vegetation indices, surface emissivity, and temperature corrections to highlight water stress levels in plants.


<img src="CWSI_ETa/figs/CWSI_ETa.png" height="350"/>

Fig: Estimated CWSI and ETa using CWSI_ETa scripts.

# Installation and Usage
Clone this repository to your local machine:

git clone https://github.com/mdfahimhasan/CWSI_EB_ETa.git

# Contributions
This repository is designed for academic purposes and can be expanded with additional functionalities such as soil moisture balance, Penman-Monteith ETref calculations, or advanced data assimilation.
