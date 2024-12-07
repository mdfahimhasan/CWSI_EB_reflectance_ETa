# Author: Md Fahim Hasan & Rahel Pommerenke
# PhD Candidate
# Colorado State University
# Fahim.Hasan@colostate.edu

# This script is a simple replication of part of SEBAL model on how to use cold-hot pixel approach to estimate ETa.
# This is an updated version of the original SEBAL script in this repo and incorporates translation of weather station
# wind speed data to the target pixel's (crop) friction velocity, followed by hot-cold pixel dT calculation,
# H, KE, ET_hourly, ET_daily calculation.

import numpy as np
import matplotlib.pyplot as plt


def calc_heights(hc):
    """
    Calculates zero plane displacement height (d), aerodynamic roughness height for momentum transfer (z_om),
    aerodynamic roughness height for heat and vapor transfer (z_oh).

    :param hc: canopy height (m).

    :return: d, z_om, and z_oh (unit in m).
    """
    d = 0.67 * hc  # m
    z_om = 0.123 * hc  # m
    z_oh = 0.1 * z_om  # m

    return d, z_om, z_oh


def calc_friction_vel_WS(u_x, z_x, z_om_station, k=0.41):
    """
    Calculates friction velocity (u*_station) at weather station.

    :param u_x: Wind speed at weather station height (z_x) (m/s).
    :param z_x: Height of wind measurements (m) at weather station.
    :param z_om_station: Roughness length governing momentum transfer (m)).
    :param k: Von Karman's constant. Default set to 0.41.

    :return: friction velocity (m/s) at weather station.
    """
    friction_vel_ws = (k * u_x) / np.log(z_x / z_om_station)

    return friction_vel_ws


def wind_speed_blending_height(fric_vel_station, z_om_station, k=0.41):
    """
    Calculates horizontal wind speed at blending height (200m).

    :param fric_vel_station: Friction velocity at weather station height (u*_station) (m/s).
    :param z_om_station: Roughness length governing momentum transfer (m)).
    :param k: Von Karman's constant. Default set to 0.41.

    :return: friction velocity (m/s) at weather station.
    """
    u200 = fric_vel_station * np.log(200 / z_om_station) / k

    return u200


def calc_friction_vel_crop_surface(u200, z_om_crop, k=0.41):
    """
    Calculates friction velocity (u*_station) at weather station.

    :param u200: Wind speed at weather station blending height of 200 m (m/s).
    :param z_om_crop: Roughness length of crop surface governing momentum transfer (m)).
    :param k: Von Karman's constant. Default set to 0.41.

    :return: friction velocity (m/s) at weather station.
    """
    friction_vel_ws = (k * u200) / np.log(200 / z_om_crop)

    return friction_vel_ws


def calc_rah(fric_vel_crop, z_1=0.1, z_2=2, k=0.41):
    """
    Calculates Surface aerodynamic resistance (rah).

    :param fric_vel_crop: friction velocity at crop location (m/s).
    :param z_1: Near surface height of 0.1 m.
    :param z_2: Near surface height of 2 m.
    :param k: Von Karman's constant. Default set to 0.41.

    :return: Surface aerodynamic resistance (s/m).
    """
    rah = np.log(z_2 / z_1) / (fric_vel_crop * k)

    return rah


def cold_hot_dT(Rn_hot, G_hot, rah_hot, air_density=1.1, Cp=1004):
    """
    Calculate dT values for cold and hot pixels.

    # SEBAL fixes dT at two “anchor” pixels (cold and hot) and uses a linear relationship between Ts and dT.
    # Using this linear relationship, dT can be estimated using Ts for any pixel in a scene.

    :param Rn_hot: Net radiation at the hot pixel (W/m2).
    :param G_hot: Soil heat flux at the hot pixel (W/m2).
    :param rah_hot: Aerodynamic resistance to heat transport at the hot pixel (s/m) (corrected for atmospheric stability).
    :param air_density: Air density (kg/m3), default is 1.1 kg/m3.
    :param Cp: Specific heat capacity of air (J/kg/K), default is 1004 J/kg/K.

    :return: dT and H values for both cold and hot pixels (dT_cold, dT_hot, H_cold, H_hot).
    """
    # cold pixel
    H_cold = 0  # considering neutral condition (recently wetted pixel)
    dT_cold = 0

    # hot pixel
    H_hot = Rn_hot - G_hot  # considering that there is crop but not transpiring (very water stressed crop)
    dT_hot = H_hot * rah_hot / (air_density * Cp)

    return dT_cold, dT_hot, H_cold, H_hot


def Ts_dT_linear(Ts_cold, dT_cold, Ts_hot, dT_hot):
    """
    Calculates the linear relationship (slope and intercept) between Ts and dT.

    :param Ts_cold: Radiometric temperature at the cold pixel (K).
    :param dT_cold: Temperature difference at the cold pixel (assumed to be 0).
    :param Ts_hot: Radiometric temperature at the hot pixel (K).
    :param dT_hot: Temperature difference at the hot pixel (K).

    :return: Intercept and slope for the Ts-dT relationship.
    """
    slope = (dT_hot - dT_cold) / (Ts_hot - Ts_cold)
    intercept = - slope * Ts_cold

    return intercept, slope


def plot_Ts_dT_eqn(Ts_cold, dT_cold, Ts_hot, dT_hot, slope, intercept):
    """
    Plots Ts-dT linear equation.
    """
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(Ts_cold, dT_cold)
    ax.scatter(Ts_hot, dT_hot)
    ax.plot([Ts_cold, Ts_hot], [dT_cold, dT_hot], 'k-')
    ax.text(0.2, 0.7, transform=ax.transAxes, s=f'dT = {round(slope, 2)} * Ts {round(intercept, 2)}')
    plt.xlabel('Ts (K)')
    plt.ylabel('dT (K)')
    plt.savefig('figs/Ts_dT_plot_SEBAL_updated.png')


def calc_dT_crop(Ts, intercept, slope):
    """
    Calculate dT for a crop pixel using the Ts-dT linear relationship.

    :param Ts: Radiometric temperature of the crop pixel (K).
    :param intercept: Intercept from Ts-dT linear relationship.
    :param slope: Slope from Ts-dT linear relationship.

    :return: dT for the crop pixel.
    """
    dT_crop = intercept + slope * Ts

    return dT_crop


def calc_H(dT_crop, rah_crop, air_density=1.1, Cp=1004):
    """
    Calculate sensible heat flux (H) for a crop pixel during satellite overpass.

    :param dT_crop: Temperature difference for the crop pixel (K).
    :param rah_crop: Aerodynamic resistance for the crop pixel (s/m).
    :param air_density: Air density (kg/m3), default is 1.1 kg/m3.
    :param Cp: Specific heat capacity of air (J/kg/K), default is 1004 J/kg/K.

    :return: Sensible heat flux (H) for the crop pixel (W/m2).
    """
    H = (air_density * Cp * dT_crop) / rah_crop

    return H


def calc_instant_LE(Rn, G, H):
    """
    Calculate instantaneous latent heat flux (LE) for a crop pixel.

    :param Rn: Net radiation at the crop pixel (W/m2).
    :param G: Soil heat flux at the crop pixel (W/m2).
    :param H: Sensible heat flux at the crop pixel (W/m2).

    :return: Instantaneous latent heat flux (LE) for the crop pixel (W/m2).
    """
    LE = Rn - G - H

    return LE


def calc_LE_vap(Ts_cold):
    """
    Calculate latent heat of vaporization based on the temperature at the cold pixel.

    :param Ts_cold: Radiometric temperature at the cold pixel (K).

    :return: Latent heat of vaporization (LE_vap) (J/kg).
    """
    LE_vap = (2.501 - 0.00236 * (Ts_cold - 273.15)) * (10 ** 6)  # unit J/Kg

    return LE_vap


def calc_ET_inst(LE, Ts_cold, t=3600):
    """
    Calculate instantaneous evapotranspiration (ET) in mm/day.

    :param LE: Instantaneous latent heat flux (W/m2).
    :param Ts_cold: Radiometric temperature at the cold pixel (K).
    :param t: Time period (seconds), default is 3600s (1 hour).

    :return: Instantaneous ET (ET_inst) in mm/day.
    """
    LE_vap = calc_LE_vap(Ts_cold)

    ET_inst = (t * LE / (LE_vap * 1000)) * 1000  # unit in mm/hr

    return ET_inst


def calc_EF(LE, Rn_crop, G_crop):
    """
    Calculate the evaporative fraction (EF), which is the fraction of available energy used for evapotranspiration.
    This approach is based on the self-preservation theory of daytime preservation theory of daytime
    fluxes which states that the ratio between the latent heat flux (LE)
    and the available energy (Rn –G) remains constant during the day) remains constant during the day.

    :param LE: Instantaneous latent heat flux (W/m2).
    :param Rn_crop: Net radiation at the crop pixel (W/m2).
    :param G_crop: Soil heat flux at the crop pixel (W/m2).

    :return: Evaporative fraction.
    """
    EF = LE / (Rn_crop - G_crop)

    return EF


def calc_ET_daily(Ts_cold, EF, Rn_daily_avg):
    """
    Calculate daily evapotranspiration (ET) for a crop pixel.

    :param Ts_cold: Radiometric temperature at the cold pixel (K).
    :param EF: Evaporative fraction, the fraction of available energy used for evapotranspiration.
    :param Rn_daily_avg: Average daily net radiation at the crop pixel (W/m2).

    :return: Daily evapotranspiration (ET_daily) in mm/day.
    """
    LE_vap = calc_LE_vap(Ts_cold)

    # the total soil heat flux for an entire day is assumed to tend to zero for vegetation and most bare soils
    ET_daily = (86400 * EF * Rn_daily_avg / (LE_vap * 1000)) * 1000  # unit in mm/day

    return ET_daily


def calc_ET_fraction(ET_inst, ETref_hour):
    """
    Calculate the evapotranspiration fraction as the ratio of instantaneous ET to ETref.

    :param ET_inst: Instantaneous ET (mm/hr).
    :param ETref_hour: Hourly alfalfa or grass reference ET (mm/hr).

    :return: Evapotranspiration fraction.
    """
    ET_frac = ET_inst / ETref_hour

    return ET_frac


if __name__ == '__main__':
    # Provided inputs for Soybean crop
    Ts_cold = 292.7  # Cold pixel radiometric surface temperature (K)
    Ts_hot = 309.9  # Hot pixel radiometric surface temperature (K)
    Rn_hot = 410.8  # Hot pixel net radiation (W/m2)
    G_hot = 176.6  # Hot pixel soil heat flux (W/m2)
    rah_crop = 20  # Crop pixel surface aerodynamic resistance (s/m)
    Cp = 1005  # Specific heat capacity of dry air (J/kg/K)
    air_density = 1.15  # Density of air (kg/m3)
    Ts_crop = (28 + 273.15)  # Crop radiometric surface temperature (K)
    Rn_sat = 591  # Crop pixel net radiation at satellite overpass (W/m2)
    G_sat = 71  # Crop pixel soil heat flux at satellite overpass (W/m2)
    Rn_daily_avg = 256  # Daily average daily net radiation (W/m2)
    hc_grass = 0.12     # weather station grass height (m)
    u_x = 3.1           # weather station wind speed at measurement height (m/s)
    measurement_height = 2  # wind speed measurement height at weather station (m)
    hc_soybean = 0.6  # Soybean crop height (m)
    ETri = 0.65  # Instantaneous alfalfa reference ET (mm/hr)
    ETrd = 7.4  # Daily alfalfa reference ET (mm/d)
    Kcb = 1.1  # Tabulated (no stress) alfalfa-based soybean basal crop coefficient

    # Step 1: Weather station grass surface friction velocity
    _, z_om_station, _ = calc_heights(hc=hc_grass)
    fric_vel_station = calc_friction_vel_WS(u_x=u_x, z_x=measurement_height, z_om_station=z_om_station, k=0.41)

    # Step 2: Horizontal wind speed at blending height
    u200 = wind_speed_blending_height(fric_vel_station=fric_vel_station, z_om_station=z_om_station, k=0.41)

    # Step 3: Friction velocity at Soybean surface
    _, z_om_crop, _ = calc_heights(hc=hc_soybean)
    fric_vel_crop = calc_friction_vel_crop_surface(u200=u200, z_om_crop=z_om_crop, k=0.41)

    # Step 4: Aerodynamic resistance for the soybean surface
    rah_crop = calc_rah(fric_vel_crop=fric_vel_crop, z_1=0.1, z_2=2, k=0.41)

    # Step 5: Hot pixel friction velocity and surface aerodynamic resistance
    # As thi simplified model considers neutral atmospheric condition
    fric_vel_hot = fric_vel_crop
    rah_hot = rah_crop

    # Step 6: cold and hot pixel dT values
    dT_cold, dT_hot, H_cold, H_hot = cold_hot_dT(Rn_hot, G_hot, rah_hot, air_density, Cp)

    # Step 7: linear relationship (slope and intercept) between Ts and dT
    intercept, slope = Ts_dT_linear(Ts_cold, dT_cold, Ts_hot, dT_hot)

    # Step 8: plot Ts-dT relationship plot
    plot_Ts_dT_eqn(Ts_cold, dT_cold, Ts_hot, dT_hot, slope, intercept)

    # Step 9: dT for the crop pixel
    dT_crop = calc_dT_crop(Ts_crop, intercept, slope)

    # Step 10: sensible heat flux (H) for the crop pixel
    H_sat = calc_H(dT_crop = dT_crop, rah_crop=rah_crop, air_density=air_density, Cp=Cp)

    # Step 11: instantaneous latent heat flux (LE) for the crop pixel
    LE_sat = calc_instant_LE(Rn_sat, G_sat, H_sat)

    # Step 12: hourly evapotranspiration rate for the crop pixel
    ET_hourly = calc_ET_inst(LE_sat, Ts_cold, t=3600)

    # Step 13: evaporative fraction (EF) for the crop pixel
    EF = calc_EF(LE_sat, Rn_sat, G_sat)

    # Step 14: Calculate daily evapotranspiration rate for the crop pixel
    ET_daily = calc_ET_daily(Ts_cold, EF, Rn_daily_avg)

    # Step 15: Comparison between EF and alfalfa-based ET fraction (ETrf)
    ETrF = calc_ET_fraction(ET_inst=ET_hourly, ETref_hour=ETri)

    # Expected Soybean ET
    Kc = Kcb  # no water stress condition
    ETc = Kc * ETrd

    # Print results using f-strings
    print(f'a. Weather station grass surface friction velocity: {fric_vel_station:.2f} m/s \n')
    print(f'b. Horizontal wind speed at a blending height of 200 m: {u200:.2f} m/s \n')
    print(f'c. Friction velocity for the soybean surface: {fric_vel_crop:.2f} m/s \n')
    print(f'd. Aerodynamic resistance for the soybean surface: {rah_crop:.2f} s/m \n')
    print('The atmosphere is assumed to be at neutral condition, no requirement for atmospheric stability correction')
    print(f'e. Hot pixel friction velocity: {fric_vel_hot:.2f} m/s')
    print(f'and Hot pixel surface aerodynamic resistance: {rah_hot:.2f} s/m \n')
    print(f'f. dT equation (slope and intercept): slope = {slope:.2f}, intercept = {intercept:.2f} \n')
    print(f'g. dT value for the soybean surface: {dT_crop} \n')
    print(f'h. Soybean sensible heat flux: {H_sat:.2f} W/m2 \n')
    print(f'i. Soybean latent heat flux: {LE_sat:.2f} W/m2 \n')
    print(f'j. Soybean instantaneous actual water use rate (ETai): {ET_hourly:.2f} mm/hr \n')
    print(f'k. Evaporative fraction for Soybean surface (EF): {EF:.2f} \n')
    print(f'l. Alfalfa-based reference ET fraction (ETrf): {ETrF:.2f} \n')
    print(f'm. Explanation of difference between EF and ETrf \n')
    print(f'n. Soybean daily actual evapotranspiration rate (ETad): {ET_daily:.2f} mm/d \n')
    # print(f'o. Soybean potential ET: {:.2f} mm/d, ETc: {:.2f} mm/d')