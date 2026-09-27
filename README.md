# G4H2O

## Status as of August 31, 2026

In the original implementation, the Tyvek reflector was modelled using a 100% diffuse (Lambertian cosine) reflection model, where every reflected optical photon was sampled from a cosine distribution independent of the incident angle. Although this approximation was sufficient for the original detector simulation, it did not accurately reproduce the optical response observed in the Module 2 detector. A noticeable discrepancy was found between the simulated optical response and the experimental detector data.

To improve the optical simulation for Module 2, I (Manoj Adhikari) developed and implemented a Data-Driven Tyvek Optical Reflection Model based on experimental measurements of Tyvek reflectivity in water reported by:

Álvaro Chavarría, *A Study on the Reflective Properties of Tyvek in Air and Underwater*, Department of Physics, Duke University (2007). https://phy.duke.edu/~schol/superk/alvaro_thesis.pdf

Instead of assuming purely diffuse reflection, the new implementation models Tyvek reflection as an angle-dependent combination of a diffused-specular (Gaussian) component and a Lambertian (cosine) component. The Gaussian width and the relative contributions of the specular and Lambertian components are determined from digitized experimental measurements reported by Chavarría. Consequently, the reflected photon direction depends explicitly on the incident angle, providing a significantly more realistic description of Tyvek optical behavior in water than the original Lambertian-only model.

The new optical model was validated using the complete Module 2 detector geometry. Simulated reflection-angle distributions were analysed and compared directly with the experimental measurements reported by Chavarría. For each incident angle, the simulated reflection distribution was fitted using a Gaussian + Lambertian model, and the ratio of the integrated Gaussian and Lambertian components was compared with the published experimental values. The implemented model shows good agreement with the measured Tyvek reflectivity over the full measured incident-angle range (0°–80°) while significantly improving the agreement between the Module 2 detector simulation and experimental data compared with the original 100% diffuse reflection model.

## Analysis Scripts

Analysis scripts are located inside the `mac` directory:

* **`RealdataMCdataComparision.py`** — fits real (experimental) data and Monte Carlo simulation data with cuts applied, allowing direct comparison between measured and simulated detector response.
* **`michel_parameter_scan.py`** — parameter scan script that minimizes the loss function between real and simulated data to find the best-fit simulation parameters.
* **`angular_analysis.py`** — angular analysis script used to examine the angular distribution of photons reflected from the Tyvek surface.
