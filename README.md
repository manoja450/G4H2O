
## Status as of August 4, 2026

This repository was originally created and developed by **Igor Bernardi** during his Ph.D. work with the COHERENT Collaboration. The overall folder and file structure is believed to have originated from **Matthew Blackstone**, while **Karla Tellez** also contributed to the development of the simulation framework during her time in the collaboration.

The complete detector geometry—including essentially all detector components except the photomultiplier tubes—was developed by Igor Bernardi and represents the detector as built using the best available engineering dimensions. The original geometry and optical simulation were developed for the **Module 1 (D₂O + H₂O)** detector.

In the original implementation, the Tyvek reflector was modelled using a **100% diffuse (Lambertian cosine)** reflection model, where every reflected optical photon was sampled from a cosine distribution independent of the incident angle. Although this approximation was sufficient for the original detector simulation, it did not accurately reproduce the optical response observed in the **Module 2** detector. A noticeable discrepancy was found between the simulated optical response and the experimental detector data.

To improve the optical simulation for Module 2, I (**Manoj Adhikari**) developed and implemented a **Data-Driven Tyvek Optical Reflection Model** based on experimental measurements of Tyvek reflectivity in water reported by:

> **Álvaro Chavarría**
> *A Study on the Reflective Properties of Tyvek in Air and Underwater*
> Department of Physics, Duke University (2007)
> https://phy.duke.edu/~schol/superk/alvaro_thesis.pdf

Instead of assuming purely diffuse reflection, the new implementation models Tyvek reflection as an **angle-dependent combination of a diffused-specular (Gaussian) component and a Lambertian (cosine) component**. The Gaussian width and the relative contributions of the specular and Lambertian components are determined from digitized experimental measurements reported by Chavarría. Consequently, the reflected photon direction depends explicitly on the incident angle, providing a significantly more realistic description of Tyvek optical behavior in water than the original Lambertian-only model.

The implementation is primarily contained in:

* `G4d2oDataDrivenReflector.cc`
* `G4d2oDataDrivenReflector.hh`
* `G4d2oCustomOpBoundary.cc`
* `G4d2oMaterialsDefinition::SetDataDrivenReflector()`

The new optical model was validated using the complete **Module 2** detector geometry. Simulated reflection-angle distributions were analysed and compared directly with the experimental measurements reported by Chavarría. For each incident angle, the simulated reflection distribution was fitted using a **Gaussian + Lambertian** model, and the ratio of the **integrated Gaussian and Lambertian components** was compared with the published experimental values. The implemented model shows good agreement with the measured Tyvek reflectivity over the full measured incident-angle range (0°–80°) while significantly improving the agreement between the Module 2 detector simulation and experimental data compared with the original 100% diffuse reflection model.
