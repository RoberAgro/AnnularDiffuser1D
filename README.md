# AnnularDiffuser1D
AnnularDiffuser1D is a one-dimensional flow model for annular diffusers that allows to use arbitrary equations of state and to account for the effects of area change, heat transfer, and friction. The model is documented in a peer-reviewed, [open-access publication](https://www.google.com) and the source code is also stored in a [Zenodo repository](https://doi.org/10.5281/zenodo.1409711).


## Requisites
AxialOpt was implemented in [MATLAB R2018a](https://nl.mathworks.com/) and requires a [REFPROP v9.1](https://dx.doi.org/10.18434/T4JS3C) installation. The folder [link_refprop_matlab](link_refprop_matlab) contains instructions about how to link REFPROP with MATLAB.


## Examples
The script [`AnnularDiffuser1D_demo.m`](AnnularDiffuser1D_examples/AnnularDiffuser1D_demo.m) is a commented example to get started with the code. Check the function [`AnnularDiffuser1D.m`](AnnularDiffuser1D_source_code/AnnularDiffuser1D.m) to see the implementation of the diffuser model.


## License
AnnularDiffuser1D is licensed under the terms of the MIT license. See the [license file](LICENSE.md) for more information.


## Contact information
AnnularDifffuser1D was developed by PhD candidate [Roberto Agromayor](https://www.ntnu.edu/employees/roberto.agromayor) and Associate Professor [Lars O. Nord](https://www.ntnu.edu/employees/lars.nord) at the [Process and Power Research Group](https://www.ntnu.edu/ept/process-power#/view/about) of the [Norwegian University of Science and Technology (NTNU)](https://www.ntnu.no/)

Please, drop us an email to [roberto.agromayor@ntnu.no](mailto:roberto.agromayor@ntnu.no) if you have questions about the code or you have a bug to report. We would also love to hear about your experiences with AxialOpt in general.

## How to cite AnnularDiffuser1D?
If you want to cite AnnularDiffuser1D in a scientific publication, please refer to the [Zenodo repository](https://doi.org/10.5281/zenodo.1409711) listed in the references below and use the DOI (permanent link).


## References
R. Agromayor, and L. O. Nord, One-Dimensional Annular Diffuser Model for Preliminary Turbomachinery Design, International Journal of Turbomachinery, Propulsion and Power (submitted).

[![DOI OneDimensional_Annular_Diffuser_Model_for_Preliminary_Turbomachinery_Design](https://img.shields.io/badge/DOI-OneDimensional_Annular_Diffuser_Model_for_Preliminary_Turbomachinery_Design-blue.svg)](https://www.google.com/)


R. Agromayor, and L. O. Nord, AnnularDiffuser1D - A One-Dimensional Model for the Performance Prediction of Annular Diffusers, Zenodo repository.

[![DOI AnnularDiffuser1D_Zenodo_repository](https://img.shields.io/badge/DOI-AnnularDiffuser1D_Zenodo_repository-blue.svg)](https://doi.org/10.5281/zenodo.1409711)



E. W. Lemmon, M. L. Huber, and M. O. McLinden, NIST Standard Reference Database 23: Reference Fluid Thermodynamic and Transport Properties (REFPROP) version 9.1, National Institute of Standards and Technology, 2013.

[![DOI REFPROP_fluid_library](https://img.shields.io/badge/DOI-REFPROP_fluid_library-blue.svg)](https://dx.doi.org/10.18434/T4JS3C)




