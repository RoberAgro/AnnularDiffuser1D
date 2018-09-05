# AnnularDiffuser1D
AnnularDiffuser1D is a one-dimensional flow model for annular diffusers that allows to use arbitrary equations of state and to account for the effects of area change, heat transfer, and friction. The model is documented in a peer-reviewed, open-access publication [1] and the source code is also stored in a Zenodo repository [2].


## Requisites
AnnularDiffuser1D was implemented in MATLAB R2018a [3] and requires a REFPROP v9.1 installation [4] (see instructions).


## Instructions
* Just run the script [`AnnularDiffuser1D_demo.m`](AnnularDiffuser1D/AnnularDiffuser1D_demo.m) to get started!
* Check the function [`AnnularDiffuser1D.m`](AnnularDiffuser1D/functions/AnnularDiffuser1D.m) within the folder [functions](AnnularDiffuser1D/functions) to see the implementation of the diffuser model.
* The folder [refprop-matlab](AnnularDiffuser1D/refprop-matlab) contains a short guide about how to link REFPROP with MATLAB.


## License
AnnularDiffuser1D is licensed under the terms of the MIT license. See the [license file](LICENSE.md) for more information.


## Contact information
AnnularDiffuser1D was developed by PhD candidate [Roberto Agromayor](https://www.ntnu.edu/employees/roberto.agromayor) and Associate Professor [Lars Olof Nord](https://www.ntnu.edu/employees/lars.nord) at the [Thermal Energy Research Group](https://www.ntnu.edu/ept/thermal-energy1
) of the [Norwegian University of Science and Technology (NTNU)](https://www.ntnu.no/)

Contact the email address [roberto.agromayor@ntnu.no](mailto:roberto.agromayor@ntnu.no) for inquiries about the code.


## References
[1] R. Agromayor, and L. O. Nord, One-Dimensional Annular Diffuser Model for Preliminary Turbomachinery Design, (not published yet), 2018.

[![DOI](https://img.shields.io/badge/DOI-Diffuser_paper_DOI-blue.svg)](https://www.google.com) (not ready yet)


[2] R. Agromayor, and L. O. Nord, AnnularDiffuser1D, Zenodo repository, 2018

[![DOI](https://zenodo.org/badge/147427825.svg)](https://zenodo.org/badge/latestdoi/147427825)

[3] The MathWorks Inc., MATLAB version R2018a, 2018.

[![URL](https://img.shields.io/badge/URL-https://nl.mathworks.com/-blue.svg)](https://nl.mathworks.com/)


[4] E. W. Lemmon, M. L. Huber, and M. O. McLinden, NIST Standard Reference Database 23: Reference Fluid Thermodynamic and Transport Properties (REFPROP) version 9.1, National Institute of Standards and Technology, 2013.

[![DOI](https://img.shields.io/badge/DOI-https://dx.doi.org/10.18434/T4JS3C-blue.svg)](https://dx.doi.org/10.18434/T4JS3C)



