# Code for Enhancement of the Solar Sailing Q-Law for Active Debris Removal

The repository contains the code used to obtain the results of the paper titled * *Enhancement of the Solar Sailing Q-Law for Active Debris Removal* * (see references).

The implemented framework allows the user to compute locally optimal circular-to-circular (C2C) solar-sail transfers around the Earth. The equations of motion (EOMs) are formulated in modified equinoctial elements (MEEs) to avoid singularities for near-zero eccentricity. The motion of the sailcraft is modeled using the perturbed two-body problem, where the considered perturbing accelerations are: 1) the oblateness of the Earth (J2 term of the geopotential), 2) the solar-sail acceleration (optical flat plate model), and 3) eclipses (cylindrical model).

The optimal control input is obtained by minimizing the time derivative of the so-called * *proximity quotient* * function Q. In this code, the original definition of Q is modified to remove the singularity that arises as the eccentricity control tends to zero, thereby improving convergence when the sun-orbit geometry is unfavorable.

The script `main.m` contains an example solar-sail C2C transfer that departs from an 800 km dawn-dusk sun-synchronous orbit and increases the orbital radius by 100 km. The script `SetTransferProb.m` defines the transfer parameters, while `mainPlot.m` produces six plots showing the control and evolution of the relevant parameters.

---

## Contents

* [Features](#features)
* [Install](#install)
* [Example](#example)
* [Project structure](#project-structure)
* [Authors](#authors)
* [License](#license)
* [References](#references)
* [Cite this repository](#cite-this-repository)
* [Would you like to contribute?](#would-you-like-to-contribute)

---

## Features

* Locally optimal circular-to-circular transfers in low Earth orbit based on the **solar sailing Q-Law**.
* The **primer vector** is derived from the minimization of the time-derivative of the **Q** function.
* The optimal **cone** and **clock** angles maximize the solar-sail acceleration along the **primer vector**.
* Rolling mean filter to remove the  short-period oscillations induced by the J2 perturbation.
* Fixed-step numerical integration in true longitude based on **ode45**.
* Script for setting the transfer parameters `SetTransferProb.m`.
* Script for testing `tests.m`.
* Script for plotting `mainPlot.m`.

---

## Install

> Requires MATLAB **2026a**.

```bash
# clone and install in editable mode
git clone https://github.com/fer202122/ssqlawpublic.git
cd ssqlawpublic
matlab main.m
```

---

## Example

Run `main.m` after install. The script:
* reads the transfer parameters from `SetTransferProb.m`;   
* obtains the average eccentricity of the target orbit;     
* displays the settings of the current run;
* computes the transfer by propagating the chaser;
* computes auxiliary quantities for plotting;
* tests the results by comparying against hardcoded values in `tests.m`;
* generates 6 figures by executing `mainPlot.m`:
    - Figure 1: Three-dimensional trajectory in the Sun frame (S);  
    - Figure 2: Evolution of the components of the primer unit vector in the Local frame (L);      
    - Figure 3: Evolution of the average/osculating eccentricity, aspect angle, and delta semi-major axis;
    - Figure 4: Evolution of the relative errors;
    - Figure 5: Evolution of the average/osculating classical orbital elements in the Earth frame (E);
    - Figure 6: Evolution of the cone and clock control angles in the Sun frame (S).  

---

## Project structure

```
.
├── .gitignore
├── LICENSES/
│   └── MIT.txt
├── README.md
├── SetTransferProb.m
├── main.m
├── mainPlot.m
└── src/
    ├── astrodynamics/
    │   ├── DynamicsOrbElem.m
    │   ├── FMmeoe_OblateEarth_E.m
    │   ├── FMmeoe_PointMass.m
    │   ├── FMmeoe_SRPctrl_Flat.m
    │   ├── GPE_ModEquLVLH.m
    │   ├── GaussEqClippedCoe.m
    │   ├── GaussPlanetaryEqs1.m
    │   ├── cart2mee.m
    │   ├── cmptTargetErrorLeo.m
    │   ├── coe2mee.m
    │   ├── coeDistance.m
    │   ├── coneOpt.m
    │   ├── getConstants.m
    │   ├── getSailParams.m
    │   ├── mEOEsE2P.m
    │   ├── mEOEsP2E.m
    │   ├── mee2cart.m
    │   ├── mee2coe.m
    │   ├── mee2r.m
    │   ├── orbPeriod.m
    │   ├── rotmatQCP.m
    │   ├── rotmatQEP.m
    │   ├── rotmatQIC.m
    │   ├── rotmatQLS.m
    │   ├── rotmatQPS.m
    │   ├── sunMotion.m
    │   ├── sunMotionEclipticCircular.m
    │   ├── sunSincInc.m
    │   └── symAxis.m
    ├── auxtools/
    │   ├── DoNothing.m
    │   ├── InterpArrays.m
    │   ├── QgainConstGradDeep.m
    │   ├── Rot_body.m
    │   ├── Rot_space.m
    │   ├── asRow.m
    │   ├── aspectAngle.m
    │   ├── cmptSign.m
    │   ├── disp_settings.m
    │   ├── linearAngle.m
    │   ├── makeUnit.m
    │   ├── matDotVec.m
    │   ├── no_qlaw.m
    │   ├── norm_vecs.m
    │   ├── qsol2acc.m
    │   ├── qsol2plot.m
    │   ├── rollingMean.m
    │   ├── rollingMeanRk.m
    │   ├── safeguard.m
    │   ├── splitPerRev.m
    │   ├── start_up.m
    │   ├── testIfEqual.m
    │   ├── tests.m
    │   ├── uInterp.m
    │   ├── vecCrossVec.m
    │   ├── vecDotVec.m
    │   ├── vecNorm.m
    │   └── wrapAa.m
    ├── plotting/
    │   ├── earth_sphere.m
    │   ├── plotControlHistory.m
    │   └── plotKeplerianOrbit.m
    └── solarsail/
        ├── COEnnAna.m
        ├── COEnnConst.m
        ├── SRPctrl_LyaFB_mEOE.m
        ├── attLyapunovFB.m
        ├── bubbleContourOptical.m
        ├── coneCrit.m
        ├── coneOptOguri.m
        ├── coneOptOp.m
        ├── dJdcone.m
        ├── dJdcone2.m
        ├── incCrit.m
        ├── optAttSrpOptOguri.m
        ├── primerConeOpt.m
        └── work_aSRP_hat.m
```

---

## Authors

This software was developed by:
- *Fernando Gámez Losada* ([@fer202122](https://github.com/fer202122)), ![ORCID logo](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png) [0009-0007-6107-8484](https://orcid.org/0009-0007-6107-8484), F.GamezLosada@tudelft.nl, Technische Universiteit Delft
- *Kenshiro Oguri* ([@OguriResearchGroup](https://github.com/OguriResearchGroup)), ![ORCID logo](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png) [0000-0003-3670-2293](https://orcid.org/0000-0003-3670-2293), koguri@purdue.edu, Purdue University   
  
  

---

## License

All source code files available in this repository (MATLAB scripts) are licensed under a MIT license (see `./LICENSES/MIT.txt`).  

Copyright notice:

Technische Universiteit Delft hereby disclaims all copyright interest in the program "Code for Enhancement of the Solar Sailing Q-Law for Active Debris Removal" written by the Author(s). 
Henri Werij, Faculty of Aerospace Engineering, Technische Universiteit Delft.

© 2026, F. Gámez Losada, K. Oguri

---

## References  

- [Enhancement of the Solar Sailing Q-Law for Active Debris Removal](BLA) 
- [Solar Sailing Q-Law for Planetocentric, Many-Revolution Sail Orbit Transfers](https://doi.org/10.2514/1.G007103)   

---

## Cite this repository

**How to cite this repository:** F. Gámez Losada, K. Oguri, 2026, Code for Enhancement of the Solar Sailing Q-Law for Active Debris Removal. 4TU.ResearchData. Software. [BLA](BLA)

---

## Would you like to contribute?

You are welcome to contribute! If you have any comments, feedback, or recommendations, feel free to reach out to the authors.  

If you want to contribute directly, you are welcome to open an issue and fork this repository.
