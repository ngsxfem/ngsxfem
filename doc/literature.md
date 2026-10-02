---
Scientific literature using `ngsxfem`
---

<!-- This file is generated from doc/literature.yaml by doc/make_literature.py. Do not edit by hand. -->

This list collects scientific works (journal articles, preprints, theses) in which `ngsxfem` has been used. Entries are sorted by category and year (newest first); where available, links to reproduction code and data are given. If you used `ngsxfem` in your work and it is missing here, please open an issue or a pull request on [GitHub](https://github.com/ngsxfem/ngsxfem/issues) (the list is generated from `doc/literature.yaml`).

### Citing `ngsxfem`

If you use `ngsxfem` for your research, please cite: C. Lehrenfeld, F. Heimann, J. Preuß and H. von Wahl. “ngsxfem: Add-on to NGSolve for geometrically unfitted finite element discretizations”. *Journal of Open Source Software* 6(64), 3237, 2021. doi: [10.21105/joss.03237](https://doi.org/10.21105/joss.03237). Software releases are archived on Zenodo: doi: [10.5281/zenodo.5081124](https://doi.org/10.5281/zenodo.5081124).

```
@article{ngsxfem,
  author  = {Lehrenfeld, Christoph and Heimann, Fabian and Preu{\ss}, Janosch and von Wahl, Henry},
  title   = {ngsxfem: Add-on to NGSolve for geometrically unfitted finite element discretizations},
  journal = {Journal of Open Source Software},
  volume  = {6},
  number  = {64},
  pages   = {3237},
  year    = {2021},
  doi     = {10.21105/joss.03237}
}
```

### Unfitted boundary value problems

Fictitious domain, CutFEM, unfitted DG and Trefftz methods for boundary value problems on level set domains, including the analysis of the higher order geometry handling.

* C. Lehrenfeld, T. van Beeck and I. Voulis. “Analysis of divergence-preserving unfitted finite element methods for the mixed Poisson problem”. *Math. Comp.* 94(354), 1667–1699, 2025. doi: [10.1090/mcom/4027](https://doi.org/10.1090/mcom/4027), arXiv: [2306.12722](https://arxiv.org/abs/2306.12722). Code: [GRO.data (doi: 10.25625/LYR2CG)](https://doi.org/10.25625/LYR2CG).

* M. Neilan, M. A. Olshanskii and H. von Wahl. “An unfitted divergence-free higher order finite element method for the Stokes problem”. preprint, 2025. arXiv: [2512.12050](https://arxiv.org/abs/2512.12050). Code: [GitHub](https://github.com/hvonwah/repro-unfitted-div-free-stokes), [Zenodo (doi: 10.5281/zenodo.21460130)](https://doi.org/10.5281/zenodo.21460130).

* E. N. Karatzas. “hp-version analysis for arbitrarily shaped elements on the boundary discontinuous Galerkin method for Stokes systems”. *Int. J. Numer. Anal. Model.* 21(4), 528–559, 2024. doi: [10.4208/ijnam2024-1021](https://doi.org/10.4208/ijnam2024-1021), arXiv: [2301.12577](https://arxiv.org/abs/2301.12577).

* F. Heimann, C. Lehrenfeld, P. Stocker and H. von Wahl. “Unfitted Trefftz discontinuous Galerkin methods for elliptic boundary value problems”. *ESAIM Math. Model. Numer. Anal.* 57(5), 2803–2833, 2023. doi: [10.1051/m2an/2023064](https://doi.org/10.1051/m2an/2023064), arXiv: [2212.12236](https://arxiv.org/abs/2212.12236). Code: [GitHub](https://github.com/hvonwah/unf-trefftz-poisson-code), [Zenodo (doi: 10.5281/zenodo.7474688)](https://doi.org/10.5281/zenodo.7474688).

* H. Liu, M. Neilan and M. A. Olshanskii. “A CutFEM divergence-free discretization for the Stokes problem”. *ESAIM Math. Model. Numer. Anal.* 57(1), 143–165, 2023. doi: [10.1051/m2an/2022072](https://doi.org/10.1051/m2an/2022072), arXiv: [2110.11456](https://arxiv.org/abs/2110.11456).

* A. Aretaki, E. N. Karatzas and G. Katsouleas. “Equal higher order analysis of an unfitted discontinuous Galerkin method for Stokes flow systems”. *J. Sci. Comput.* 91(2), 2022. doi: [10.1007/s10915-022-01823-w](https://doi.org/10.1007/s10915-022-01823-w), arXiv: [2006.00435](https://arxiv.org/abs/2006.00435).

* H. Liu. “Unfitted finite element methods for the Stokes problem using the Scott-Vogelius pair”. *PhD thesis* University of Pittsburgh, 2022. url: <http://d-scholarship.pitt.edu/43584/>.

* F. Heimann and C. Lehrenfeld. “Numerical integration on hyperrectangles in isoparametric unfitted finite elements”. *Numerical Mathematics and Advanced Applications ENUMATH 2017* Lecture Notes in Computational Science and Engineering 126, Springer, 193–202, 2019. doi: [10.1007/978-3-319-96415-7_16](https://doi.org/10.1007/978-3-319-96415-7_16).

* C. Lehrenfeld. “A higher order isoparametric fictitious domain method for level set domains”. *Geometrically Unfitted Finite Element Methods and Applications* Lecture Notes in Computational Science and Engineering 121, Springer, 65–92, 2017. doi: [10.1007/978-3-319-71431-8_3](https://doi.org/10.1007/978-3-319-71431-8_3).

* C. Lehrenfeld. “High order unfitted finite element methods on level set domains using isoparametric mappings”. *Comput. Methods Appl. Mech. Engrg.* 300, 716–733, 2016. doi: [10.1016/j.cma.2015.12.005](https://doi.org/10.1016/j.cma.2015.12.005).

### Interface problems and two-phase flows

Nitsche-XFEM / CutFEM discretizations of scalar and Stokes interface problems, two-phase flows and the corresponding preconditioners.

* E. Burman and J. Preuß. “Unique continuation for an elliptic interface problem using unfitted isoparametric finite elements”. *SMAI J. Comput. Math.* 2025. doi: [10.5802/smai-jcm.122](https://doi.org/10.5802/smai-jcm.122), arXiv: [2307.05210](https://arxiv.org/abs/2307.05210). Code: [GitHub](https://github.com/UCL/interface-uc-unfitted-iso), [Zenodo (doi: 10.5281/zenodo.13328144)](https://doi.org/10.5281/zenodo.13328144).

* D. Capatina, F. Caubet, M. Dambrine and R. Zelada. “Nitsche extended finite element method of a Ventcel transmission problem with discontinuities at the interface”. *ESAIM Math. Model. Numer. Anal.* 59(2), 999–1021, 2025. doi: [10.1051/m2an/2025014](https://doi.org/10.1051/m2an/2025014).

* S. Groß and A. Reusken. “Analysis of optimal preconditioners for CutFEM”. *Numer. Linear Algebra Appl.* 30(5), e2486, 2023. doi: [10.1002/nla.2486](https://doi.org/10.1002/nla.2486), arXiv: [2202.09069](https://arxiv.org/abs/2202.09069). Code: [Zenodo (doi: 10.5281/zenodo.7249209)](https://doi.org/10.5281/zenodo.7249209).

* M. A. Olshanskii, A. Quaini and Q. Sun. “A finite element method for two-phase flow with material viscous interface”. *Comput. Methods Appl. Math.* 22(2), 443–464, 2022. doi: [10.1515/cmam-2021-0185](https://doi.org/10.1515/cmam-2021-0185), arXiv: [2106.02922](https://arxiv.org/abs/2106.02922).

* M. A. Olshanskii, A. Quaini and Q. Sun. “An unfitted finite element method for two-phase Stokes problems with slip between phases”. *J. Sci. Comput.* 89(2), 41, 2021. doi: [10.1007/s10915-021-01658-x](https://doi.org/10.1007/s10915-021-01658-x), arXiv: [2101.09627](https://arxiv.org/abs/2101.09627).

* T. Ludescher. “Multilevel preconditioning of stabilized unfitted finite element discretizations”. *PhD thesis* RWTH Aachen University, 2020. doi: [10.18154/RWTH-2020-07305](https://doi.org/10.18154/RWTH-2020-07305).

* C. Lehrenfeld and A. Reusken. “L2-error analysis of an isoparametric unfitted finite element method for elliptic interface problems”. *J. Numer. Math.* 27(2), 85–99, 2019. doi: [10.1515/jnma-2017-0109](https://doi.org/10.1515/jnma-2017-0109).

* C. Lehrenfeld and A. Reusken. “Analysis of a high-order unfitted finite element method for elliptic interface problems”. *IMA J. Numer. Anal.* 38(3), 1351–1387, 2018. doi: [10.1093/imanum/drx041](https://doi.org/10.1093/imanum/drx041), arXiv: [1602.02970](https://arxiv.org/abs/1602.02970).

* C. Lehrenfeld and A. Reusken. “Optimal preconditioners for Nitsche-XFEM discretizations of interface problems”. *Numer. Math.* 135(2), 313–332, 2017. doi: [10.1007/s00211-016-0801-6](https://doi.org/10.1007/s00211-016-0801-6), arXiv: [1408.2940](https://arxiv.org/abs/1408.2940).

* P. L. Lederer, C.-M. Pfeiler, C. Wintersteiger and C. Lehrenfeld. “Higher order unfitted FEM for Stokes interface problems”. *Proc. Appl. Math. Mech.* 16(1), 7–10, 2016. doi: [10.1002/pamm.201610003](https://doi.org/10.1002/pamm.201610003).

* C. Lehrenfeld. “High order unfitted finite element methods on level set domains using isoparametric mappings”. *Comput. Methods Appl. Mech. Engrg.* 300, 716–733, 2016. doi: [10.1016/j.cma.2015.12.005](https://doi.org/10.1016/j.cma.2015.12.005).

* C. Lehrenfeld. “Removing the stabilization parameter in fitted and unfitted symmetric Nitsche formulations”. *Proceedings of the VII European Congress on Computational Methods in Applied Sciences and Engineering (ECCOMAS Congress 2016)* 2016. doi: [10.7712/100016.1820.4573](https://doi.org/10.7712/100016.1820.4573), arXiv: [1603.00617](https://arxiv.org/abs/1603.00617).

### Moving domains: space-time methods

Unfitted space-time finite element methods for PDEs on moving domains.

* E. Burman and F. Heimann. “Higher order unfitted space-time methods for transport problems”. preprint, 2025. arXiv: [2509.02253](https://arxiv.org/abs/2509.02253). Code: [Zenodo (doi: 10.5281/zenodo.17185244)](https://doi.org/10.5281/zenodo.17185244).

* F. Heimann and C. Lehrenfeld. “A higher order unfitted space-time finite element method for coupled surface-bulk problems”. *Numerical Mathematics and Advanced Applications ENUMATH 2023* Lecture Notes in Computational Science and Engineering, Springer, 2025. doi: [10.1007/978-3-031-86173-4_43](https://doi.org/10.1007/978-3-031-86173-4_43), arXiv: [2401.07807](https://arxiv.org/abs/2401.07807). Code: [GitLab](https://gitlab.gwdg.de/fabian.heimann/repro-ho-unf-space-time-coupled).

* F. Heimann, C. Lehrenfeld and J. Preuß. “Discretization error analysis of a high-order unfitted space-time method for moving domain problems”. *IMA J. Numer. Anal.* 2025. doi: [10.1093/imanum/draf084](https://doi.org/10.1093/imanum/draf084), arXiv: [2504.08608](https://arxiv.org/abs/2504.08608). (advance article)

* F. Heimann and C. Lehrenfeld. “Geometry error analysis of a parametric mapping for higher order unfitted space-time methods”. *IMA J. Numer. Anal.* 45(6), 3643–3697, 2025. doi: [10.1093/imanum/drae098](https://doi.org/10.1093/imanum/drae098), arXiv: [2311.02348](https://arxiv.org/abs/2311.02348). Code: [GitLab](https://gitlab.gwdg.de/fabian.heimann/repro-ho-unf-space-time-fem-geom).

* F. Heimann. “Higher order unfitted space-time finite element methods for moving domain problems”. *PhD thesis* Georg-August-Universität Göttingen, 2025. doi: [10.53846/goediss-11003](https://doi.org/10.53846/goediss-11003).

* A. Reusken and H. Sass. “Analysis of a space-time unfitted finite element method for PDEs on evolving surfaces”. preprint, 2024. arXiv: [2401.01215](https://arxiv.org/abs/2401.01215).

* F. Heimann, C. Lehrenfeld and J. Preuß. “Geometrically higher order unfitted space-time methods for PDEs on moving domains”. *SIAM J. Sci. Comput.* 45(2), B139–B165, 2023. doi: [10.1137/22M1476034](https://doi.org/10.1137/22M1476034), arXiv: [2202.02216](https://arxiv.org/abs/2202.02216). Code: [GitLab](https://gitlab.gwdg.de/fabian.heimann/repro-ho-unf-space-time-fem).

* H. Sass. “Space-time trace finite element methods for partial differential equations on evolving surfaces”. *PhD thesis* RWTH Aachen University, 2022. doi: [10.18154/RWTH-2022-09895](https://doi.org/10.18154/RWTH-2022-09895).

* A. C. Wendler. “Monolithic unfitted space-time FEM for an osmotic cell swelling problem”. *Master's thesis* Georg-August-Universität Göttingen, 2022. doi: [10.25625/0KPEON](https://doi.org/10.25625/0KPEON).

* F. Heimann. “On discontinuous- and continuous-in-time unfitted space-time methods for PDEs on moving domains”. *Master's thesis* Georg-August-Universität Göttingen, 2020. doi: [10.25625/CDCMYT](https://doi.org/10.25625/CDCMYT).

* J. Preuß. “Higher order unfitted isoparametric space-time FEM on moving domains”. *Master's thesis* Georg-August-Universität Göttingen, 2018. doi: [10.25625/UACWXS](https://doi.org/10.25625/UACWXS).

### Moving domains: Eulerian time-stepping

Unfitted finite element methods on moving domains with time-stepping on a fixed background mesh (ghost-penalty extension, BDF schemes, narrow band methods).

* M. A. Olshanskii and H. von Wahl. “A conservative Eulerian finite element method for transport and diffusion in moving domains”. *Comput. Methods Appl. Math.* 25(4), 961–979, 2025. doi: [10.1515/cmam-2024-0055](https://doi.org/10.1515/cmam-2024-0055), arXiv: [2404.07130](https://arxiv.org/abs/2404.07130). Code: [GitHub](https://github.com/hvonwah/conserv-eulerian-moving-domian-repro), [Zenodo (doi: 10.5281/zenodo.10951768)](https://doi.org/10.5281/zenodo.10951768).

* M. A. Olshanskii, A. Reusken and P. Schwering. “A narrow band finite element method for the level set equation”. *SIAM J. Sci. Comput.* 47(2), 2025. doi: [10.1137/24M1674182](https://doi.org/10.1137/24M1674182), arXiv: [2407.02950](https://arxiv.org/abs/2407.02950).

* M. A. Olshanskii and H. von Wahl. “Stability of instantaneous pressures in an Eulerian finite element method for moving boundary flow problems”. *J. Comput. Phys.* 2025. arXiv: [2412.17657](https://arxiv.org/abs/2412.17657). Code: [GitHub](https://github.com/hvonwah/stable_inst_pressure_moving_domain_repro), [Zenodo (doi: 10.5281/zenodo.14548166)](https://doi.org/10.5281/zenodo.14548166).

* M. Neilan and M. A. Olshanskii. “An Eulerian finite element method for the linearized Navier–Stokes problem in an evolving domain”. *IMA J. Numer. Anal.* 44(6), 3234–3258, 2024. doi: [10.1093/imanum/drad105](https://doi.org/10.1093/imanum/drad105), arXiv: [2308.01444](https://arxiv.org/abs/2308.01444).

* H. von Wahl and T. Richter. “An Eulerian time-stepping scheme for a coupled parabolic moving domain problem using equal order unfitted finite elements”. *Proc. Appl. Math. Mech.* 22(1), e202200003, 2023. doi: [10.1002/pamm.202200003](https://doi.org/10.1002/pamm.202200003).

* H. von Wahl and T. Richter. “Error analysis for a parabolic PDE model problem on a coupled moving domain in a fully Eulerian framework”. *SIAM J. Numer. Anal.* 61(1), 286–314, 2023. doi: [10.1137/21M1458417](https://doi.org/10.1137/21M1458417), arXiv: [2111.05607](https://arxiv.org/abs/2111.05607). Code: [GitHub](https://github.com/hvonwah/repro-eulerian-coupled-heat), [Zenodo (doi: 10.5281/zenodo.6505243)](https://doi.org/10.5281/zenodo.6505243).

* Y. Lou and C. Lehrenfeld. “Isoparametric unfitted BDF–finite element method for PDEs on evolving domains”. *SIAM J. Numer. Anal.* 60(4), 2069–2098, 2022. doi: [10.1137/21M142126X](https://doi.org/10.1137/21M142126X), arXiv: [2105.09162](https://arxiv.org/abs/2105.09162). Code: [GitLab](https://gitlab.gwdg.de/lehrenfeld/repro-isop-unf-bdf-fem).

* H. von Wahl, T. Richter and C. Lehrenfeld. “An unfitted Eulerian finite element method for the time-dependent Stokes problem on moving domains”. *IMA J. Numer. Anal.* 42(3), 2505–2544, 2022. doi: [10.1093/imanum/drab044](https://doi.org/10.1093/imanum/drab044), arXiv: [2002.02352](https://arxiv.org/abs/2002.02352). Code: [Zenodo (doi: 10.5281/zenodo.3647571)](https://doi.org/10.5281/zenodo.3647571).

* C. Lehrenfeld and M. A. Olshanskii. “An Eulerian finite element method for PDEs in time-dependent domains”. *ESAIM Math. Model. Numer. Anal.* 53(2), 585–614, 2019. doi: [10.1051/m2an/2018068](https://doi.org/10.1051/m2an/2018068).

### Fluid-structure interaction, contact and fracture

Fully Eulerian unfitted methods for fluid-rigid body and fluid-structure interaction, including contact and phase-field fracture.

* S. Lee, H. von Wahl and T. Wick. “A thermo-flow-mechanics-fracture model coupling a phase-field interface approach and thermo-fluid-structure interaction”. *Int. J. Numer. Methods Eng.* 2025. doi: [10.1002/nme.7646](https://doi.org/10.1002/nme.7646), arXiv: [2409.03416](https://arxiv.org/abs/2409.03416). Code: [GitHub](https://github.com/hvonwah/repro-tfsi-pff), [Zenodo (doi: 10.5281/zenodo.13685486)](https://doi.org/10.5281/zenodo.13685486).

* H. von Wahl and T. Wick. “A coupled high-accuracy phase-field fluid-structure interaction framework for Stokes fluid-filled fracture surrounded by an elastic medium”. *Results Appl. Math.* 22, 100455, 2024. doi: [10.1016/j.rinam.2024.100455](https://doi.org/10.1016/j.rinam.2024.100455), arXiv: [2308.15400](https://arxiv.org/abs/2308.15400). Code: [GitHub](https://github.com/hvonwah/repro-coupled-phase-field-fsi), [Zenodo (doi: 10.5281/zenodo.10362611)](https://doi.org/10.5281/zenodo.10362611).

* H. von Wahl and T. Wick. “A high-accuracy framework for phase-field fracture interface reconstructions with application to Stokes fluid-filled fracture surrounded by an elastic medium”. *Comput. Methods Appl. Mech. Engrg.* 415, 116202, 2023. doi: [10.1016/j.cma.2023.116202](https://doi.org/10.1016/j.cma.2023.116202), arXiv: [2212.07982](https://arxiv.org/abs/2212.07982). Code: [GitHub](https://github.com/hvonwah/stationary_phase_field_stokes_fsi), [Zenodo (doi: 10.5281/zenodo.7443025)](https://doi.org/10.5281/zenodo.7443025).

* M. Kemper. “Pure Eulerian unfitted FEM for biological fluid-structure interaction problems”. *Master's thesis* Georg-August-Universität Göttingen, 2022. doi: [10.25625/DYUGCA](https://doi.org/10.25625/DYUGCA).

* H. von Wahl, T. Richter, S. Frei and T. Hagemeier. “Falling balls in a viscous fluid with contact: Comparing numerical simulations with experimental data”. *Phys. Fluids* 33(3), 033304, 2021. doi: [10.1063/5.0037971](https://doi.org/10.1063/5.0037971), arXiv: [2011.08691](https://arxiv.org/abs/2011.08691). Code: [Zenodo (doi: 10.5281/zenodo.3989604)](https://doi.org/10.5281/zenodo.3989604).

* H. von Wahl. “Unfitted finite elements for fluid-rigid body interaction problems”. *PhD thesis* Otto-von-Guericke-Universität Magdeburg, 2021. doi: [10.25673/40013](https://doi.org/10.25673/40013).

* H. von Wahl and T. Richter. “Using a deep neural network to predict the motion of under-resolved triangular rigid bodies in an incompressible flow”. *Int. J. Numer. Methods Fluids* 93(12), 2021. doi: [10.1002/fld.5037](https://doi.org/10.1002/fld.5037), arXiv: [2102.11636](https://arxiv.org/abs/2102.11636).

### PDEs on surfaces (TraceFEM)

Trace finite element methods for scalar and vector-valued PDEs on stationary and evolving surfaces.

* T. Alemán and A. Reusken. “Numerical analysis of a constrained strain energy minimization problem”. *SIAM J. Sci. Comput.* 48(3), A1564–A1586, 2026. doi: [10.1137/24M1713545](https://doi.org/10.1137/24M1713545), arXiv: [2411.19089](https://arxiv.org/abs/2411.19089).

* L. Bouck, R. H. Nochetto, M. Shakipov and V. Yushutin. “Inf-sup stability of parabolic TraceFEM”. *Found. Comput. Math.* 2026. doi: [10.1007/s10208-026-09757-7](https://doi.org/10.1007/s10208-026-09757-7), arXiv: [2409.13944](https://arxiv.org/abs/2409.13944).

* M. Neilan and H. Wan. “A TraceFEM C0 interior penalty method for the surface biharmonic equation”. *J. Sci. Comput.* 108(3), 66, 2026. doi: [10.1007/s10915-026-03384-8](https://doi.org/10.1007/s10915-026-03384-8), arXiv: [2512.18949](https://arxiv.org/abs/2512.18949).

* F. Heimann and C. Lehrenfeld. “A higher order unfitted space-time finite element method for coupled surface-bulk problems”. *Numerical Mathematics and Advanced Applications ENUMATH 2023* Lecture Notes in Computational Science and Engineering, Springer, 2025. doi: [10.1007/978-3-031-86173-4_43](https://doi.org/10.1007/978-3-031-86173-4_43), arXiv: [2401.07807](https://arxiv.org/abs/2401.07807). Code: [GitLab](https://gitlab.gwdg.de/fabian.heimann/repro-ho-unf-space-time-coupled).

* P. Schwering. “Surface Navier–Stokes equations: numerical methods and analysis”. *PhD thesis* RWTH Aachen University, 2025. Code: [Zenodo (doi: 10.5281/zenodo.17120084)](https://doi.org/10.5281/zenodo.17120084).

* E. Bachini, P. Brandner, T. Jankuhn, M. Nestler, S. Praetorius, A. Reusken and A. Voigt. “Diffusion of tangential tensor fields: numerical issues and influence of geometric properties”. *J. Numer. Math.* 32(1), 55–75, 2024. doi: [10.1515/jnma-2022-0088](https://doi.org/10.1515/jnma-2022-0088), arXiv: [2205.12581](https://arxiv.org/abs/2205.12581). Code: [Zenodo (doi: 10.5281/zenodo.7096487)](https://doi.org/10.5281/zenodo.7096487).

* M. A. Olshanskii, A. Reusken and P. Schwering. “An Eulerian finite element method for tangential Navier–Stokes equations on evolving surfaces”. *Math. Comp.* 93, 2031–2065, 2024. arXiv: [2302.00779](https://arxiv.org/abs/2302.00779).

* A. Reusken and H. Sass. “Analysis of a space-time unfitted finite element method for PDEs on evolving surfaces”. preprint, 2024. arXiv: [2401.01215](https://arxiv.org/abs/2401.01215).

* S. Lu and X. Xu. “Numerical investigations on trace finite element methods for the Laplace–Beltrami eigenvalue problem”. *J. Sci. Comput.* 97(1), 12, 2023. doi: [10.1007/s10915-023-02326-y](https://doi.org/10.1007/s10915-023-02326-y), arXiv: [2108.02434](https://arxiv.org/abs/2108.02434). Code: [GitHub](https://github.com/lusongno1/surface_PDEs).

* H. Sass and A. Reusken. “An accurate and robust Eulerian finite element method for partial differential equations on evolving surfaces”. *Comput. Math. Appl.* 146, 253–270, 2023. arXiv: [2212.12030](https://arxiv.org/abs/2212.12030).

* E. Schlesinger. “Embedded Trefftz trace DG methods for PDEs on unfitted surfaces”. *Master's thesis* Georg-August-Universität Göttingen, 2023. doi: [10.25625/QTOPWD](https://doi.org/10.25625/QTOPWD). Code: [GRO.data (doi: 10.25625/L6J3DM)](https://doi.org/10.25625/L6J3DM).

* P. Brandner, T. Jankuhn, S. Praetorius, A. Reusken and A. Voigt. “Finite element discretization methods for velocity-pressure and stream function formulations of surface Stokes equations”. *SIAM J. Sci. Comput.* 44(4), A1807–A1832, 2022. doi: [10.1137/21M1403126](https://doi.org/10.1137/21M1403126), arXiv: [2103.03843](https://arxiv.org/abs/2103.03843).

* P. Brandner. “Numerical methods for surface Navier–Stokes equations in stream function formulation”. *PhD thesis* RWTH Aachen University, 2022. doi: [10.18154/RWTH-2022-04531](https://doi.org/10.18154/RWTH-2022-04531).

* A. Reusken. “Analysis of finite element methods for surface vector-Laplace eigenproblems”. *Math. Comp.* 91, 1587–1623, 2022. doi: [10.1090/mcom/3728](https://doi.org/10.1090/mcom/3728), arXiv: [2011.02851](https://arxiv.org/abs/2011.02851).

* H. Sass. “Space-time trace finite element methods for partial differential equations on evolving surfaces”. *PhD thesis* RWTH Aachen University, 2022. doi: [10.18154/RWTH-2022-09895](https://doi.org/10.18154/RWTH-2022-09895).

* T. Jankuhn and A. Reusken. “Trace finite element methods for surface vector-Laplace equations”. *IMA J. Numer. Anal.* 41(1), 48–83, 2021. doi: [10.1093/imanum/drz062](https://doi.org/10.1093/imanum/drz062), arXiv: [1904.12494](https://arxiv.org/abs/1904.12494).

* P. Brandner and A. Reusken. “Finite element error analysis of surface Stokes equations in stream function formulation”. *ESAIM Math. Model. Numer. Anal.* 54(6), 2069–2097, 2020. doi: [10.1051/m2an/2020044](https://doi.org/10.1051/m2an/2020044).

* T. Jankuhn and A. Reusken. “Higher order trace finite element methods for the surface Stokes equation”. preprint, 2019. arXiv: [1909.08327](https://arxiv.org/abs/1909.08327).

* J. Grande, C. Lehrenfeld and A. Reusken. “Analysis of a high-order trace finite element method for PDEs on level set surfaces”. *SIAM J. Numer. Anal.* 56(1), 228–255, 2018. doi: [10.1137/16M1102203](https://doi.org/10.1137/16M1102203).

* F. Heimann. “Higher order discontinuous Galerkin methods for the Laplace-Beltrami problem on unfitted smooth surfaces”. *Bachelor's thesis* Georg-August-Universität Göttingen, 2018. doi: [10.25625/OIBRT4](https://doi.org/10.25625/OIBRT4).

### Optimization, inverse problems and model order reduction

Shape optimization, PDE-constrained optimal control, unique continuation / data assimilation and reduced order modelling with unfitted discretizations.

* E. Burman, J. Preuß and T. van Beeck. “Variational data assimilation for the wave equation in heterogeneous media: numerical investigation of stability”. *Commun. Appl. Math. Comput.* 2026. doi: [10.1007/s42967-026-00583-w](https://doi.org/10.1007/s42967-026-00583-w), arXiv: [2509.13108](https://arxiv.org/abs/2509.13108). Code: [GitHub](https://github.com/TimvanBeeck/waveUC_discCoefs).

* E. Burman and J. Preuß. “Unique continuation for an elliptic interface problem using unfitted isoparametric finite elements”. *SMAI J. Comput. Math.* 2025. doi: [10.5802/smai-jcm.122](https://doi.org/10.5802/smai-jcm.122), arXiv: [2307.05210](https://arxiv.org/abs/2307.05210). Code: [GitHub](https://github.com/UCL/interface-uc-unfitted-iso), [Zenodo (doi: 10.5281/zenodo.13328144)](https://doi.org/10.5281/zenodo.13328144).

* E. Burman, L. Oksanen, J. Preuß and Z. Zhao. “Unique continuation for the wave equation: the stability landscape”. preprint, 2025. Code: [GitHub](https://github.com/janoschpreuss/wave-uc-stability-landscape-repro), [Zenodo (doi: 10.5281/zenodo.17370507)](https://doi.org/10.5281/zenodo.17370507).

* W. Gong and Z. Zhang. “A novel shape optimization approach for source identification in elliptic equations”. *Inverse Probl. Imaging* 2025. doi: [10.3934/ipi.2025037](https://doi.org/10.3934/ipi.2025037), arXiv: [2407.02909](https://arxiv.org/abs/2407.02909).

* G. Katsouleas, E. N. Karatzas and F. Travlopanos. “Discrete empirical interpolation and unfitted mesh FEMs: application in PDE-constrained optimization”. *Optimization* 72(6), 1609–1642, 2023. doi: [10.1080/02331934.2022.2032697](https://doi.org/10.1080/02331934.2022.2032697), arXiv: [2010.09059](https://arxiv.org/abs/2010.09059).

* A. Aretaki and E. N. Karatzas. “Random geometries for optimal control PDE problems based on fictitious domain FEMs and cut elements”. *J. Comput. Appl. Math.* 412, 114286, 2022. doi: [10.1016/j.cam.2022.114286](https://doi.org/10.1016/j.cam.2022.114286), arXiv: [2003.00352](https://arxiv.org/abs/2003.00352).

* E. N. Karatzas, M. Nonino, F. Ballarin and G. Rozza. “A reduced order cut finite element method for geometrically parametrized steady and unsteady Navier–Stokes problems”. *Comput. Math. Appl.* 116, 140–160, 2022. doi: [10.1016/j.camwa.2021.07.016](https://doi.org/10.1016/j.camwa.2021.07.016), arXiv: [2010.04953](https://arxiv.org/abs/2010.04953).

* E. N. Karatzas and G. Rozza. “A reduced order model for a stable embedded boundary parametrized Cahn–Hilliard phase-field system based on cut finite elements”. *J. Sci. Comput.* 89, 9, 2021. doi: [10.1007/s10915-021-01623-8](https://doi.org/10.1007/s10915-021-01623-8), arXiv: [2009.01596](https://arxiv.org/abs/2009.01596).

* E. N. Karatzas, F. Ballarin and G. Rozza. “Projection-based reduced order models for a cut finite element method in parametrized domains”. *Comput. Math. Appl.* 79(3), 833–851, 2020. doi: [10.1016/j.camwa.2019.08.003](https://doi.org/10.1016/j.camwa.2019.08.003), arXiv: [1901.03846](https://arxiv.org/abs/1901.03846).

* H.-G. Raumer. “Shape optimization for interface problems using unfitted finite elements”. *Master's thesis* Georg-August-Universität Göttingen, 2018. url: <http://cpde.math.uni-goettingen.de/data/Rau18_Ma.pdf>.

### Further applications and benchmarks

Works in which ngsxfem is used for other applications, for fitted space-time discretizations or as one of several benchmarked codes.

* M. Loibl, G. H. Teixeira, T. Toprak, I. Shishkina, C. Miao, J. Kiendl, F. Kummer and B. Marussig. “Comparative study of different quadrature methods for cut elements”. *Arch. Comput. Methods Eng.* 2026. doi: [10.1007/s11831-026-10619-2](https://doi.org/10.1007/s11831-026-10619-2). Code: [GitHub](https://github.com/B2-M/CutElementIntegration), [Zenodo (doi: 10.5281/zenodo.14961680)](https://doi.org/10.5281/zenodo.14961680). (ngsxfem is one of the benchmarked cut-element quadrature codes)

* T. Toprak, M. Loibl, G. H. Teixeira, I. Shishkina, C. Miao, J. Kiendl, B. Marussig and F. Kummer. “Employing continuous integration inspired workflows for benchmarking of scientific software – a use case on numerical cut element quadrature”. *Adv. Eng. Softw.* 213, 104087, 2026. doi: [10.1016/j.advengsoft.2025.104087](https://doi.org/10.1016/j.advengsoft.2025.104087), arXiv: [2503.17192](https://arxiv.org/abs/2503.17192). Code: [GitHub](https://github.com/B2-M/CutElementIntegration). (ngsxfem is one of the benchmarked cut-element quadrature codes)

* G. Fu and Y. Yang. “A hybridizable discontinuous Galerkin method on unfitted meshes for single-phase Darcy flow in fractured porous media”. *Adv. Water Resour.* 173, 104390, 2023. doi: [10.1016/j.advwatres.2023.104390](https://doi.org/10.1016/j.advwatres.2023.104390), arXiv: [2209.05445](https://arxiv.org/abs/2209.05445). Code: [GitHub](https://github.com/gridfunction/fracturedPorousMedia).

* J. Smith-Roberge. “Microcolony dynamics: motion from growth, order, and incompressibility”. *PhD thesis* University of Waterloo, 2023. hdl: [10012/19340](https://hdl.handle.net/10012/19340).

* G. Fu and Z. Xu. “High-order space-time finite element methods for the Poisson–Nernst–Planck equations: positivity and unconditional energy stability”. *Comput. Methods Appl. Mech. Engrg.* 395, 115031, 2022. doi: [10.1016/j.cma.2022.115031](https://doi.org/10.1016/j.cma.2022.115031), arXiv: [2105.01163](https://arxiv.org/abs/2105.01163). (fitted space-time discretization using the space-time tools of ngsxfem)

### Software building on ngsxfem

Libraries and solver collections that build on top of ngsxfem.

* ngsxfem developers. “ngsxditto: a high-level Python library for PDEs on moving domains and two-phase flows built on NGSolve and ngsxfem”. software, 2025. Code: [GitHub](https://github.com/ngsxfem/ngsxditto), [Documentation](https://ngsxfem.github.io/ngsxditto/).

* M. Shakipov. “surf-pde: a collection of TraceFEM solvers for surface PDEs implemented with ngsxfem”. software, 2023. Code: [GitHub](https://github.com/chromomons/surf-pde).

* P. Stocker. “NGSTrefftz: Add-on to NGSolve for Trefftz methods”. *Journal of Open Source Software* 7(71), 4135, 2022. doi: [10.21105/joss.04135](https://doi.org/10.21105/joss.04135). Code: [GitHub](https://github.com/PaulSt/NGSTrefftz). (NGSTrefftz uses ngsxfem for unfitted (embedded) Trefftz methods)
