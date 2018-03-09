<h1>UppASD</h1>

<b>Upp</b>sala <b>A</b>tomistic <b>S</b>pin <b>D</b>ynamics software
<!--![logo][logo]-->


The `UppASD` software package is a simulation suite to study magnetization dynamics by means of the atomistic version of the Landau-Lifshitz-Gilbert (LLG) equation.

***Applications:***
- [Critical temperature][11] determination
- [Magnon spectra][7] simulations
- [Micromagnetic exchange and DMI stiffness][6] calculations
- [Current][9], [field][8], and [temperature-driven][5] dynamics

***Features:***
- [LLG][12], Metropolis, and Heat-bath algorithms
- Generalized Hamiltonian including 
	-	bilinear and [biquadratic Heisenberg][10] exchange
	-	[Dzyaloshinskii-Moriya][15] exchange 
	-	magnetocrystalline anisotropies.
- Support for [arbitrary crystal structures][14] and [random alloys][13]
- [Spin-transfer][9] and spin-orbit torques
- [Temperature gradients][4]
- [Quantum statistics][3]

Detailed information about the method can be found in   
[**Atomistic Spin Dynamics: Foundations and Applications**  
O. Eriksson et. al,  Oxford University Press 2017][1]


---
Quick-start:
```python
./setup_UppASD.sh   
make <profile>
```

Where `<profile>` is a suitable compiler profile, i.e. `ifort`, `gfortran`, `gfortran-osx` and so on.   

The binary is compiled as `./source/sd`

Examples are provided in `./examples_revision_controlled/`

The manual is found at `./docs/UppASDmanual.pdf` 

A GUI for VTK-based visualization of output data is available at `./ASD_Viz/ASD_Visualizer.py`

---
(C) 2008-2018 [UppASD group][2]

[1]:https://global.oup.com/academic/product/atomistic-spin-dynamics-9780198788669
[2]:http://www.physics.uu.se/research/materials-theory/ongoing-research/uppasd/
[3]:https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.2.013802
[4]:https://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.014434
[5]:https://www.nature.com/articles/ncomms12430
[6]:https://journals.aps.org/prb/abstract/10.1103/PhysRevB.92.214424
[7]:http://iopscience.iop.org/article/10.1088/0953-8984/27/24/243202/meta
[8]:https://journals.aps.org/prb/abstract/10.1103/PhysRevB.86.224401
[9]:https://www.nature.com/articles/srep25685
[10]:https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.111.127204
[11]:https://journals.aps.org/prb/abstract/10.1103/PhysRevB.93.214439
[12]:http://iopscience.iop.org/article/10.1088/0953-8984/20/31/315203
[13]:https://journals.aps.org/prb/abstract/10.1103/PhysRevB.94.214410
[14]:https://journals.aps.org/prb/abstract/10.1103/PhysRevB.92.094411
[15]:https://www.nature.com/articles/ncomms5815
[logo]:./docs/uppasd_rot.png
