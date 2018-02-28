Computational Polymer Tg Workflow Implemented using pysimm
==========================================================

This code is a high-level representation of a computational workflow that can be used to calculate the glass transition of polymer materials using molecular simulations. The workflow covers each stage of the methdology laid out in Chapter 4 of the dissertation of Michael Fortunato from the Chemistry Department at the University of Florida.

Usage
=====

The workflow uses many features present in the pysimm software package. Details on installation and configuration with LAMMPS should can be found in its documentation.

The following code is used to create and run an simulation workflow to study polymethyl methacrylate. All parameters exposed to the user are defined explicitly below, however only the monomer, chain_length, and force field are required; the remaining keyword arguments have default values defined.

```python
from pysimm import forcefield
from pysimm.models.monomers.gaff2.pmma import monomer

from tg.workflow import TgWorkflow

pmma = monomer()
gaff = forcefield.Gaff2()
gaff.assign_charges(pmma, charges='gasteiger')
workflow = TgWorkflow(
    monomer=pmma, chain_length=700, ff=gaff, nproc=32, workdir=None,
    growth_density=0.3, equil_temp=600, cool_step_time=100000, 
    cool_output=None, random_walk_sim=None, equil_step_length=100000, 
    equil_npt_length=1000000, equil_output=None, cool_temp_range=None, 
    calc_voronoi=True
)
workflow.run()
```

The simulations defined by the workflow (polymer growth, equilibration, and stepwise cooling) should take ~5 hours on 32 processors.

The following code will parse the output files produced by those simulations to extract a glass transition temperature.

```python
from tg.fit import Tg
tg = Tg(
    logfiles=['log.cool'], temperature_rangenp.arange(600, 99, -10),
    split='auto1d'
)
print(tg.fit.tg, tg.fit.tg_density)
```