import os
from random import randint
from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import random_walk, copolymer
from pysimm.apps.equilibrate import equil

class TgWorkflow(object):
    def __init__(self, monomer, chain_length, ff, growth_density=0.3, 
            equil_temp=600, cool_step_time=100000,
            cool_output=None, random_walk_sim=None, workdir=None, 
            equil_output=None, nproc=1, cool_temp_range=None):
        self.monomer = monomer
        self.polymer = None
        self.chain_length = chain_length
        self.ff = ff
        self.growth_density = growth_density
        self.equil_temp = equil_temp
        self.cool_step_time = cool_step_time
        self.workdir = workdir
        if self.workdir is None:
            self.workdir = os.getcwd()
        self.equil_output = equil_output
        if self.equil_output is None:
            self.equil_output = lmps.OutputSettings(
                dump={'freq': 10000, 'filename': 'dump.equil.*'}, 
                thermo={'freq': 1000})
        self.cool_output = cool_output
        if self.cool_output is None:
            self.cool_output = lmps.OutputSettings(
                dump={'freq': 10000, 'filename': 'dump.cool.*'}, 
                thermo={'freq': 1000})
        self.cool_temp_range = cool_temp_range
        if self.cool_temp_range is None:
            self.cool_temp_range = reversed(range(100, 601, 10))
        self.nproc = self.nproc
    
    def prepare_monomer(self):
        # preparation for monomer before simulation
        pass
    
    def grow_polymer(self):
        polymer = random_walk(self.monomer, self.chain_length, 
            forcefield=self.ff, density=self.grow_polymer, 
            settings={'np': self.nproc})
        
    def cap_polymer(self):
        # cap polymer chain ends
        pass
        
    def equilibrate(self):
        lmps.quick_min(self.polymer)
        self.polymer = equil(polymer, tfinal=self.equil_temp, np=self.nproc,
            output_settings=self.equil_output, log='log.equil')
        
    def equilibrate_npt(self):
        sim = lmps.Simulation(self.polymer)
        sim.add(lmps.OutputSettings(dump={'freq': 10000, 'filename': 'dump.equil.npt.*'}, thermo={'freq': 1000}))
        sim.add_md(ensemble='npt', temperature=self.equil_temp, length=1000000)
        sim.run(np=self.nproc)
        
    def stepwise_cooling(self):
        sim = lmps.Simulation(self.polymer, name='cool', log='log.cool')
        sim.add(self.cool_output)
        for temp in self.cool_temp_range:
            velocity = lmps.Velocity(style='scale', temperature=temp)
            md = lmps.MolecularDynamics(ensemble='npt', temperature=temp, pressure=1., run=self.cool_step_time, timestep=1)
            sim.add(velocity)
            sim.add(md)
        sim.run(np=nproc)
        
    def run():
        os.chdir(self.workdir)
        self.prepare_monomer()
        self.grow_polymer()
        self.cap_polymer()
        self.equilibrate()
        self.equilibrate_npt()
        self.stepwise_cooling()