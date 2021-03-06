import os
import shutil
from random import randint
from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import random_walk, copolymer
from pysimm.apps.equilibrate import equil

class TgWorkflow(object):
    def __init__(self, monomer, chain_length, ff, growth_density=0.3, 
            equil_temp=600, cool_step_time=100000,
            cool_output=None, random_walk_sim=None, workdir=None, 
            equil_step_length=100000, equil_npt_length=1000000,
            equil_output=None, nproc=1, cool_temp_range=None, calc_voronoi=True):
        self.monomer = monomer
        self.polymer = None
        self.chain_length = chain_length
        self.ff = ff
        self.growth_density = growth_density
        self.equil_temp = equil_temp
        self.equil_step_length = equil_step_length
        self.equil_npt_length = equil_npt_length
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
            args = [
                'id', 'type', 'mol', 'x', 'y', 'z', 
                'vx', 'vy', 'vz'
            ]
            if self.calc_voronoi:
                args += ['c_voronoi[1]']
            self.cool_output = lmps.OutputSettings(
                dump={
                    'freq': 10000, 'filename': 'dump.cool.*',
                    'args': args
                }, 
                thermo={'freq': 1000})
        self.cool_temp_range = cool_temp_range
        if self.cool_temp_range is None:
            self.cool_temp_range = reversed(range(100, 601, 10))
        self.nproc = nproc
        self.calc_voronoi = calc_voronoi
    
    def prepare_monomer(self):
        # preparation for monomer before simulation
        pass
    
    def grow_polymer(self):
        self.polymer = random_walk(self.monomer, self.chain_length, 
            forcefield=self.ff, density=self.growth_density, 
            settings={'np': self.nproc})
        self.polymer.forcefield = self.ff.name
        
    def growth_cleanup(self):
        # cap polymer chain ends
        pass
        
    def equilibrate(self):
        lmps.quick_min(self.polymer)
        equil(self.polymer, tfinal=self.equil_temp, np=self.nproc,
            output_settings=self.equil_output, log='log.equil',
            length_list=[self.equil_step_length for _ in range(7)])
        
    def equilibrate_npt(self):
        sim = lmps.Simulation(self.polymer)
        sim.add(lmps.OutputSettings(dump={'freq': 10000, 'filename': 'dump.equil.npt.*'}, thermo={'freq': 1000}))
        sim.add_md(ensemble='npt', temperature=self.equil_temp, length=self.equil_npt_length)
        sim.run(np=self.nproc)
        
    def stepwise_cooling(self):
        sim = lmps.Simulation(self.polymer, name='cool', log='log.cool')
        if self.calc_voronoi:
            sim.add('compute voronoi all voronoi/atom')
        sim.add(self.cool_output)
        for temp in self.cool_temp_range:
            velocity = lmps.Velocity(style='scale', temperature=temp)
            md = lmps.MolecularDynamics(ensemble='npt', temperature=temp, pressure=1., run=self.cool_step_time, timestep=1)
            sim.add(velocity)
            sim.add(md)
        sim.run(np=self.nproc)
        
    def run(self):
        if not os.path.exists(self.workdir):
            os.mkdir(self.workdir)
        os.chdir(self.workdir)
        self.prepare_monomer()
        self.grow_polymer()
        self.growth_cleanup()
        self.equilibrate()
        self.equilibrate_npt()
        self.stepwise_cooling()