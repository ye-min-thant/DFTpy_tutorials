from ase.calculators.interface import Calculator
from ase.lattice.cubic import FaceCenteredCubic
from ase.optimize import BFGS, LBFGS, FIRE
from ase.optimize.sciopt import SciPyFminBFGS, SciPyFminCG
from ase.optimize import FIRE
from ase.constraints import StrainFilter, UnitCellFilter
from ase.io.trajectory import Trajectory
from ase import units
import ase.io

from dftpy.config.config import DefaultOption, PrintConf, OptionFormat
from dftpy.api.api4ase import DFTpyCalculator
import pathlib
dftpy_data_path = pathlib.Path(__file__).resolve().parent
logfile = dftpy_data_path/ "opt_ase.log"
energy_logfile = dftpy_data_path/"final_energy.log"
############################## initial config ##############################
conf = DefaultOption()
conf['PATH']['pppath'] = dftpy_data_path
conf['PP']['Al'] = 'al.gga.recpot'
conf['JOB']['calctype'] = 'Energy Force Stress'
conf['OPT']['method'] = 'CG-HS'
conf['GRID']['ecut'] = '1400'
conf['KEDF']['kedf'] = 'WT'
conf["EXC"]["xc"] = "PBE"
conf['OUTPUT']['stress'] = True
conf = OptionFormat(conf)
PrintConf(conf)
#-----------------------------------------------------------------------
path = conf['PATH']['pppath']
atoms = ase.io.read(path+'/'+'Al.vasp')
trajfile = 'opt.traj'

calc = DFTpyCalculator(config = conf)
atoms.set_calculator(calc)

############################## Relaxation type ##############################
'''
Ref :
    https ://wiki.fysik.dtu.dk/ase/ase/optimize.html#module-optimize
    https ://wiki.fysik.dtu.dk/ase/ase/constraints.html
'''
af = atoms
af = StrainFilter(atoms)
af = UnitCellFilter(atoms)
############################## Relaxation method ##############################
# opt = BFGS(af, trajectory = trajfile)
# opt = LBFGS(af, trajectory = trajfile, memory = 10, use_line_search = True)
# opt = LBFGS(af, trajectory = trajfile, memory = 10, use_line_search = False)
# opt = SciPyFminCG(af, trajectory = trajfile)
# opt = SciPyFminBFGS(af, trajectory = trajfile)
opt = FIRE(af, logfile=logfile, trajectory = trajfile)
opt.run(fmax = 1e-2, steps = 1000)

traj = Trajectory(trajfile)
ase.io.write('opt.vasp', traj[-1], direct = True, long_format=True, vasp5 = True)

energy_eV = atoms.get_potential_energy()
print(f"Final energy (eV): {energy_eV:.6f}\n")

with open(energy_logfile, "w") as f:
    f.write(f"Final energy (eV): {energy_eV:.6f}\n")
