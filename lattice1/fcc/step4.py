
import csv
import os
from ase import Atom, Atoms
from ase.lattice.cubic import *
from ase.visualize import view
from numpy import *
from emto import *
from ase.utils.eos import EquationOfState
import matplotlib.pyplot as plt
from ase.lattice import bulk
from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

print rank, size

name = 'fccx4-ortho'

curr_dir = os.getcwd()
os.system('mkdir eos')
os.system('mkdir result')

result = '{0}/result/result-{1}.txt'.format(curr_dir,name)
os.system('rm {0}'.format(result))

result_sum = '{0}/result/summary-{1}.csv'.format(curr_dir,name)
os.system('rm {0}'.format(result_sum))

save(result, '{0}'.format(name))
save(result_sum, '{0}'.format(name))


OPTIONS = np.linspace(0.98, 1.02, 1)
volumes = []
energies = []

mn = 0.2
fe = 1.0-mn

for opt in OPTIONS:

    l = 3.59 * opt
    atoms = Atoms('Fe4',
                 scaled_positions=[
                                (0.5, 0.5, 0),
                                (0.5, 0, 0.5),
                                (0, 0.5, 0.5),
                                (0, 0, 0)],
                 magmoms = [5,5,5,5],
                 cell=[l, l, l],
                 pbc=(1,1,1))

    atoms.set_tags([1,1,1,1])

    alloys = []
    alloys.append(Alloy(1, 'Fe', fe / 2, 1.0))
    alloys.append(Alloy(1, 'Fe', fe / 2, -1.0))
    alloys.append(Alloy(1, 'Mn', mn / 2, 1.0))
    alloys.append(Alloy(1, 'Mn', mn / 2, -1.0))

    calc = EMTO()
    calc.set(dir='work/{0}/opt-{1}'.format(name, opt),
             lat=8,
             ncpa=20,
             amix=0.05,
             afm='F',
             kpts=[13,13,13]
             )
    calc.set_alloys(alloys)

    atoms.set_calculator(calc)
    nm_e = atoms.get_potential_energy() / atoms.get_number_of_atoms()
    nm_v = atoms.get_volume() / atoms.get_number_of_atoms()

    if nm_e < -0.001:
        volumes.append(nm_v)
        energies.append(nm_e)

    save(result, '{3} result : {0} {1} {2}'.format(opt, nm_v, nm_e, name))

print volumes, energies

temp_volumes = []
temp_energies = []
pivot = energies[0]
for v, e in zip(volumes, energies):
    if e - pivot > -0.06 and e - pivot < 0.01:
        temp_volumes.append(v)
        temp_energies.append(e)

if len(temp_volumes) > 1:
    eos = EquationOfState(temp_volumes, temp_energies)
    v0, e0, B = eos.fit()
    eos.plot('eos/{0}.png'.format(name))

    save(result, '{0} {1} {2} {3}'.format(v0, e0, B, (4.0 * v0) ** (1.0 / 3.0)))

save(result, OPTIONS)
save(result, volumes)
save(result, energies)

save(result, '------------------------')

save(result_sum, '{0}, {1}, {2}, {3}, {4}, {5}'.format(name, e0, v0, B, volumes, energies))







