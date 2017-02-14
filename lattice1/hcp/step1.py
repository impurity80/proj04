
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

name = 'dhcp'

comp = 0

curr_dir = os.getcwd()
os.system('mkdir eos')
os.system('mkdir result')

result = '{0}/result/result-{1}-{2}.txt'.format(curr_dir,name,comp)
os.system('rm {0}'.format(result))

result_sum = '{0}/result/summary-{1}.csv'.format(curr_dir,name)
os.system('rm {0}'.format(result_sum))

save(result, '{0} {1} '.format(name, comp))
save(result_sum, '{0} {1}'.format(name, comp))

RATIO_OPT = np.linspace(1.59, np.sqrt(8/3.0), 1)
ratio_energies = []

for ratio in RATIO_OPT:

    save(result_sum, '{0} {1}'.format(name, ratio))

    OPTIONS = np.linspace(0.97, 1.01, 9)
    volumes = []
    energies = []

    for opt in OPTIONS:
        #  a0 = 3.59 * opt / np.sqrt(2)
        a0 = 3.21 * opt
        c0 = ratio * a0
        # c0 = a0 * 1.585
        atoms = bulk('Mg', 'hcp', a=a0, c=c0)
        atoms.set_tags([1, 1])
        atoms = atoms*(1,1,2)

        alloys = []
        alloys.append(Alloy(1, 'Mg', 1.0, .0))

        calc = EMTO()
        calc.set(dir='work/{0}/comp-{1}/opt-{2}'.format(name, comp, opt),
                 lat=4,
                 ncpa=20,
                 amix=0.05,
                 afm='F',
                 kpts=[13, 13, 13]
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
        eos.plot('eos/{0}-{1}.png'.format(name,ratio))

        save(result, '{0} {1} {2} {3}'.format(v0, e0, B, (4.0 * v0) ** (1.0 / 3.0)))

    ratio_energies.append(e0)

    save(result, OPTIONS)
    save(result, volumes)
    save(result, energies)

    save(result, '------------------------')

    save(result_sum, '{0}, {1}, {2}, {3}, {4}, {5}'.format(ratio, e0, v0, B, volumes, energies))

eos = EquationOfState(RATIO_OPT, ratio_energies)
min_v, min_e, min_B = eos.fit()
eos.plot('eos/final.png')

save(result_sum, '{0}, {1}, {2}'.format(RATIO_OPT, ratio_energies, min_e))






