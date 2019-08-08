'''
Dimerization model:
 A + A <> AA
'''
from pysb import Model, Parameter, Monomer, Rule, Observable, Initial
import numpy as np
from swarm_it import SwarmParam

swarm_param = SwarmParam()

Model()
#######
V = 100.
#######
swarm_param(Parameter('kf',   0.001))
swarm_param(Parameter('kr',   1.), loc=np.log10(1.)-1., width=2.)


Monomer('A', ['d'])

# Rules
Rule('ReversibleBinding', A(d=None) + A(d=None) | A(d=1) % A(d=1), kf, kr)

#Observables
Observable("A_free", A(d=None))
Observable("A_dimer", A(d=1) % A(d=1))

# Inital Conditions
Parameter("A_0", 20.*V)
Initial(A(d=None), A_0)
