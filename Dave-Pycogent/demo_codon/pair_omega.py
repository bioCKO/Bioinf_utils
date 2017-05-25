from warnings import filterwarnings

filterwarnings("ignore", "Not using MPI as mpi4py not found")

from cogent import LoadSeqs, LoadTree
from cogent.evolve.models import CNFGTR
from cogent.maths.stats import chisqprob
import sys


# change the filename
aln = LoadSeqs(sys.argv[1])

assert len(aln.Names) <= 3, "Need to modify the script to handle case of more than 3 taxa"

tree = LoadTree(tip_names=aln.Names[:])

sm = CNFGTR()
lf = sm.makeLikelihoodFunction(tree)
lf.setAlignment(aln)

# following required for pairs
lf.setParamRule('length', is_independent=True)

# we test neutrality by setting omega to 1 first off

lf.setParamRule('omega', value=1.0, is_constant=True)

# limit_action='raise' will cause an error if the functiuon cannot be optimise
# you can set that to 'warn', but I suggest raise
opt_args = dict(local=True, max_restarts=5, max_evaluations=100000, limit_action='raise') 
lf.optimise(**opt_args)
null_lnL = lf.getLogLikelihood()
null_nfp = lf.getNumFreeParams()

# now do alternate
lf.setParamRule('omega', is_constant=False)
lf.optimise(**opt_args)
alt_lnL = lf.getLogLikelihood()
alt_nfp = lf.getNumFreeParams()

# get the results tables
stats = lf.getStatistics(with_titles=True)
for stats_table in stats:
    print stats_table

# demo just getting out omega
omega = lf.getParamValue('omega')

# compute the likelihood ratio
lr = 2 * (alt_lnL - null_lnL)
df = alt_nfp - null_nfp
p = chisqprob(lr, df)

print 'LR=%.4f; df = %d; p=%.4f; omega=%.4f' % (lr, df, p, omega)
