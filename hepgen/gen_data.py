from hepgen.decays import decay_table 
from hepgen.data_tree import data_tree
from scipy.stats import expon

def run_gen(decay_id, tree=None, nevts=1):
    _tree  = data_tree('DecayTree_{}'.format(decay_id) if not tree else tree)
    _decay = decay_table(decay_id)
    assert _decay
    _tree.add_branch('{}_TAU'.format(_decay.Mother.name), 
                     expon.rvs(size=nevts, scale=_decay.Mother.lifetime))
    for daughter in _decay.Daughters:
        _tree.add_branch('{}_TAU'.format(daughter.name),
                     expon.rvs(size=nevts, scale=daughter.lifetime))

    return _tree
