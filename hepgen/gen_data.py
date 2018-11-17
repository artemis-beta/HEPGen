from hepgen.decays import decay_table 
from scipy.stats import expon

def run_gen(decay_id):
    _decay = decay_table(decay_id)
    assert _decay
    _res_dict = {'{}_TAU'.format(_decay.Mother.name) : expon.rvs(scale=_decay.Mother.lifetime)}
    for daughter in _decay.Daughters:
        _res_dict['{}_TAU'.format(daughter.name)] = expon.rvs(scale=daughter.lifetime)

    print(_res_dict)


if __name__ in "__main__":
    run_gen('SK020001')
