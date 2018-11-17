from hepgen.decays import decay_table 
from hepgen.data_tree import data_tree
from hepgen import __version__
from pktools.PKLorentzVector import PKLorentzVector
from scipy.stats import expon
from random import uniform
from math import pi, atan

import logging
logging.basicConfig()

def sq_rt(number):
    _num = uniform(-1,1)
    _vec = _num/abs(_num)
    return _vec*pow(number, 0.5)

class HEPGen(object):
    def __init__(self, decay_id, tree=None, nevts=1):
        self._logger  = logging.getLogger(__class__.__name__)
        self._logger.setLevel('INFO')
        self._init()
        self._nevts   = nevts
        self._dec     = self._get_decay(decay_id)
        self._tree    = self._prepare_tree(tree)

    def _init(self):
        self._logger.info("\tRunning {}/v{}".format(__class__.__name__, __version__))

    def _get_decay(self, decid):
        try:
            _dec = decay_table(decid)
            assert _dec
            return _dec
        except AssertionError:
            self._logger.error("\tCould not find Decay '{}' in DecayTable".format(decid))
            exit(1)

    def _prepare_tree(self, treename):
        self._logger.info("\tCreating new data tree '{}'".format(treename))
        _tree = data_tree('DecayTree_{}'.format(self._dec.ID) if not treename else treename)
        for var in ['TAU', 'PX', 'THETA', 'PHI', 'P', 'PE', 'PY', 'PZ', 'PT']:
            _tree.add_branch('{}_{}'.format(self._dec.Mother.name, var))
            for daughter in self._dec.Daughters:
                _tree.add_branch('{}_{}'.format(daughter.name, var))
        return _tree

    def _gen_time(self):
        self._tree._fill_branch('{}_TAU'.format(self._dec.Mother.name), 
                               expon.rvs(scale=self._dec.Mother.lifetime))
        for daughter in self._dec.Daughters:
            self._tree._fill_branch('{}_TAU'.format(daughter.name),
                         expon.rvs(scale=daughter.lifetime))

    def _gen_momentum(self):

        P_m = PKLorentzVector(self._dec.Mother.mass, 0, 0, 0)

        self._tree._fill_branch('{}_PX'.format(self._dec.Mother.name), P_m.X[1].value)
        
        self._tree._fill_branch('{}_PY'.format(self._dec.Mother.name), P_m.X[2].value)
        self._tree._fill_branch('{}_PZ'.format(self._dec.Mother.name), P_m.X[3].value)
        self._tree._fill_branch('{}_P'.format(self._dec.Mother.name), pow(P_m.X[3].value**2+P_m.X[1].value**2+P_m.X[2].value**2, 0.5))
        _M_pt = pow(P_m.X[1].value**2+P_m.X[2].value**2, 0.5)
 
        _M_theta = atan(_M_pt/float(P_m.X[3].value)) if P_m.X[3].value != 0 else -9999

        _M_phi = atan(P_m.X[1].value/float(P_m.X[2].value)) if P_m.X[2].value != 0 else -9999

        _totE = self._dec.Mother.mass #Start with mass of the mother as total available energy

        tot_p_x_sq = uniform(0, _totE**2)
        _totE = pow(_totE**2-tot_p_x_sq**2, 0.5)
        tot_p_y_sq = uniform(0, _totE**2)
        tot_p_z_sq = pow(_totE**2-tot_p_y_sq**2, 0.5)

        p_x_sq = uniform(0, tot_p_x_sq)
        p_y_sq = uniform(0, tot_p_y_sq)
        p_z_sq = uniform(0, tot_p_z_sq)

        self._tree._fill_branch('{}_PT'.format(self._dec.Mother.name), _M_pt)
        self._tree._fill_branch('{}_THETA'.format(self._dec.Mother.name), _M_theta)
        self._tree._fill_branch('{}_PHI'.format(self._dec.Mother.name), _M_phi)

        _D0 = PKLorentzVector(pow(self._dec.Daughters[0].mass**2+p_x_sq+p_y_sq+p_z_sq, 0.5), sq_rt(p_x_sq), sq_rt(p_y_sq), sq_rt(p_z_sq)) #First Daughter can take any values from the range
        self._tree._fill_branch('{}_PX'.format(self._dec.Daughters[0].name), _D0.X[1].value)
        self._tree._fill_branch('{}_PY'.format(self._dec.Daughters[0].name), _D0.X[2].value)
        self._tree._fill_branch('{}_PZ'.format(self._dec.Daughters[0].name), _D0.X[3].value)
        self._tree._fill_branch('{}_P'.format(self._dec.Daughters[0].name), pow(_D0.X[3].value**2+_D0.X[1].value**2+_D0.X[2].value**2, 0.5))
        _D0_pt = pow(_D0.X[1].value**2+_D0.X[2].value**2, 0.5)
 
        _D0_theta = atan(_D0_pt/float(_D0.X[3].value)) if _D0.X[3].value != 0 else -9999

        _D0_phi = atan(_D0.X[1].value/float(_D0.X[2].value)) if _D0.X[2].value != 0 else -9999

        self._tree._fill_branch('{}_PT'.format(self._dec.Daughters[0].name), _D0_pt)
        self._tree._fill_branch('{}_THETA'.format(self._dec.Daughters[0].name), _D0_theta)
        self._tree._fill_branch('{}_PHI'.format(self._dec.Daughters[0].name), _D0_phi)


        P_m -= _D0 #Subtract first daughter 4-vector from Mother 4-vector
 
        counter = 1

        while counter < len(self._dec.Daughters)-1: #Iterate through all other daughters bar the last
            tot_p_x_sq -= p_x_sq
            tot_p_y_sq -= p_y_sq
            tot_p_z_sq -= p_z_sq

            p_x_sq = uniform(0, tot_p_x_sq)
            p_y_sq = uniform(0, tot_p_y_sq)
            p_z_sq = uniform(0, tot_p_z_sq) 

            _D = PKLorentzVector(pow(self._dec.Daughters[counter].mass**2+p_x_sq+p_y_sq+p_z_sq, 0.5),
                                    sq_rt(p_x_sq), sq_rt(p_y_sq), sq_rt(p_z_sq))
 
            self._tree._fill_branch('{}_PX'.format(self._dec.Daughters[counter].name), _D.X[1].value)
            self._tree._fill_branch('{}_PY'.format(self._dec.Daughters[counter].name), _D.X[2].value)
            self._tree._fill_branch('{}_PZ'.format(self._dec.Daughters[counter].name), _D.X[3].value)
            self._tree._fill_branch('{}_P'.format(self._dec.Daughters[counter].name), pow(_D.X[3].value**2+_D.X[1].value**2+_D.X[2].value**2, 0.5))
            _D_pt = pow(_D.X[1].value**2+_D.X[2].value**2, 0.5)
 
            _D_theta = atan(_D_pt/float(_D.X[3].value)) if _D.X[3].value != 0 else -9999

            _D_phi = atan(_D.X[1].value/float(_D.X[2].value)) if _D.X[2].value != 0 else -9999

            self._tree._fill_branch('{}_PT'.format(self._dec.Daughters[counter].name), _D_pt)
            self._tree._fill_branch('{}_THETA'.format(self._dec.Daughters[counter].name), _D_theta)
            self._tree._fill_branch('{}_PHI'.format(self._dec.Daughters[counter].name), _D_phi)
  
            P_m -= _D #Remove from remaining 4-momentum

            counter += 1

        self._tree._fill_branch('{}_PX'.format(self._dec.Daughters[-1].name), P_m.X[1].value)
        self._tree._fill_branch('{}_PY'.format(self._dec.Daughters[-1].name), P_m.X[2].value)
        self._tree._fill_branch('{}_PZ'.format(self._dec.Daughters[-1].name), P_m.X[3].value)
        self._tree._fill_branch('{}_P'.format(self._dec.Daughters[-1].name), pow(P_m.X[3].value**2+P_m.X[1].value**2+P_m.X[2].value**2, 0.5))
        _D_pt = pow(P_m.X[1].value**2+P_m.X[2].value**2, 0.5)
 
        _D_theta = atan(_D_pt/float(P_m.X[3].value)) if P_m.X[3].value != 0 else -9999

        _D_phi = atan(P_m.X[1].value/float(P_m.X[2].value)) if P_m.X[2].value != 0 else -9999

        self._tree._fill_branch('{}_PT'.format(self._dec.Daughters[-1].name), _D_pt)
        self._tree._fill_branch('{}_THETA'.format(self._dec.Daughters[-1].name), _D_theta)
        self._tree._fill_branch('{}_PHI'.format(self._dec.Daughters[-1].name), _D_phi)
  

    def __call__(self):
        self._logger.info("\tWill generate {} Events of type '{}'".format(self._nevts, self._dec.ID))
        for i in range(self._nevts):
            if i % 1000 == 0:
                self._logger.info("\tGenerating Event {}/{}".format(i, self._nevts)) 
            self._gen_time()
            self._gen_momentum()
        return self._tree
