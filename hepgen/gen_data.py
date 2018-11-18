from hepgen.decays import decay_table 
from hepgen.data_tree import data_tree
from hepgen import __version__
from pktools.PKLorentzVector import PKLorentzVector
from scipy.stats import expon, norm
from random import uniform
from math import pi, atan, log, tan

import logging
logging.basicConfig()

def sq_rt(number):
    _num = uniform(-1,1)
    _vec = _num/abs(_num)
    return _vec*pow(number, 0.5)

class HEPGen(object):
    def __init__(self, decay_id, tree=None, nevts=1, energy=0):
        self._logger  = logging.getLogger(__class__.__name__)
        self._logger.setLevel('INFO')
        self._init()
        self._energy  = energy
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
        for var in ['TAU', 'PX', 'THETA', 'PHI', 'P', 'PE', 'PY', 'PZ', 'PT', 'ETA', 'M',
                    'FDX', 'FDY', 'FDZ', 'FD']:
            _tree.add_branch('{}_{}'.format(self._dec.Mother.name, var))
            for daughter in self._dec.Daughters:
                _tree.add_branch('{}_{}'.format(daughter.name, var))
        return _tree


    def _pseudorapidity(self, theta):
        try:
            if theta == -9999:
                raise ValueError
            return -log(tan(abs(theta/2.)))
        except ValueError:
            return -9999

    def _gen_data(self, boosted=0):
        p_x = uniform(0, boosted) if boosted != 0 else 0
        p_y = uniform(0, pow(boosted**2-p_x**2, 0.5)) if boosted != 0 else 0
        p_z = pow(boosted**2-p_x**2-p_y**2, 0.5) if boosted != 0 else 0

        P_m = PKLorentzVector(pow(self._dec.Mother.mass**2+p_x**2+p_y**2+p_z**2, 0.5), p_x, p_y, p_z)

        _meas_tau = expon.rvs(scale=self._dec.Mother.lifetime)

        _fac = 5.729E-29*1E-12/(1E6*1.911E-43*1E-3)

        _meas_px = P_m.X[1].value
        _meas_py = P_m.X[2].value
        _meas_pz = P_m.X[3].value

        _meas_pe = P_m.X[0].value

        _m = P_m.getMagnitude().value

        _gamma = _meas_pe/_m
        _dx  = _gamma*_meas_tau*(_meas_px/_m)*_fac
        _dy  = _gamma*_meas_tau*(_meas_py/_m)*_fac
        _dz  = _gamma*_meas_tau*(_meas_pz/_m)*_fac

        self._tree._fill_branch('{}_PE'.format(self._dec.Mother.name), _meas_pe)
        self._tree._fill_branch('{}_FDX'.format(self._dec.Mother.name), _dx)
        self._tree._fill_branch('{}_FDY'.format(self._dec.Mother.name), _dy)
        self._tree._fill_branch('{}_FDZ'.format(self._dec.Mother.name), _dz)
        self._tree._fill_branch('{}_FD'.format(self._dec.Mother.name), pow(_dx**2+_dy**2+_dz**2, 0.5))

        self._tree._fill_branch('{}_PX'.format(self._dec.Mother.name), _meas_px)
        self._tree._fill_branch('{}_M'.format(self._dec.Mother.name), _m)
        
        self._tree._fill_branch('{}_PY'.format(self._dec.Mother.name), _meas_py)
        self._tree._fill_branch('{}_PZ'.format(self._dec.Mother.name), _meas_pz)
        self._tree._fill_branch('{}_P'.format(self._dec.Mother.name), pow(_meas_px**2+_meas_py**2+_meas_pz**2, 0.5))
        _M_pt = pow(_meas_px**2+_meas_py**2, 0.5)
 
        _M_theta = atan(_M_pt/float(_meas_pz)) if _meas_pz != 0 else -9999

        _M_phi = atan(_meas_px/float(_meas_py)) if _meas_py != 0 else -9999

        _totE = self._dec.Mother.mass #Start with mass of the mother as total available energy

        tot_p_x_sq = uniform(0, _totE**2)
        _totE = pow(_totE**2-tot_p_x_sq, 0.5)
        tot_p_y_sq = uniform(0, _totE**2)
        tot_p_z_sq = pow(_totE**2-tot_p_y_sq, 0.5)

        p_x_sq = uniform(0, tot_p_x_sq)
        p_y_sq = uniform(0, tot_p_y_sq)
        p_z_sq = uniform(0, tot_p_z_sq)

        self._tree._fill_branch('{}_PT'.format(self._dec.Mother.name), _M_pt)
        self._tree._fill_branch('{}_THETA'.format(self._dec.Mother.name), _M_theta)
        self._tree._fill_branch('{}_ETA'.format(self._dec.Mother.name), self._pseudorapidity(_M_theta))
        self._tree._fill_branch('{}_PHI'.format(self._dec.Mother.name), _M_phi)

        _D0 = PKLorentzVector(pow(self._dec.Daughters[0].mass**2+p_x_sq+p_y_sq+p_z_sq, 0.5), sq_rt(p_x_sq), sq_rt(p_y_sq), sq_rt(p_z_sq)) #First Daughter can take any values from the range
        _meas_px = _D0.X[1].value
        _meas_py = _D0.X[2].value
        _meas_tau = expon.rvs(scale=self._dec.Daughters[1].lifetime)
        _meas_pz = _D0.X[3].value
        _meas_pe = _D0.X[0].value
        _m = _D0.getMagnitude().value

        _gamma = _meas_pe/_m
        _dx  = _gamma*_meas_tau*(_meas_px/_m)*_fac
        _dy  = _gamma*_meas_tau*(_meas_py/_m)*_fac
        _dz  = _gamma*_meas_tau*(_meas_pz/_m)*_fac
        self._tree._fill_branch('{}_FDX'.format(self._dec.Daughters[0].name), _dx)
        self._tree._fill_branch('{}_FDY'.format(self._dec.Daughters[0].name), _dy)
        self._tree._fill_branch('{}_FDZ'.format(self._dec.Daughters[0].name), _dz)
        self._tree._fill_branch('{}_FD'.format(self._dec.Daughters[0].name), pow(_dx**2+_dy**2+_dz**2, 0.5))
        self._tree._fill_branch('{}_PX'.format(self._dec.Daughters[0].name), _meas_px)
        self._tree._fill_branch('{}_M'.format(self._dec.Daughters[0].name), _m)
        self._tree._fill_branch('{}_PY'.format(self._dec.Daughters[0].name), _meas_py)
        self._tree._fill_branch('{}_PZ'.format(self._dec.Daughters[0].name), _meas_pz)
        self._tree._fill_branch('{}_PE'.format(self._dec.Daughters[0].name), _meas_pe)
        self._tree._fill_branch('{}_P'.format(self._dec.Daughters[0].name), pow(_meas_pz**2+_meas_px**2+_meas_py**2, 0.5))
        _D0_pt = pow(_meas_px**2+_meas_py**2, 0.5)
 
        _D0_theta = atan(_D0_pt/float(_meas_pz)) if _meas_pz != 0 else -9999

        _D0_phi = atan(_meas_px/float(_meas_py)) if _meas_py != 0 else -9999

        self._tree._fill_branch('{}_PT'.format(self._dec.Daughters[0].name), _D0_pt)
        self._tree._fill_branch('{}_THETA'.format(self._dec.Daughters[0].name), _D0_theta)
        self._tree._fill_branch('{}_ETA'.format(self._dec.Daughters[0].name), self._pseudorapidity(_D0_theta))
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
            _meas_px = _D.X[1].value
            _meas_py = _D.X[2].value
            _meas_pz = _D.X[3].value
            _meas_tau = expon.rvs(scale=self._dec.Daughters[counter].lifetime)
            _meas_pe = _D.X[0].value
            _m = _D.getMagnitude().value
            _gamma = _meas_pe/_m
            _dx  = _gamma*_meas_tau*(_meas_px/_m)*_fac
            _dy  = _gamma*_meas_tau*(_meas_py/_m)*_fac
            _dz  = _gamma*_meas_tau*(_meas_pz/_m)*_fac
            self._tree._fill_branch('{}_FDX'.format(self._dec.Daughters[counter].name), _dx)
            self._tree._fill_branch('{}_FDY'.format(self._dec.Daughters[counter].name), _dy)
            self._tree._fill_branch('{}_FDZ'.format(self._dec.Daughters[counter].name), _dz)
            self._tree._fill_branch('{}_FD'.format(self._dec.Daughters[counter].name), pow(_dx**2+_dy**2+_dz**2, 0.5))
 
            self._tree._fill_branch('{}_PX'.format(self._dec.Daughters[counter].name), _meas_px)
            self._tree._fill_branch('{}_PY'.format(self._dec.Daughters[counter].name), _meas_py)
            self._tree._fill_branch('{}_M'.format(self._dec.Daughters[counter].name), _m)
            self._tree._fill_branch('{}_PZ'.format(self._dec.Daughters[counter].name), _meas_pz)
            self._tree._fill_branch('{}_PE'.format(self._dec.Daughters[counter].name), _meas_pe)
            self._tree._fill_branch('{}_P'.format(self._dec.Daughters[counter].name), pow(_meas_pz**2+_meas_px**2+_meas_py**2, 0.5))
            _D_pt = pow(_meas_px**2+_meas_py**2, 0.5)
 
            _D_theta = atan(_D_pt/float(_meas_pz)) if _meas_pz != 0 else -9999

            _D_phi = atan(_meas_px/float(_meas_py)) if _meas_py != 0 else -9999

            self._tree._fill_branch('{}_PT'.format(self._dec.Daughters[counter].name), _D_pt)
            self._tree._fill_branch('{}_THETA'.format(self._dec.Daughters[counter].name), _D_theta)
            self._tree._fill_branch('{}_ETA'.format(self._dec.Daughters[counter].name), self._pseudorapidity(_D_theta))
            self._tree._fill_branch('{}_PHI'.format(self._dec.Daughters[counter].name), _D_phi)
  
            P_m -= _D #Remove from remaining 4-momentum

            counter += 1

        _meas_px = P_m.X[1].value
        _meas_py = P_m.X[2].value
        _meas_pz = P_m.X[3].value
        _meas_pe = P_m.X[0].value
        _meas_tau = expon.rvs(scale=self._dec.Daughters[-1].lifetime)
        _m = P_m.getMagnitude().value
        _gamma = _meas_pe/_m
        _dx  = _gamma*_meas_tau*(_meas_px/_m)*_fac
        _dy  = _gamma*_meas_tau*(_meas_py/_m)*_fac
        _dz  = _gamma*_meas_tau*(_meas_pz/_m)*_fac
        self._tree._fill_branch('{}_FDX'.format(self._dec.Daughters[-1].name), _dx)
        self._tree._fill_branch('{}_FDY'.format(self._dec.Daughters[-1].name), _dy)
        self._tree._fill_branch('{}_FDZ'.format(self._dec.Daughters[-1].name), _dz)
        self._tree._fill_branch('{}_FD'.format(self._dec.Daughters[-1].name), pow(_dx**2+_dy**2+_dz**2, 0.5))
        self._tree._fill_branch('{}_PX'.format(self._dec.Daughters[-1].name), _meas_px)
        self._tree._fill_branch('{}_PY'.format(self._dec.Daughters[-1].name), _meas_py)
        self._tree._fill_branch('{}_M'.format(self._dec.Daughters[-1].name), _m)
        self._tree._fill_branch('{}_PZ'.format(self._dec.Daughters[-1].name), _meas_pz)
        self._tree._fill_branch('{}_PE'.format(self._dec.Daughters[-1].name), _meas_pe)
        self._tree._fill_branch('{}_P'.format(self._dec.Daughters[-1].name), pow(_meas_pz**2+_meas_px**2+_meas_py**2, 0.5))
        _D_pt = pow(_meas_px**2+_meas_py**2, 0.5)
 
        _D_theta = atan(_D_pt/float(_meas_pz)) if _meas_pz != 0 else -9999

        _D_phi = atan(_meas_px/float(_meas_py)) if _meas_py != 0 else -9999

        self._tree._fill_branch('{}_PT'.format(self._dec.Daughters[-1].name), _D_pt)
        self._tree._fill_branch('{}_THETA'.format(self._dec.Daughters[-1].name), _D_theta)
        self._tree._fill_branch('{}_ETA'.format(self._dec.Daughters[-1].name), self._pseudorapidity(_D_theta))
        self._tree._fill_branch('{}_PHI'.format(self._dec.Daughters[-1].name), _D_phi)
  

    def __call__(self):
        self._logger.info("\tWill generate {} Events of type '{}'".format(self._nevts, self._dec.ID))
        for i in range(self._nevts):
            if i % 1000 == 0:
                self._logger.info("\tGenerating Event {}/{}".format(i, self._nevts)) 
            self._gen_data(self._energy)
        return self._tree
