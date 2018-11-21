import logging
logging.basicConfig()

from Herne.AppBase import HerneApp

import copy

class Generator(HerneApp):
    __version__ = 'v1r0p1'
    def __init__(self):
        self._members = copy.deepcopy(HerneApp._members)
        self._members += ['OutFile', 'DecID', 'nEvts', 'OutTree', 'Energy']
        HerneApp.__init__(self, 'Generator')

    def __call__(self):
        HerneApp.__call__(self)
        import yaml
        from hepgen.gen_data import HEPGen
        self._logger.info(" Will write to '{}'".format(self.OutFile))

        _gen_inst = HEPGen(self.DecID, tree=self.OutTree, nevts=self.nEvts, energy=self.Energy)

        yaml.dump(_gen_inst(), open(self.OutFile, 'w'))
