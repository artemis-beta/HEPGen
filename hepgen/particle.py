import pypdt

_all_particles = list(pypdt.PDT().ids())

def convert_name(name):
    _symbols = {')(' : '_', '(' : '_', ')' : '_',
                '*' : 'star', '+' : 'plus',
                '-' : 'minus'}
    for s in _symbols:
        name = name.replace(s, _symbols[s])
    return name
         
class particle(object):
     def __init__(self, _id):
          _part = pypdt.get(_id)
          self.symbol   = _part.name
          self.name     = convert_name(_part.name)
          self.charge   = _part.charge
          self.mass     = _part.mass
          self.lifetime = _part.get_lifetime()

class ParticleList(object):
    def __init__(self):
        for p in _all_particles:
            _particle = particle(p)
            setattr(self, _particle.name, particle(p))
            _anti = _particle.name.replace('plus', 'minus') if 'plus' in _particle.name else _particle.name.replace('minus', 'plus')
            _particle = particle(p)
            _particle.symbol = _particle.symbol.replace('+', '-')
            _particle.charge = _particle.charge*-1
            _particle.name = _anti
            setattr(self, _anti, _particle)

pdg = ParticleList()
