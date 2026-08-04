import openmm
from openmm.app import Simulation

class DipoleReporter(object):
    def __init__(self, file, reportInterval):
        self._reportInterval = reportInterval
        self._out = open(file, "w")

    def describeNextReport(self, simulation):
        steps= self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, False, False, False, False, None)
    
    def report(self, simulation, state):
        system: openmm.System = simulation.system
        multipoles = [f for f in system.getForces() if isinstance(f, openmm.AmoebaMultipoleForce)][0]
        total = multipoles.getTotalDipoles(simulation.context)
        
        for f in total:
            self._out.write('%g %g %g\n' % (f[0], f[1], f[2]))
    
    def __del__(self):
        self._out.close()

class IndDipoleReporter(object):
    def __init__(self, file, reportInterval):
        self._reportInterval = reportInterval
        self._out = open(file, "w")

    def describeNextReport(self, simulation):
        steps= self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, False, False, False, False, None)
    
    def report(self, simulation, state):
        system: openmm.System = simulation.system
        multipoles = [f for f in system.getForces() if isinstance(f, openmm.AmoebaMultipoleForce)][0]
        induced = multipoles.getInducedDipoles(simulation.context)
        
        for f in induced:
            self._out.write('%g %g %g\n' % (f[0], f[1], f[2]))
    
    def __del__(self):
        self._out.close()

class PermDipoleReporter(object):
    def __init__(self, file, reportInterval):
        self._reportInterval = reportInterval
        self._out = open(file, "w")

    def describeNextReport(self, simulation):
        steps= self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, False, False, False, False, None)
    
    def report(self, simulation, state):
        system: openmm.System = simulation.system
        multipoles = [f for f in system.getForces() if isinstance(f, openmm.AmoebaMultipoleForce)][0]
        permanent = multipoles.getLabFramePermanentDipoles(simulation.context)
        
        for f in permanent:
            self._out.write('%g %g %g\n' % (f[0], f[1], f[2]))
    
    def __del__(self):
            self._out.close()