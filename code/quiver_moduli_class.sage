# Don't I first have to import quiver_class.sage? That doesn't work for me somehow. The only way I can make it work is load both quiver_class.sage and this file in sage and go from there. There must be a way to import this but I'm too stupid.

class Quiver_moduli:

    """
    A quiver moduli space is implemented by storing a quiver Q, a dimension vector d and a stability parameter theta. We want to distinguish between moduli spaces and moduli stacks (stack = False/True) and if we're considering the stable or the semi-stable version.
    We do not require theta*d = 0 but work with slope stability instead.

    Variables:
    quiver
    dimensionVector
    stabilityParameter
    stack = False
    version = 'sst'
    """

    def __init__(self, quiver, dimensionVector, stabilityParameter, stack=False, version='sst'):
        assert ( (quiver.number_of_vertices() == dimensionVector.length()) & (quiver.number_of_vertices() == stabilityParameter.length()) )
        self._quiver = quiver
        self._dimensionVector = dimensionVector
        self._stabilityParameter = stabilityParameter
        self._stack = stack
        assert (version in ['st','sst'])
        self._version = version

    def __repr__(self):
        output = ""
        if (self._version == 'sst'):
            output += "Semi-stable "
        else:
            output += "Stable "
        output += "quiver moduli "
        if self._stack:
            output += "stack "
        else:
            output += "space "
        output += "with:\n"+"Q = "+str(self._quiver)+"\n"+"d = "+str(self._dimensionVector)+"\n"+"theta = "+str(self._stabilityParameter)
        return output

    def quiver(self):
        return self._quiver

    def dimension_vector(self):
        return self._dimensionVector

    def stability_parameter(self):
        return self._stabilityParameter

    def slope(self, e):
        """The slope of e is defined as theta*e/(sum_i e_i). We need to ensure that e is non-negative and at least one entry is positive."""
        theta = self.stability_parameter()
        return (theta*e)/(sum(list(e)))
