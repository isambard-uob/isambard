import warnings


class NoncanonicalWarning(RuntimeWarning):
    pass


class NotParameterisedWarning(RuntimeWarning):
    pass


class MalformedPDBWarning(RuntimeWarning):
    pass


class DependencyNotFoundWarning(RuntimeWarning):
    pass

warnings.simplefilter('always', DependencyNotFoundWarning)
