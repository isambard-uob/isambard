from functools import wraps
import warnings


def check_availability(program, test_func, global_settings):
    def function_grabber(f):
        @wraps(f)
        def function_with_check(*args, **kwargs):
            if program in global_settings:
                if 'available' not in global_settings[program]:
                    global_settings[program]['available'] = test_func()
                if global_settings[program]['available']:
                    return f(*args,  **kwargs)
            warning_string = ('{0} not found, side chains have not been packed.\n'
                              'Check that the path to the {0} binary in `.isambard_settings` is correct.\n'
                              'You might want to try rerunning `configure.py`').format(program)
            warnings.warn(warning_string, DependencyNotFoundWarning)
            return
        return function_with_check
    return function_grabber


class NoncanonicalWarning(RuntimeWarning):
    pass


class NotParameterisedWarning(RuntimeWarning):
    pass


class MalformedPDBWarning(RuntimeWarning):
    pass


class DependencyNotFoundWarning(RuntimeWarning):
    pass


warnings.simplefilter('always', DependencyNotFoundWarning)
warnings.simplefilter('once', PendingDeprecationWarning)
