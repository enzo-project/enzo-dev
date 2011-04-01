# This is the object that locates all *.enzo_test directories

import glob

known_variables = dict(
    name = str,
    answer_testing_script = str,
    nprocs = int,
    runtime = str,
    critical = bool,
    cadence = str,
    hydro = bool,
    mhd = bool,
    gravity = bool,
    cosmology = bool,
    AMR = bool,
    dimensionality = int
)

class EnzoTestCollection(object):
    def __init__(self):
        # Now we look for all our *.enzo_test files
        fns = glob.glob("**/*.enzo_test")
        self.tests = []
        for fn in sorted(fns):
            print "HANDLING", fn
            self.add_test(fn)

    def add_test(self, fn):
        # We now do something dangerous: we exec the file directly and grab
        # its environment variables from it.
        local_vars = {}
        execfile(fn, {}, local_vars)
        test_spec = {}
        for var, val in local_vars.items():
            if var in known_variables:
                caster = known_variables[var]
                test_spec[var] = caster(val)
            else:
                print "%s UNRECOGNIZED VARIABLE %s" % ( fn, var)
        self.tests.append(test_spec)

if __name__ == "__main__":
    etc = EnzoTestCollection()
