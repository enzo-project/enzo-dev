# This is the object that locates all *.enzo_test directories

import os.path

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
    dimensionality = int,
    fullpath = str,
    fulldir = str
)

def add_files(my_list, dirname, fns):
    my_list += [os.path.join(dirname, fn) for
                fn in fns if fn.endswith(".enzotest")]

class EnzoTestCollection(object):
    def __init__(self, tests = None):
        if tests is None:
            # Now we look for all our *.enzo_test files
            fns = []
            os.path.walk(".", add_files, fns)
            self.tests = []
            for fn in sorted(fns):
                print "HANDLING", fn
                self.add_test(fn)
        else:
            self.tests = tests

    def add_test(self, fn):
        # We now do something dangerous: we exec the file directly and grab
        # its environment variables from it.
        local_vars = {}
        execfile(fn, {}, local_vars)
        test_spec = dict(fullpath = fn,
                         fulldir = os.path.dirname(fn))
        for var, val in local_vars.items():
            if var in known_variables:
                caster = known_variables[var]
                test_spec[var] = caster(val)
                if val == "None": test_spec[var] = None
            else:
                print "%s UNRECOGNIZED VARIABLE %s" % ( fn, var)
        self.tests.append(test_spec)

    def unique(self, param):
        pp = set()
        for t in self.tests:
            pp.add(t.get(param, "Key Missing"))
        return pp

    def params(self):
        pp = set()
        for t in self.tests:
            pp.update(set(t.keys()))
        return pp

    def select(self, **kwargs):
        pp = []
        for param, value in kwargs.items():
            if value == "None": value = None
            for t in self.tests:
                if t.get(param, "Key Missing") == value:
                    pp.append(t)
        return EnzoTestCollection(tests = pp)

    def summary(self):
        for param in sorted(self.params()):
            if param.startswith("full"): continue
            print param
            for v in self.unique(param):
                print "     %s" % (v)
            print
        print
        print "NUMBER OF TESTS", len(self.tests)

if __name__ == "__main__":
    etc = EnzoTestCollection()
    etc2 = etc.select(runtime="long")
