# This is the object that locates all *.enzo_test directories

import os.path
import optparse
import sys
import shutil

varspec = dict(
    name = (str, ''),
    answer_testing_script = (str, None),
    nprocs = (int, 1),
    runtime = (str, 'short'),
    critical = (bool, True),
    cadence = (str, 'nightly'),
    hydro = (bool, False),
    mhd = (bool, False),
    gravity = (bool, False),
    cosmology = (bool, False),
    chemistry = (bool, False),
    cooling = (bool, False),
    AMR = (bool, False),
    dimensionality = (int, 1),
    author = (str, ''),
    fullpath = (str, '.'),
    run_par_file = (str, None),
    fulldir = (str, '.'),
    max_time_minutes = (float, 1),
    radiation = (str, None),
)

known_variables = dict( [(k, v[0]) for k, v in varspec.items()] )
variable_defaults = dict( [(k, v[1]) for k, v in varspec.items()] )

run_template_dir = 'run_templates'
machines = {'local':       dict(script = 'local.run',
                                command = 'bash'),

            'nics-kraken': dict(script = 'nics-kraken.run',
                                command = 'qsub')}

template_vars = {'N_PROCS': 'nprocs',
                 'PAR_FILE': 'run_par_file',
                 'TEST_NAME': 'name'}

def add_files(my_list, dirname, fns):
    my_list += [os.path.join(dirname, fn) for
                fn in fns if fn.endswith(".enzotest")]

class EnzoTestCollection(object):
    def __init__(self, tests = None, machine = 'local'):
        if tests is None:
            # Now we look for all our *.enzotest files
            fns = []
            os.path.walk(".", add_files, fns)
            self.tests = []
            for fn in sorted(fns):
                print "HANDLING", fn
                self.add_test(fn)
        else:
            self.tests = tests
        self.test_container = []

    def go(self, output_dir, interleaved, machine):
        if interleaved:
            for my_test in self.tests:
                print "Preparing test: %s." % my_test['name']
                self.test_container.append(EnzoTestRun(output_dir, my_test, machine))
                self.test_container[-1].run_sim()
                self.test_container[-1].run_test()
        else:
            self.prepare_all_tests(output_dir, machine)
            self.run_all_sims()
            self.run_all_tests()

    def prepare_all_tests(self, output_dir, machine):
        print "Preparing all tests."
        for my_test in self.tests:
            print "Preparing test: %s." % my_test['name']
            self.test_container.append(EnzoTestRun(output_dir, my_test, machine))

    def run_all_sims(self):
        print "Running all simulations."
        for my_test in self.test_container:
            my_test.run_sim()

    def run_all_tests(self):
        print "Running all tests."
        for my_test in self.test_container:
            my_test.run_test()

    def add_test(self, fn):
        # We now do something dangerous: we exec the file directly and grab
        # its environment variables from it.
        local_vars = {}
        execfile(fn, {}, local_vars)
        test_spec = variable_defaults.copy()
        test_spec['fullpath'] = fn
        test_spec['fulldir'] = os.path.dirname(fn)
        test_spec['run_par_file'] = os.path.basename(test_spec['fulldir']) + ".enzo"
        for var, val in local_vars.items():
            if var in known_variables:
                caster = known_variables[var]
                if val == "False": val = False
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
        for t in self.tests:
            include = True
            for param, value in kwargs.items():
                if value == "None": value = None
                if value == "False": value = False
                if t.get(param, "Key Missing") != value:
                    include = False
                    break
            if include == True: pp.append(t)
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

class EnzoTestRun(object):
    def __init__(self, test_dir, test_data, machine, exe_path="../src/enzo/enzo.exe"):
        self.machine = machine
        self.test_dir = test_dir
        self.test_data = test_data
        self.exe_path = exe_path
        if self.exe_path is None:
            self.local_exe = None
        else:
            self.local_exe = os.path.basename(exe_path)

        self.run_dir = os.path.join(self.test_dir, self.test_data['fulldir'])

        self._copy_test_files()
        self._create_run_script()

    def _copy_test_files(self):
        # Check for existence
        if os.path.exists(self.run_dir):
            print "%s already exists. Skipping directory." % self.test_data['name']
        else:
            shutil.copytree(self.test_data['fulldir'], self.run_dir)
            if self.exe_path is not None:
                shutil.copy(self.exe_path, os.path.join(self.run_dir, self.local_exe))

    def _create_run_script(self):
        template_path = os.path.join(os.path.dirname(__file__), 
                                     run_template_dir, machines[self.machine]['script'])
        template_dest = os.path.join(self.run_dir, machines[self.machine]['script'])
        f = open(template_path, 'r')
        template = f.read()
        f.close()
        for var in template_vars.keys():
            template = template.replace(('${%s}' % var), 
                                        str(self.test_data[template_vars[var]]))
        template = template.replace('${EXECUTABLE}', "./%s" % self.local_exe)
        f = open(template_dest, 'w')
        f.write(template)
        f.close()

    def run_sim(self):
        print "Running test simulation: %s." % self.test_data['name']
        cur_dir = os.getcwd()
        # Check for existence
        if os.path.exists(self.run_dir+'/RunFinished'):
            print "%s run already completed, continuing..." % self.test_data['name']
            return
        
        os.chdir(self.run_dir) 
        command = "%s %s" % (machines[self.machine]['command'], 
                             machines[self.machine]['script'])
        print "Executing \"%s\"." % command
        os.system(command)
        os.chdir(cur_dir)

    def run_test(self):
        cur_dir = os.getcwd()
        os.chdir(self.run_dir)
        print "Running test: %s" % self.test_data['name']
        os.chdir(cur_dir)

class UnspecifiedParameter(object):
    pass
unknown = UnspecifiedParameter()

if __name__ == "__main__":
    etc = EnzoTestCollection()
    parser = optparse.OptionParser()
    parser.add_option("-o", "--output-dir", dest='output_dir',
                      help="Where to place the run directory")
    parser.add_option("--interleaved", action='store_true', dest='interleaved', default=False,
                      help="Option to interleave preparation, running, and testing.")
    parser.add_option("-m", "--machine", dest='machine', default='local', 
                      help="Machine to run tests on.")
    for var, caster in sorted(known_variables.items()):
        parser.add_option("", "--%s" % (var),
                          type=str, default = unknown)
    options, args = parser.parse_args()
    construct_selection = {}
    for var, caster in known_variables.items():
        if getattr(options, var) != unknown:
            val = getattr(options, var)
            if val == 'None': val = None
            if val == "False": val = False
            construct_selection[var] = caster(val)
    print
    print "Selecting with:"
    for k, v in sorted(construct_selection.items()):
        print "     %s = %s" % (k, v)
    etc2 = etc.select(**construct_selection)
    print
    print "\n".join(list(etc2.unique('name')))
    print "Total: %s" % len(etc2.tests)

    # Make it happen
    if not os.path.exists(options.output_dir): os.mkdir(options.output_dir)
    etc2.go(options.output_dir, options.interleaved, options.machine)
