#!/usr/bin/env python
# This is the object that locates all *.enzo_test directories

import imp
import optparse
import os.path
import shutil
import signal
import subprocess
import sys
import time
import logging
sys.path.insert(0,'/Users/dccollins/local/src/yt-yt')

known_categories = [
    "Cooling",
    "Cosmology",
    "DrivenTurbulence3D",
    "FLD",
    "GravitySolver",
    "Hydro",
    "MHD",
    "RadiationTransport",
    "RadiationTransportFLD",
]

try:
    from yt.config import ytcfg
    from yt.utilities.answer_testing.api import \
        RegressionTestRunner, clear_registry, create_test, \
        TestFieldStatistics, TestAllProjections
    from yt.utilities.logger import ytLogger as mylog
    from yt.utilities.logger import \
        disable_stream_logging, ufstring
    disable_stream_logging()
    ytcfg["yt","suppressStreamLogging"] = "True"
except:
    raise

    RegressionTestRunner = None

try:
    import numpy
    numpy.seterr(all = "ignore")
except ImportError:
    pass

# Test keyword types and default values.
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
    max_time_minutes = (float, 60),
    radiation = (str, None),
)

known_variables = dict( [(k, v[0]) for k, v in varspec.items()] )
variable_defaults = dict( [(k, v[1]) for k, v in varspec.items()] )

# Machine specific job scripts.
run_template_dir = 'run_templates'
machines = {'local':       dict(script = 'local.run',
                                command = 'bash'),

            'nics_kraken': dict(script = 'nics_kraken.run',
                                command = 'qsub')}

# Map between job script variables and test keywords.
template_vars = {'N_PROCS'   : 'nprocs',
                 'PAR_FILE'  : 'run_par_file',
                 'TEST_NAME' : 'name',
                 'WALL_TIME' : 'run_walltime'}

results_filename = 'test_results.txt'
version_filename = 'version.txt'

# If we are able to, let's grab the ~/.enzo/machine_config.py file.
try:
    f, filename, desc = imp.find_module("machine_config",
                          [os.path.expanduser("~/.enzo/")])
    machine_config = imp.load_module("machine_config", f, filename, desc)
except ImportError:
    machine_config = None

def _get_hg_version(path):
    print "Getting current revision."
    from mercurial import hg, ui, commands 
    u = ui.ui() 
    u.pushbuffer() 
    repo = hg.repository(u, path) 
    commands.identify(u, repo) 
    return u.popbuffer()

def _to_walltime(ts):
    return "%02d:%02d:%02d" % \
        ((ts / 3600), ((ts % 3600) / 60), 
         (ts % 60))

def add_files(my_list, dirname, fns):
    my_list += [os.path.join(dirname, fn) for
                fn in fns if fn.endswith(".enzotest")]

class EnzoTestCollection(object):
    def __init__(self, tests = None, verbose=True):
        self.verbose = verbose
        if tests is None:
            # Now we look for all our *.enzotest files
            fns = []
            for cat in known_categories:
                os.path.walk(cat, add_files, fns)
            self.tests = []
            for fn in sorted(fns):
                if self.verbose: print "HANDLING", fn
                self.add_test(fn)
        else:
            self.tests = tests
        self.test_container = []

    def go(self, output_dir, interleaved, machine, exe_path, compare_dir,
           sim_only=False, test_only=False):
        go_start_time = time.time()
        self.output_dir = output_dir
        total_tests = len(self.tests)
        if interleaved:
            for i, my_test in enumerate(self.tests):
                print "Preparing test: %s." % my_test['name']
                self.test_container.append(EnzoTestRun(output_dir, my_test, 
                                                       machine, exe_path))
                if not test_only:
                    print "Running simulation: %d of %d." % (i, total_tests)
                    self.test_container[i].run_sim()
                if not sim_only:
                    print "Running test: %d of %d." % (i, total_tests)
                    self.test_container[-1].run_test(compare_dir)
        else:
            self.prepare_all_tests(output_dir, machine, exe_path)
            if not test_only: self.run_all_sims()
            if not sim_only: self.run_all_tests(compare_dir)
        if not sim_only: self.save_test_summary()
        go_stop_time = time.time()
        print "\n\nComplete!"
        print "Total time: %f seconds." % (go_stop_time - go_start_time)
        print "See %s/%s for a summary of all tests." % \
            (self.output_dir, results_filename)

    def prepare_all_tests(self, output_dir, machine, exe_path):
        print "Preparing all tests."
        for my_test in self.tests:
            print "Preparing test: %s." % my_test['name']
            self.test_container.append(EnzoTestRun(output_dir, my_test, 
                                                   machine, exe_path))

    def run_all_sims(self):
        total_tests = len(self.test_container)
        print "Running all simulations."
        for i, my_test in enumerate(self.test_container):
            print "Running simulation: %d of %d." % (i, total_tests)
            my_test.run_sim()

    def run_all_tests(self, compare_dir):
        total_tests = len(self.test_container)
        print "Running all tests."
        for i, my_test in enumerate(self.test_container):
            print "Running test: %d of %d." % (i, total_tests)
            my_test.run_test(compare_dir)

    def add_test(self, fn):
        # We now do something dangerous: we exec the file directly and grab
        # its environment variables from it.
        local_vars = {}
        execfile(fn, {}, local_vars)
        test_spec = variable_defaults.copy()
        test_spec['fullpath'] = fn
        test_spec['fulldir'] = os.path.dirname(fn)
        test_spec['run_par_file'] = os.path.basename(test_spec['fulldir']) + ".enzo"
        test_spec['run_walltime'] = _to_walltime(60 * test_spec['max_time_minutes'])
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

    def save_test_summary(self):
        all_passes = all_failures = 0
        run_passes = run_failures = 0
        dnfs = default_test = 0
        f = open(os.path.join(self.output_dir, results_filename), 'w')
        for my_test in self.test_container:
            default_only = False
            if my_test.run_finished:
                if my_test.test_data['answer_testing_script'] == 'None' or \
                        my_test.test_data['answer_testing_script'] is None:
                    default_only = True
                    default_test += 1
                t_passes = 0
                t_failures = 0
                for t_result in my_test.results.values():
                    t_passes += int(t_result)
                    t_failures += int(not t_result)
                f.write("%-70sPassed: %4d, Failed: %4d" % (my_test.test_data['fulldir'], 
                                                           t_passes, t_failures))
                if default_only:
                    f.write(" (default tests).\n")
                else:
                    f.write(".\n")
                all_passes += t_passes
                all_failures += t_failures
                run_passes += int(not (t_failures > 0))
                run_failures += int(t_failures > 0)
            else:
                dnfs += 1
                f.write("%-70sDID NOT FINISH\n" % my_test.test_data['fulldir'])

        f.write("\n")
        f.write("%-70sPassed: %4d, Failed: %4d.\n" % ("Total", 
                                                      all_passes, all_failures))
        f.write("Runs finished with all tests passed: %d.\n" % run_passes)
        f.write("Runs finished with at least one failure: %d.\n" % run_failures)
        f.write("Runs failed to complete: %d.\n" % dnfs)
        f.write("Runs finished with only default tests available: %d.\n" % default_test)
        f.close()

class EnzoTestRun(object):
    def __init__(self, test_dir, test_data, machine, exe_path):
        self.machine = machine
        self.test_dir = test_dir
        self.test_data = test_data
        self.exe_path = exe_path
        self.results = {}
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
            if options.clobber:
                print "%s exists, but clobber == True, so overwriting it." % self.test_data['name']
                shutil.rmtree(self.run_dir)
                shutil.copytree(self.test_data['fulldir'], self.run_dir)
                # Copy version file into run directory
                shutil.copy(os.path.join(self.test_dir, version_filename),
                            os.path.join(self.run_dir, version_filename))
                if self.exe_path is not None:
                    shutil.copy(self.exe_path, os.path.join(self.run_dir, self.local_exe))
            else:
                print "%s already exists. Skipping directory." % self.test_data['name']
        else:
            shutil.copytree(self.test_data['fulldir'], self.run_dir)
            # Copy version file into run directory
            shutil.copy(os.path.join(self.test_dir, version_filename),
                        os.path.join(self.run_dir, version_filename))
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
        for var, value in sorted(getattr(machine_config, self.machine, {}).items()):
            template = template.replace('${%s}' % var, value)
        f = open(template_dest, 'w')
        f.write(template)
        f.close()

    def run_sim(self):
        print "Running test simulation: %s." % self.test_data['fulldir']
        cur_dir = os.getcwd()
        # Check for existence
        if os.path.exists(os.path.join(self.run_dir, 'RunFinished')):
            print "%s run already completed, continuing..." % self.test_data['name']
            return
        
        os.chdir(self.run_dir)
        command = "%s %s" % (machines[self.machine]['command'], 
                             machines[self.machine]['script'])
        sim_start_time = time.time()
        # Run the command.
        proc = subprocess.Popen(command, shell=True, close_fds=True, 
                                preexec_fn=os.setsid)

        print "Simulation started on %s with maximum run time of %d seconds." % \
            (time.ctime(), (self.test_data['max_time_minutes'] * 60))
        running = 0
        # Kill the script if the max run time exceeded.
        while proc.poll() is None:
            if running > (self.test_data['max_time_minutes'] * 60):
                print "Simulation exceeded maximum run time."
                os.killpg(proc.pid, signal.SIGUSR1)
            running += 1
            time.sleep(1)
        
        sim_stop_time = time.time()
        if os.path.exists(os.path.join(self.run_dir, 'RunFinished')):
            f = open(os.path.join(self.run_dir, 'run_time'), 'w')
            f.write("%f seconds.\n" % (sim_stop_time - sim_start_time))
            f.close()
            print "Simulation completed in %f seconds." % \
                (sim_stop_time - sim_start_time)
        os.chdir(cur_dir)

    def run_test(self, compare_dir):
        cur_dir = os.getcwd()
        if compare_dir is None:
            compare_id = None
        else:
            compare_id = ""
            compare_dir = os.path.join(cur_dir, compare_dir,
                            self.test_data['fulldir'])
        os.chdir(self.run_dir)
        print "Running test: %s" % self.test_data['fulldir']
        self.run_finished = os.path.exists("RunFinished")

        if os.path.exists(results_filename):
            if self.run_finished:
                print "Reading test results from file."
                res_lines = file(results_filename)
                for line in res_lines:
                    if len(line.split()) == 2:
                        this_test, this_result = line.split()
                        self.results[this_test] = bool(this_result)

        else:
            fn = self.test_data['answer_testing_script']
            if RegressionTestRunner is None:
                print "This installation of yt does not support testing, please update to the branch 'yt'."
                return
            clear_registry()

            handler = logging.FileHandler("testing.log")
            f = logging.Formatter(ufstring)
            handler.setFormatter(f)
            mylog.addHandler(handler)
            if self.run_finished:
                if fn != 'None' and fn is not None:
                    if fn.endswith(".py"): fn = fn[:-3]
                    print "Loading module %s" % (fn)
                    f, filename, desc = imp.find_module(fn, ["."])
                    project = imp.load_module(fn, f, filename, desc)
                if fn is None or fn == "None":
                    create_test(TestFieldStatistics, "field_stats", tolerance = 1e-10)
                    create_test(TestAllProjections, "all_projs", tolerance = 1e-10)
                rtr = RegressionTestRunner("", compare_id,
                            compare_results_path = compare_dir)
                rtr.run_all_tests()
                self.results = rtr.passed_tests.copy()
            mylog.removeHandler(handler)
            handler.close()

        os.chdir(cur_dir)
        self.save_results()

    def save_results(self):
        f = open(os.path.join(self.run_dir, results_filename), 'w')
        if self.run_finished:
            if self.test_data['answer_testing_script'] == 'None' \
                    or self.test_data['answer_testing_script'] is None:
                f.write("Ran default binary tests since no others were available.\n")
            my_tests = self.results.keys()
            my_tests.sort()
            for my_test in my_tests:
                f.write("%-70s%s\n" % (my_test, self.results[my_test]))
        else:
            f.write("All tests failed because simulation did not finish.\n")
        f.close()

class UnspecifiedParameter(object):
    pass
unknown = UnspecifiedParameter()

if __name__ == "__main__":
    parser = optparse.OptionParser()
    parser.add_option("-c", "--compare-dir", dest='compare_dir',
                      default=None,
                      help="The directory structure to compare against")
    parser.add_option("--clobber", dest='clobber', default=False,
                      action="store_true", 
                      help="Recopies tests and tests from scratch.")
    parser.add_option("--interleave", action='store_true', dest='interleave', default=False,
                      help="Option to interleave preparation, running, and testing.")
    parser.add_option("-m", "--machine", dest='machine', default='local', 
                      help="Machine to run tests on.")
    parser.add_option("-o", "--output-dir", dest='output_dir',
                      help="Where to place the run directory")
    parser.add_option("--repo", dest='repository', default="../",
                      help="Path to repository being tested.")
    parser.add_option("--sim-only", dest='sim_only', action="store_true", 
                      default=False, help="Only run simulations.")
    parser.add_option("--test-only", dest='test_only', action="store_true", 
                      default=False, help="Only perform tests.")
    parser.add_option("-v", "--verbose", dest='verbose', action="store_true",
                      default=False, help="Slightly more verbose output.")
    for var, caster in sorted(known_variables.items()):
        parser.add_option("", "--%s" % (var),
                          type=str, default = unknown)
    options, args = parser.parse_args()

    etc = EnzoTestCollection(verbose=options.verbose)

    # Break out if output directory not specified.
    if options.output_dir is None:
        print 'Please enter an output directory with -o option'
        sys.exit(1)
    
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

    # get current revision
    options.repository = os.path.expanduser(options.repository)
    if options.compare_dir is not None:
        options.compare_dir = os.path.expanduser(options.compare_dir)
    hg_current = _get_hg_version(options.repository)
    rev_hash = hg_current.split()[0]
    options.output_dir = os.path.join(options.output_dir, rev_hash)
    if not os.path.exists(options.output_dir): os.makedirs(options.output_dir)
    f = open(os.path.join(options.output_dir, version_filename), 'w')
    f.write(hg_current)
    f.close()

    # the path to the executable we're testing
    exe_path = os.path.join(options.repository, "src/enzo/enzo.exe")

    # Make it happen
    etc2.go(options.output_dir, options.interleave, options.machine, exe_path,
            options.compare_dir, sim_only=options.sim_only, 
            test_only=options.test_only)
    try:
        import json
    except ImportError:
        json = None
    if json is not None and options.compare_dir is not None:
        f = open("results.js", "w")
        results = []
        for test in etc2.test_container:
            # This is to avoid any sorting code in JS
            vals = test.results.items()
            vals.sort()
            results.append( dict(name = test.test_data['name'],
                             results = vals) )
        f.write("test_data = %s;\n" % (json.dumps(results, indent=2)))
        f.write("compare_set = '%s';\ncurrent_set = '%s';\n" % (
                  options.compare_dir.strip(), hg_current.strip()))
