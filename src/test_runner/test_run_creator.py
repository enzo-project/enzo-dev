# This object creates a directory of tests and their appropriate submission

import glob
import os
import shutil

run_template_dir = 'run_templates'
machines = {'local': dict(script = 'local.sh',
                          command = 'bash'),

            'nics-kraken': dict(script = 'nics-kraken.sh',
                                command = 'qsub')}

template_vars = {'N_PROCS': 'nprocs',
                 'PAR_FILE': 'run_par_file'}
                 

class EnzoTestRun(object):
    def __init__(self, test_dir, test_data, machine='local', enzo_exe='./enzo.exe'):
        self.machine = machine
        self.test_dir = test_dir
        self.test_data = test_data
        self.enzo_exe = enzo_exe
        self.run_dir = os.path.join(self.test_dir, self.test_data['fulldir'])

        if not os.path.exists(self.test_dir): os.mkdir(self.test_dir)

        self._copy_test_files()
        self._create_run_script()

    def _copy_test_files(self):
        shutil.copytree(self.test_data['fulldir'], self.run_dir)

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
        template = template.replace('${ENZO_EXE}', self.enzo_exe)
        f = open(template_dest, 'w')
        f.write(template)
        f.close()

    def run_sim(self):
        os.chdir(self.run_dir)
        command = "%s %s" % (machines[self.machine]['command'], 
                             machines[self.machine]['script'])
        os.system(command)

if __name__ == "__main__":
    from test_finder import *
    etc = EnzoTestCollection()
    etc2 = etc.select(runtime="long")
    my_test = EnzoTestRun("../my_tests", etc2.tests[0])

