#!/usr/bin/env python
import shutil, os

# Do not run the standard tests on these test problems.
ignore_list = ('GravityTest',)

template = open("test_type.py.template").read()

for root, dirs, files in os.walk("."):
    for fn in files:
        if fn.endswith(".enzotest") and \
          os.path.basename(fn)[:-9] not in ignore_list:
            simname = os.path.splitext(fn)[0]
            simpath = root
            testname = os.path.basename(fn)[:-9]
            oname = os.path.join(root, testname + "__test_standard.py")
            output = template % dict(filename = fn[:-4], simpath = simpath)
            open(oname, "w").write(output)
