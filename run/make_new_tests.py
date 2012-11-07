import shutil, os

template = open("test_type.py.template").read()

for root, dirs, files in os.walk("."):
    for fn in files:
        if fn.endswith(".enzo"):
            simname = os.path.splitext(fn)[0]
            simpath = root
            oname = os.path.join(root, "test_standard.py")
            output = template % dict(filename = fn, simpath = simpath)
            open(oname, "w").write(output)
