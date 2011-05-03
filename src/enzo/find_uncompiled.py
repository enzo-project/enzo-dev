import os, glob
import itertools

exts = ["C","F","c"]
objs = [line.strip().split(".")[0]
        for line in open("Make.config.objects")]

# Make.config.assemble also adds:

objs += ["enzo", "ProblemType_Python", "lcaperf", "python_bridge/problemtype_handler"]

uncompiled = []
for fn in itertools.chain(*(glob.glob("*.%s" % ext) for ext in exts)):
    if fn.split(".")[0] not in objs: uncompiled.append(fn)

for fn in sorted(uncompiled): print fn
