import os, sys

MAX_WIDTH = 35
TMPFILE = "all_tests"
OUTPUT = "TestProblems.rst"
gcode_base = "http://code.google.com/p/enzo/source/browse/run"

outf = open(OUTPUT, "w")

# Get all parameter files
os.chdir("../../../../run")
os.system("find . -name \"*.enzo\" > %s" % TMPFILE)
lines = open(TMPFILE).readlines()

data = {}
cwidth = {}
keys = ['dir', 'file', 'gcode']
for k in keys:
    data[k] = []
    cwidth[k] = 0

for l in lines:
    _d = os.path.dirname(l[:-1])
    data['dir'].append(os.path.relpath(_d))
    data['file'].append(os.path.basename(l[:-1]))
    gcode = gcode_base + '/' + os.path.relpath(l[:-1])
    data['gcode'].append(gcode)

for k in keys:    
    for d in data[k]:
        cwidth[k] = max(cwidth[k], len(d))
    cwidth[k] = min(cwidth[k], MAX_WIDTH)
cwidth['gcode'] = 7

# HEADER
for k in keys:
    outf.write(cwidth[k]*"=")
    outf.write(" ")
outf.write("\n")

headers = ["Directory", "Parameter File", "Source"]
for i,k in enumerate(keys):
    outf.write(headers[i])
    outf.write(" " * (cwidth[k] - len(headers[i]) + 1))
outf.write("\n")

for k in keys:
    outf.write(cwidth[k]*"=")
    outf.write(" ")
outf.write("\n")

for i in range(len(data['dir'])):

    if len(data['dir'][i]) > MAX_WIDTH:
        outf.write(data['dir'][i][:MAX_WIDTH-3])
        nchar = MAX_WIDTH
    else:
        outf.write(data['dir'][i])
        nchar = len(data['dir'][i])
    if len(data['dir'][i]) > MAX_WIDTH-3:
        if len(data['dir'][i]) <= MAX_WIDTH:
            _len = MAX_WIDTH - len(data['dir'][i])
            nchar += _len
            outf.write("." * _len)
        else:
            outf.write("...")
            nchar += 3
    outf.write(" " * (cwidth['dir'] - nchar))
    outf.write(" ")

    if len(data['file'][i]) > MAX_WIDTH:
        outf.write(data['file'][i][:MAX_WIDTH-3])
        nchar = MAX_WIDTH
    else:
        outf.write(data['file'][i])
        nchar = len(data['file'][i])
    if len(data['file'][i]) > MAX_WIDTH-3:
        if len(data['file'][i]) <= MAX_WIDTH:
            _len = MAX_WIDTH - len(data['file'][i])
            nchar += _len
            outf.write("." * _len)
        else:
            outf.write("...")
            nchar += 3
    outf.write(" " * (cwidth['file'] - nchar))
    outf.write(" ")

    link = "|link%d|_" % i
    outf.write(link)
    outf.write("\n")

# BOTTOM BAR
for k in keys:
    outf.write(cwidth[k]*"=")
    outf.write(" ")
outf.write("\n")

# all of the google code links

outf.write("\n")
for i,d in enumerate(data['gcode']):
    outf.write(".. |link%d| replace:: X\n" % i)
    outf.write(".. _link%d: %s\n" % (i,d))

outf.close()
os.remove(TMPFILE)
os.chdir("../doc/manual/source/user_guide/")
