import os

filename = "Make.config.objects"
exts = ["C", "c", "F", "F90"]

# Read object names from Makefile
lines = open(filename).readlines()
files = []
for l in lines:
    if l.find(".o") > 0 and not(l.startswith("#")) and l.find("auto_show") < 0:
        source = l.strip().split(".")[0]
        files.append(source)

# Check if their sources exists
AllFound = True
for f in files:
    found = False
    for ext in exts:
        testfile = f + "." + ext
        if os.path.exists(testfile):
            found = True
            break
    if found == False:
        print "%s not found." % (f)
        AllFound = False

if AllFound:
    print "All source files in Makefile accounted for!"
