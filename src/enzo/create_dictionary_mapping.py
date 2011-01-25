import sys, glob

# Enzo's configuration simply does not scale.  There's no schema, there's no
# correlation between writing and reading, and there's no mechanism for
# programmatically accessing any value other than through global symbol lookup.
#
# This file will attempt to parse global_data and create the parameter lookup
# method used in InitializePythonInterface.C.  It doesn't do anything with
# pointers or with arrays.
#
# As a note, this routine will probably miss a few parameters, as it neither
# implements a proper parser nor makes any effort to look for line
# continuation.
#
#       -- Matt Turk

stack = []

def sanitize_parameter_name(vn):
    return vn.replace(";","").replace(",","")

def parse_file(fn):
    all_vars = []
    for line in open(fn):
        line = line[:line.rfind("//")]
        if line.startswith("#include") and "<" not in line:
            all_vars += parse_file(line.split('"')[1])
            continue
        elif line.startswith("#ifdef"):
            vn = line.split(None, 1)[1]
            all_vars.append(("IFDEF", vn))
            stack.append(vn)
            continue
        elif line.startswith("#endif"):
            if len(stack) > 0: all_vars.append(("ENDIF", stack.pop(-1)))
            continue
        elif not line.startswith("EXTERN"):
            continue
        extern, vtype, vars = line.split(None, 2)
        vars = vars.split()
        pointers = [v[0] == "*" for v in vars]
        arrays = ["[" in v for v in vars]
        if vtype.endswith("*"): pointers[0] = [True]
        for vn, vp, va in zip(vars, pointers, arrays):
            vn = vn.replace(";","")
            if vtype == "char" and va and vn.count("*") < 2:
                all_vars.append((vn[:vn.rfind("[")].replace("*",""), "string"))
            if va or vp: continue
            all_vars.append((sanitize_parameter_name(vn), vtype))
    return all_vars

# Note that we only use TEMP_PYFLOAT, even for "FLOAT" variables, because
# there is no long double variable type in Python, even though natively it can
# upcast to arbitrary precision as necessary.  This is a deficiency in the Enzo
# interface that I wrote, not in Python.
python_functions = {
    "FLOAT" : ("PyFloat_FromDouble", "double"),
    "float" : ("PyFloat_FromDouble", "double"),
    "int" : ("PyLong_FromLong", "long"),
    "string" : ("PyString_FromString", "char *"),
}

template = """
    if (strncmp(parameter_name, "%(parameter_name)s", %(name_length)i) == 0) {
        PyObject *rv = %(py_func)s((%(cast_type)s) %(parameter_name)s);
        return rv;
    }
"""

def finder_function(vname, vtype):
    if vname == "IFDEF": return "#ifdef %s\n" % vtype
    elif vname == "ENDIF": return "#endif /* %s */ \n" % (vtype)
    if vtype not in python_functions: return ""
    pf, cn = python_functions[vtype]
    kwargs = dict(parameter_name = vname,
                  name_length = len(vname) + 1,
                  py_func = pf,
                  cast_type = cn)
    return template % kwargs

if __name__ == "__main__":
    results = parse_file("global_data.h")
    output = open("InitializePythonInterface_finderfunctions.inc", "w")
    for vn, vt in results:
        a = finder_function(vn, vt)
        output.write(a)
    output.close()
