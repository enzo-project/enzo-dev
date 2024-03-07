# This parses global_data.h using the library robot-cppheaderparser. It only
# needs to be run when a new parameter is created that also needs to be
# supplied to libyt, which will be very rarely.

import CppHeaderParser

libyt_functions = {
    #"int" : "yt_set_UserParameterInt",
    "int" : "yt_set_UserParameterLongLong",
    "long" : "yt_set_UserParameterLong",
    "uint" : "yt_set_UserParameterUint",
    "Ulong" : "yt_set_UserParameterUlong",
    #"float" : "yt_set_UserParameterFloat",
    "float" : "yt_set_UserParameterDouble",
    "double" : "yt_set_UserParameterDouble",
    "FLOAT" : "yt_set_UserParameterDouble",
    "char" : "yt_set_UserParameterString", # This will need special case
    "hydro_method" : "yt_set_UserParameterLong"
}

skip_variables = [
    "LevelLookupTable",
    "ProblemTypeName",
    "GlobalCommunication",
    "RecvComm",
    "WaitComm",
    "timer",
    "counter",
    "starttime",
    "endtime",
    "Start_Wall_Time",
    "End_Wall_Time",
    "WallTime",
    "flagging_count",
    "in_count",
    "out_count",
    "moving_count",
    "flagging_pct",
    "moving_pct",
    "traceMEM",
    "StarParticlesOnProcOnLvl_Mass",
    "StarParticlesOnProcOnLvl_Type",
    "StageInput",
]

def preprocess_header_file(fn = "global_data.h"):
    # This removes enzo-isms, specifically EXTERN.
    lines = []
    for line in open(fn):
        if line.startswith("EXTERN "): line = line[7:]
        lines.append(line)
    text = "\n".join(lines)
    return text

def handle_loop(param_name, loop_length):
    prefix = []
    prefix.append(f"for (i = 0; i < {loop_length}; i++)" "{")
    prefix.append(f"    snprintf(tempname, 255, \"{param_name}_%d\", i);")
    suffix = ["}"]
    return prefix, suffix

def handle_char(var):
    # We need to special case this because we sometimes will have a string and
    # sometimes an array of strings.
    f = libyt_functions[var['raw_type']]
    lines = []#f"fprintf(stderr, \"Handling variable {var['name']}\\n\");"]
    if var['array']:
        if not var['pointer']: return []
        prefix, suffix = handle_loop(var['name'], var['array_size'])
        lines.extend(prefix)
        lines.append(f"    if (!{var['name']}[i]) continue;")
        lines.append(f"    {f}(tempname, {var['name']}[i]);")
        lines.extend(suffix)
    else:
        lines.append(f"if ({var['name']})")
        lines.append(f"    {f}(\"{var['name']}\", {var['name']});")
    return lines

def handle_variable(var):
    if var['name'] in skip_variables: return []
    if var.get('multi_dimensional_array'):
        # At present, multi_dimensional_arrays are not supported in this
        # library.  This touches a number of reasonably useful pieces of
        # information, but none are currently *used* within yt.
        return []
    if var['raw_type'] not in libyt_functions:
        # We don't know how to do typedefs quite yet, but we will soon.
        return []
    if var['raw_type'] == 'char':
        return handle_char(var)
    if not var['array']:
        reference = f"&{var['name']}"
        count = 1
    else:
        reference = var['name']
        count = var['array_size']
    f = libyt_functions[var['raw_type']]
    line = f"{f}(\"{var['name']}\", {count}, {reference});"
    return [line]


if __name__ == "__main__":
    text = preprocess_header_file('global_data.h')
    gd = CppHeaderParser.CppHeader(text, argType="string")
    lines = []
    for var in gd.variables:
        if var['name'] == 'AvoidRefineRegionLeftEdge':
            AVR = var
        lines.extend(handle_variable(var))
    with open("InitializeLibytInterface_finderfunctions.inc", "w") as f:
        f.write("\n".join(lines))
    #print("\n".join(lines))
