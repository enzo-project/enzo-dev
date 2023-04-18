# This parses global_data.h using the library robot-cppheaderparser. It only
# needs to be run when a new parameter is created that also needs to be
# supplied to libyt, which will be very rarely.

import CppHeaderParser

# We 

def preprocess_header_file(fn = "global_data.h"):
    # This removes enzo-isms, specifically EXTERN.
    lines = []
    for line in open(fn):
        if line.startswith("EXTERN "): line = line[7:]
        lines.append(line)
    text = "\n".join(lines)
    return text

def handle_variable(var):
    if var['type'] == 'char *':
        if var['array']:
            print(var['name'], var['array_size'])
        else:
            print(var['name'])

if __name__ == "__main__":
    text = preprocess_header_file()
    gd = CppHeaderParser.CppHeader(text, argType="string")
