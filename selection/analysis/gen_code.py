import re
import sys

# Parse variable declarations from the header file
def parse_header(filename):
    with open(filename) as f:
        lines = f.readlines()

    vars = []
    blacklist = {"fCurrent"}  # Exclude MakeClass internal variables

    for line in lines:
        line = line.strip()

        # Only match ROOT leaf types
        if not re.match(r'^(Int_t|Float_t|Double_t|Bool_t)\s+', line):
            continue
        # Skip comments and special markers
        if line.startswith("//") or "//! " in line:
            continue

        # Match variable definitions like: Float_t Bmass[1104];
        m = re.match(r'(\w+_t|\w+)\s+(\w+)(\[[^\]]+\])?;', line)
        if m:
            vartype = m.group(1)
            name = m.group(2)
            arr = m.group(3)

            if name in blacklist:
                continue  # skip unwanted variables

            vars.append((vartype, name, arr))

    return vars

# Map C++ types to ROOT branch type specifiers
def root_type(vartype):
    if vartype in ["Int_t", "int"]:
        return "/I"
    if vartype in ["Float_t", "float", "Double_t"]:
        return "/F"
    if vartype in ["Bool_t", "bool"]:
        return "/O"
    return "/D"   # default double

def main(filename):
    vars = parse_header(filename)

    # struct definition
    print("struct Out {")
    for t, n, arr in vars:
        print(f"  {t} {n};")
    print("};")
    print("\nOut out;")

    # Branch definitions
    print("\n// Branch definitions")
    for t, n, arr in vars:
        print(f'tree->Branch("{n}", &out.{n}, "{n}{root_type(t)}");')

    # Assignment code
    print("\n// Assignment code")
    for t, n, arr in vars:
        if arr:  # candidate-level array
            print(f"out.{n} = {n}[i_B];")
        else:    # event-level scalar
            print(f"out.{n} = {n};")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python gen_code.py analysis_DATA_pbpb.h")
    else:
        main(sys.argv[1])
