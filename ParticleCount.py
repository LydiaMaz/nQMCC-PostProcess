import re
import sys
from decimal import Decimal, ROUND_HALF_UP
#-----------------------------------------------------------------------
# USER SHOULD CHANGE THESE PATHS. For this code, the paths MUST be in order. 
# If there are 2 isotopes, for example 13C and 13N, the order must be 
# 13C, 13C, 13C, ... , 13N, ... /NOT/ 13C, 13N, 13C, ... 
# If inputting multiple isotopes, go to def build_final_table(grouped_data) and uncomment the notes line. 
FILE_PATHS = [ 
    "/Users/camillejohnson/Desktop/nQMCC/Testing/outputs/c13_nv2+3_106.out",
    "/Users/camillejohnson/Desktop/nQMCC/Testing/outputs/c13_nv2+3_116.out",
    "/Users/camillejohnson/Desktop/nQMCC/Testing/outputs/c13_av18ux.out", 
    "/Users/camillejohnson/Desktop/nQMCC/Testing/outputs/n13_nv2+3_106.new.1.out", 
    "/Users/camillejohnson/Desktop/nQMCC/Testing/outputs/n13_nv2+3_116.out",
    "/Users/camillejohnson/Desktop/nQMCC/Testing/outputs/n13_av18ux.out", 
]
#-----------------------------------------------------------------------
PARTICLE_NAMES = ["PROTON-UP", "PROTON-DOWN", "NEUTRON-UP", "NEUTRON-DOWN"]
#-----------------------------------------------------------------------
LATEX_LABELS = {
    "PROTON-UP":    r"$N_{\uparrow p}$",
    "PROTON-DOWN":  r"$N_{\downarrow p}$",
    "NEUTRON-UP":   r"$N_{\uparrow n}$",
    "NEUTRON-DOWN": r"$N_{\downarrow n}$",
}
#-----------------------------------------------------------------------
# ISOTOPE_MAP can be extended if more isotopes are needed. 
ISOTOPE_MAP = {"h": "H", "he": "He", "li": "Li", "be": "Be", "b": "B", "c": "C", "n": "N"}
#-----------------------------------------------------------------------
# RegEx to match and find isotope and spin strings 
# ISOTOPE_REGEX: matches .out line of form like "work/c13/106/c13_opt.calc"
ISOTOPE_REGEX = re.compile(r'work[/\\]([a-zA-Z]+)(\d+)[/\\]')
# DATA_REGEX: matches .out line of form like 
DATA_REGEX = re.compile(r'([\d\.-]+)\s*\(([\d\.E\+-]+)\)')
#-----------------------------------------------------------------------
def get_isotope(filepath):
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if m := ISOTOPE_REGEX.search(line):
                    element, mass = m.group(1), m.group(2)
                    symbol = ISOTOPE_MAP.get(element, element.capitalize())
                    return f"$^{{{mass}}}\\mathrm{{{symbol}}}$"
    except IOError:
        pass
    return "Unknown"

#-----------------------------------------------------------------------
def rounding(val, err, base_decimals=3):
    """
    Rounds values and errors. The last input is the number of decimals desired. 
    Calling rounding(matches[0][0], matches[0][1], n) rounds to n decimals. 
    """
    n = base_decimals
    d_err = Decimal(str(err))
    
    if d_err != 0:
        while round(d_err * 10**n) == 0:
            n += 1
            
    val_rounded = Decimal(str(val)).quantize(Decimal(f'1.{"0"*n}'), rounding=ROUND_HALF_UP)
    err_int = int((d_err * 10**n).quantize(Decimal('1'), rounding=ROUND_HALF_UP))
    return f"${val_rounded}({err_int})$"

#-----------------------------------------------------------------------
def get_potential(filepath):
    """
    Finds the potential from the line (for example): 
    "pots/nv2_Ia.2bp            pots/nv3_Ia.3bp                  2B_POT_FILE, 3B_POT_FILE"
    and splits the line into four parts. Then, if an isotope name match is found, it replaces it. 
    If the 2 and 3 body potentials are the same, it joins them and tacks the suffix on to save spacing. 
    """
    two_bp, three_bp = "Unknown2BP", "Unknown3BP"

    two_bp_map = {
        "nv2_Ia_star": "NV2-Ia*", "nv2_Ia": "NV2-Ia", "nv2_Ib": "NV2-Ib", "nv2_IIa": "NV2-IIa", "nv2_IIb": "NV2-IIb",
        "nv8p_Ia": "NV8-Ia", "nv8p_Ib": "NV8-Ib", "nv8p_IIa": "NV8-IIa", "nv8p_IIb": "NV8-IIb",
        "av8prime": "AV8", "av14": "AV14", "av18": "AV18", "mt": "MT"
    }
    three_bp_map = {
        "nv3_Ia_star": "NV3-Ia*", "nv3_Ib_star": "NV3-Ib*", "nv3_IIa_star": "NV3-IIa*", "nv3_IIb_star": "NV3-IIb*",
        "nv3_Ia": "NV3-Ia", "nv3_Ib": "NV3-Ib", "nv3_IIa": "NV3-IIa", "nv3_IIb": "NV3-IIb",
        "none": "", "uix": "UIX", "ux": "UX"
    }

    try:
        with open(filepath, 'r') as f:
            for line in f:
                if "2B_POT_FILE" in line and "3B_POT_FILE" in line:
                    parts = line.split()
                    raw_two = parts[0] if parts else ""
                    raw_three = parts[1] if len(parts) > 1 else ""

                    for sub, rep in two_bp_map.items():
                        if sub in raw_two: two_bp = rep; break
                    for sub, rep in three_bp_map.items():
                        if sub in raw_three: three_bp = rep; break

                    if "NV2" in two_bp and "NV3" in three_bp:
                        two_suffix = two_bp.split('-')[-1]
                        three_suffix_star = three_bp.split('-')[-1].replace('*', '')
                        
                        if two_suffix == three_suffix_star:
                            star = "*" if "*" in three_bp else ""
                            return f"NV2+3-{two_suffix}{star}"
                    
                    return f"{two_bp}+{three_bp}" if three_bp and three_bp != "Unknown3BP" else two_bp
    except IOError:
        print(f"ERROR: Cannot read file to isolate potential: {filepath}")
    return f"{two_bp}{three_bp}"

#-----------------------------------------------------------------------
def get_particle_count(filepath):
    """
    Parses the output file from the bottom up to find the "TOTAL AVERAGE" line. 
    Then it parses through the block following the line to identify the "PARTICLE COUNTS" section, particularly:
    NEUTRON-DOWN            NEUTRON-UP             NEUTRONS
    3.6520 (0.9265E-02)     3.3480 (0.8368E-02)     7.0000 (0.1248E-01)... (or similar)
    then looks at the data line and extracts the values. 
    """
    results = {key: "NaN" for key in PARTICLE_NAMES}
    
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
    except IOError:
        print(f"ERROR: File not found: {filepath}")
        return results

    tot_avg = next((i for i in range(len(lines)-1, -1, -1) if "TOTAL AVERAGE" in lines[i]), None)
    if tot_avg is None:
        return results

    for i in range(tot_avg, len(lines)):
        line = lines[i]
        #matches[0] is neutron down, matches[1] is neutron up (or proton)
        for pair in [("NEUTRON-DOWN", "NEUTRON-UP"), ("PROTON-DOWN", "PROTON-UP")]:
            if pair[0] in line and pair[1] in line and (i + 1) < len(lines):
                matches = DATA_REGEX.findall(lines[i + 1])
                if len(matches) >= 2:
                    results[pair[0]] = rounding(matches[0][0], matches[0][1], 3)
                    results[pair[1]] = rounding(matches[1][0], matches[1][1], 3)
    return results

#-----------------------------------------------------------------------
def build_final_table(grouped_data):
    """
    Compiles all data and writes the LaTeX table
    """
    headers = ["Potential"] + [LATEX_LABELS[name] for name in PARTICLE_NAMES]
    
    out = [r"\begin{tabular}{lcccc}", r"\hline\hline", f"{' & '.join(headers)} \\\\"]
    
    for isotope, entries in grouped_data.items():
        #out.extend([r"\hline"]) 
        #For a multi-isotope header, swap with 
        out.extend([r"\hline", rf"\multicolumn{{5}}{{c}}{{\textbf{{{isotope}}}}} \\[2pt]", r"\hline"])
        
        for potential, data in entries:
            row = [potential] + [data.get(name, "NaN") for name in PARTICLE_NAMES]
            out.append(f"{' & '.join(row)} \\\\")
        
    out.extend([r"\hline\hline", r"\end{tabular}"])
    return "\n".join(out)

#-----------------------------------------------------------------------
def main():
    if not FILE_PATHS:
        sys.exit("Error: No valid output files specified.")

    grouped_data = {}
    for path in FILE_PATHS:
        isotope = get_isotope(path)
        pot = get_potential(path)
        metrics = get_particle_count(path)
        
        grouped_data.setdefault(isotope, []).append((pot, metrics))
       
    print("\n-------------------- LaTeX output --------------------\n")
    print("\n------------ Copy and paste the following ------------\n")
    print(build_final_table(grouped_data))
#-----------------------------------------------------------------------
if __name__ == "__main__":
    main()
#-----------------------------------------------------------------------