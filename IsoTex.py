import re
import sys
import pandas as pd
from decimal import Decimal, ROUND_HALF_UP

#-----------------------------------------------------------------------
# USER SHOULD CHANGE THESE PATHS 
FILE_PATHS = [ 
    "/Users/camillejohnson/Desktop/nQMCC/Testing/outputs/c13_nv2+3_106.out",
    "/Users/camillejohnson/Desktop/nQMCC/Testing/outputs/c13_nv2+3_116.out",
    "/Users/camillejohnson/Desktop/nQMCC/Testing/outputs/c13_av18ux.out", 
    "/Users/camillejohnson/Desktop/nQMCC/Testing/outputs/n13_nv2+3_106.new.1.out", 
    "/Users/camillejohnson/Desktop/nQMCC/Testing/outputs/n13_nv2+3_116.out",
    "/Users/camillejohnson/Desktop/nQMCC/Testing/outputs/n13_av18ux.out", 
]
#-----------------------------------------------------------------------
ENERGY_NAMES =  ["H", "Ti", "Vij", "Vijk", "Vem"]
RADIUS_NAMES =  ["POINT NEUTRON RMS RADIUS", "POINT PROTON  RMS RADIUS"]
MOMENT_NAMES =  ["NEUTRON MAGNETIC MOMENT CONTRIBUTION", "PROTON  MAGNETIC MOMENT CONTRIBUTION", "ORBITAL MAGNETIC MOMENT CONTRIBUTION", "MAGNETIC MOMENT"]
#-----------------------------------------------------------------------
RADIUS_KEY_MAP = {
    "POINT NEUTRON RMS RADIUS": "rn",
    "POINT PROTON  RMS RADIUS": "rp"
}
MOMENT_KEY_MAP = {
    "PROTON  MAGNETIC MOMENT CONTRIBUTION": "mup",
    "NEUTRON MAGNETIC MOMENT CONTRIBUTION": "mun",
    "ORBITAL MAGNETIC MOMENT CONTRIBUTION": "mul",
    "MAGNETIC MOMENT": "mu"
}
LATEX_LABELS = {
    "H":    r"$\langle H \rangle$", "Ti":   r"$\langle T \rangle$", "Vij":  r"$\langle V_{ij}\rangle $",
    "Vijk": r"$\langle V_{ijk} \rangle$", "Vem":  r"$\langle V_{\mathrm{em}} \rangle$",
    "rp": r"$r_p$", "rn": r"$r_n$", "mup": r"$\mu_p$", 
    "mun": r"$\mu_n$", "mul": r"$\mu_L$", "mu": r"$\mu$"
}
#-----------------------------------------------------------------------
# ISOTOPE_MAP can be extended if more isotopes are needed. 
ISOTOPE_MAP = {"h": "H", "he": "He", "li": "Li", "be": "Be", "b": "B", "c": "C", "n": "N"}
#-----------------------------------------------------------------------
# ISOTOPE_REGEX: matches .out line of form like "work/c13/106/c13_opt.calc"
ISOTOPE_REGEX = re.compile(r'work[/\\]([a-zA-Z]+)(\d+)[/\\]')
ENERGY_REGEX = re.compile(r'([-+]?\d+\.\d+(?:E[-+]?\d+)?)\s*\(([\d.E+-]+)\)')
#-----------------------------------------------------------------------
def get_isotope(lines):
    for line in lines:
        if m := ISOTOPE_REGEX.search(line):
            element, mass = m.group(1), m.group(2)
            symbol = ISOTOPE_MAP.get(element, element.capitalize())
            return f"$^{{{mass}}}\\mathrm{{{symbol}}}$"
    return "Unknown"
#-----------------------------------------------------------------------
def rounding(val, err, base_decimals=3):
    """
    Rounds values and errors, based on first non-zero digit of the error.
    It adds decimal places to the error until (error*10^n) isn't 0. 
    And assigns the decimals of the value to be n.
    The last input is automatically n=3, but can be adjusted. 
    
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
def get_potential(lines):
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

    for line in lines:
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
    return f"{two_bp}{three_bp}"
#-----------------------------------------------------------------------
def get_energy(lines):
    """
    Parses the output file from the bottom up to find the "TOTAL AVERAGE" line. 
    Then it parses through the block following the line to identify the " ONE-BODY R-SPACE DISTRIBUTIONS" section
    and finds all the times when the value expression is found, extracts it + error, and saves them.
    """
    results = {}
    remaining = set(ENERGY_NAMES)
    
    tot_avg = next((i for i in range(len(lines)-1, -1, -1) if "TOTAL AVERAGE" in lines[i]), None)
    if tot_avg is None:
        return None
    target_lines = lines[tot_avg:] if tot_avg is not None else lines 

    for line in target_lines:
        for key in list(remaining):
            if re.match(r'^\s*' + re.escape(key) + r'\s*=', line):
                e = ENERGY_REGEX.search(line)
                #if we find the correct pattern of characters, then we add that key to our results
                #group 1 is the main value, and group 2 is the MC error (placed in parenthesis)
                if e: 
                    try: 
                            value_1 = float(e.group(1))
                            value_2 = float(e.group(2))
                            base = 3 if key == "Vem" else 2
                            results[key] = rounding(value_1, value_2, base)
                            remaining.remove(key)
                    except (ValueError, IndexError) as err:
                        #Error message
                        print(f"Warning: not supported data input for key {key}: {err}")
        if not remaining:
            break
    for key in ENERGY_NAMES:
        if key not in results:
            results[key] = "NaN"
    return results

#-----------------------------------------------------------------------
def get_radius(lines):
    radius = {}
    remaining = set(RADIUS_NAMES)

    tot_avg = next((i for i in range(len(lines)-1, -1, -1) if "TOTAL AVERAGE" in lines[i]), None)
    if tot_avg is None:
        return None
    target_lines = lines[tot_avg:] if tot_avg is not None else lines 

    for line in target_lines:
        for key in list(remaining):
            if re.match(r'^\s*' + re.escape(key) + r'\s*=', line):
                m = ENERGY_REGEX.search(line)
                if m:
                    short_key = RADIUS_KEY_MAP[key]
                    val_float = float(m.group(1))
                    err_float = float(m.group(2))
                    radius[short_key] = rounding(val_float, err_float, 3)
                    remaining.remove(key)
    for key in RADIUS_NAMES:
        short_key = RADIUS_KEY_MAP[key]
        if short_key not in radius:
            radius[short_key] = "NaN"
    return radius

#-----------------------------------------------------------------------
def get_moments(lines):
    moment = {}
    remaining = set(MOMENT_NAMES)

    tot_avg = next((i for i in range(len(lines)-1, -1, -1) if "TOTAL AVERAGE" in lines[i]), None)
    if tot_avg is None:
        return None
    target_lines = lines[tot_avg:] if tot_avg is not None else lines 

    for line in lines:
        for key in list(remaining):
            if re.match(r'^\s*' + re.escape(key) + r'\s*=', line):
                m = ENERGY_REGEX.search(line)
                if m:
                    short_key = MOMENT_KEY_MAP[key]
                    val_float = float(m.group(1))
                    err_float = float(m.group(2))
                    moment[short_key] = rounding(val_float, err_float, 3)
                    remaining.remove(key)
    for key in MOMENT_NAMES:
        short_key = MOMENT_KEY_MAP[key]
        if short_key not in moment:
            moment[short_key] = "NaN"
    return moment

#-----------------------------------------------------------------------
def build_final_table(grouped_data):
    global_potentials = []
    for files in grouped_data.values(): #groups the potentials in order, ignoring isotopes.
        for potential, _ in files:
            if potential not in global_potentials:
                global_potentials.append(potential)
                
    num_cols = len(global_potentials) + 1
    col_formatting = "l" + "c" * (num_cols - 1)
    
    master_table = []
    master_table.append(r"\begin{tabular}{" + col_formatting+ r"}")
    master_table.append(r"\hline\hline")
    master_table.append(" & "+" & ".join(global_potentials)+r" \\")
    
    for isotope, files in grouped_data.items():
        master_table.append(r"\hline")
        master_table.append(rf"\multicolumn{{{num_cols}}}{{c}}{{\textbf{{{isotope}}}}} \\[2pt]")
        master_table.append(r"\hline")
        
        all_row_keys = ENERGY_NAMES + ["rp", "rn"] + ["mup", "mun", "mul", "mu"]
        table = {"": [LATEX_LABELS[k] for k in all_row_keys]}
        
        potential_map = {pot: data for pot, data in files}
        for pot in global_potentials:
            if pot in potential_map:
                table[pot] = [potential_map[pot][k] for k in all_row_keys]
            else:
                table[pot] = [""] * len(all_row_keys) 

        df = pd.DataFrame(table)
        latex_output = df.to_latex(columns=None, index=False, escape=None, column_format=col_formatting, caption=None, label=None)
        
        lines = latex_output.strip().split("\n")
        data_rows = [line for line in lines if line.strip().startswith("$")]
        
        for line in data_rows:
            if line.startswith(r"$r_p$ &"):
                master_table.append(r"\hline " + line)
            elif line.startswith(r"$\mu_p$ &"):
                master_table.append(r"\hline " + line)
            else:
                master_table.append(line)
        
    master_table.append(r"\hline\hline")
    master_table.append(r"\end{tabular}")

    return "\n".join(master_table)
#-----------------------------------------------------------------------
def main():
    if not FILE_PATHS:
        print("Error: Missing path(s). Please add file paths.")
        sys.exit(1)

    grouped_data = {}

    print("Processing files...")
    for path in FILE_PATHS:
        try:
            with open(path, 'r') as f:
                file_lines = f.readlines()
        except IOError:
            print(f"ERROR: Skipping unreadable file path: {path}")
            continue

        isotope = get_isotope(file_lines) 
        potential = get_potential(file_lines)
        energies = get_energy(file_lines)
        radii = get_radius(file_lines)
        moments = get_moments(file_lines)
        combined_data = {**energies, **radii, **moments}
        
        if isotope not in grouped_data:
            grouped_data[isotope] = []
            
        grouped_data[isotope].append((potential, combined_data))
        print(f"    Found {isotope} with potential: {potential}")
            
    print("\n-------------------- LaTeX output --------------------\n")
    print("\n------------ Copy and paste the following ------------\n")
    
    print(build_final_table(grouped_data))
#-----------------------------------------------------------------------
if __name__ == "__main__":
    main()
#-----------------------------------------------------------------------