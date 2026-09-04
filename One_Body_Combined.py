import re
import sys
import math
import matplotlib.pyplot as plt
#-----------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "figure.autolayout": False
})
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
# ISOTOPE_MAP can be extended if more isotopes are needed. 
ISOTOPE_MAP = {"h": "H", "he": "He", "li": "Li", "be": "Be", "b": "B", "c": "C", "n": "N"}
#-----------------------------------------------------------------------
# ISOTOPE_REGEX: matches .out line of form like "work/c13/106/c13_opt.calc"
ISOTOPE_REGEX = re.compile(r'work[/\\]([a-zA-Z]+)(\d+)[/\\]')
#VALUE_REGEX: matches .out line of the form "0.2977E-01 (0.1215E-01)"
VALUE_REGEX = re.compile(r'([0-9.E+-]+)\s*\(([0-9.E+-]+)\)')
#-----------------------------------------------------------------------
def get_isotope(lines):
    for line in lines:
        if m := ISOTOPE_REGEX.search(line):
            element, mass = m.group(1), m.group(2)
            symbol = ISOTOPE_MAP.get(element, element.capitalize())
            return f"$^{{{mass}}}\\mathrm{{{symbol}}}$"
    return "Unknown"
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
def get_density(lines):
    """
    Parses the output file from the bottom up to find the "TOTAL AVERAGE" line. 
    Then it parses through the block following the line to identify the " ONE-BODY R-SPACE DISTRIBUTIONS" section
    and finds all the times when the value expression is found, extracts, adds neutron and protons, writes the error, 
    and saves them.
    """

    r_vals, n_vals, n_errs, p_vals, p_errs = [], [], [], [], []
    in_block = False

    tot_avg = next((i for i in range(len(lines)-1, -1, -1) if "TOTAL AVERAGE" in lines[i]), None)
    if tot_avg is None:
        return None

    for i in range(tot_avg, len(lines)):
        line = lines[i]
        if "ONE-BODY R-SPACE DISTRIBUTIONS" in line:
            in_block = True
            continue
        if in_block:
            if "---" in line and r_vals:
                break
                
        matches = VALUE_REGEX.findall(line)
        if len(matches) == 4:
            parts = line.split()
            try:
                r_vals.append(float(parts[0]))

                n_down, n_down_err = float(matches[0][0]), float(matches[0][1])
                n_up, n_up_err     = float(matches[1][0]), float(matches[1][1])
                p_down, p_down_err = float(matches[2][0]), float(matches[2][1])
                p_up, p_up_err     = float(matches[3][0]), float(matches[3][1])

                n_vals.append(n_down + n_up)
                p_vals.append(p_down + p_up)

                n_errs.append(math.sqrt(n_down_err**2 + n_up_err**2))
                p_errs.append(math.sqrt(p_down_err**2 + p_up_err**2))
            except ValueError:
                continue

    return {"r": r_vals, "neutron": n_vals, "neutron_err": n_errs, "proton": p_vals, "proton_err": p_errs} if r_vals else None
#-----------------------------------------------------------------------
def main():
    if not FILE_PATHS:
        print("Error: Missing path(s). Please add file paths.")
        sys.exit(1)

    grouped_isotope = {} 
    markers = ['o', 's', '^', 'd', 'x', '*']
    colors = ['red', 'blue', 'green', 'purple', 'orange', 'pink']

    for path in FILE_PATHS:
        print(f"Reading file: {path}")
        try:
            with open(path, 'r') as f:
                file_lines = f.readlines()
        except IOError:
            print(f"ERROR: Skipping unreadable file path: {path}")
            continue

        isotope = get_isotope(file_lines)
        potential = get_potential(file_lines)
        data = get_density(file_lines)

        if data is None:
            continue

        grouped_isotope.setdefault(isotope, []).append({
            "potential": potential,
            "data": data
        })

    isotopes = list(grouped_isotope.keys())
    if len(isotopes) != 2:
        print(f"ERROR: Expected 2 isotopes for a 4-panel figure, found {len(isotopes)}: {isotopes}")
        sys.exit(1)

    fig, axes = plt.subplots(2, 2, figsize=(12, 11), sharey='row')
    #top row = neutron densities (both isotopes)
    #bottom row = proton densities (both isotopes)
    ax_n_iso1, ax_n_iso2 = axes[0]
    ax_p_iso1, ax_p_iso2 = axes[1] 
    

    for ax in [ax_n_iso1, ax_n_iso2, ax_p_iso1, ax_p_iso2]:
        ax.tick_params(direction='in', top=True, right=True,  length=10, labelsize=18)
        ax.set_xlabel(r'$\mathrm{r\ (fm)}$', fontsize=20)
        ax.set_xlim(0, 4)
        ax.locator_params(axis='y', nbins=10)

    x_ticks = [1, 2, 3, 4]
    ax_p_iso2.set_xticks(x_ticks)

    ax_n_iso1.tick_params(labelbottom=False)
    ax_n_iso2.tick_params(labelbottom=False)

    ax_n_iso1.set_ylabel(r'$\mathrm{\rho_n(r)}$', fontsize=20)
    ax_p_iso1.set_ylabel(r'$\mathrm{\rho_p(r)}$', fontsize=20)

    ax_n_iso1.tick_params(right=False)
    ax_p_iso1.tick_params(right=False, top=False)  
    ax_n_iso2.tick_params(labelleft=False)
    ax_p_iso2.tick_params(top=False, labelleft=False)

    #ax_n_iso1.set_title(isotopes[0], fontsize=20, pad=7)
    ax_n_iso1.text(0.45, 0.92, isotopes[0], transform=axes[0][0].transAxes, 
                     fontsize=20, ha='left', va='top')
    #ax_n_iso2.set_title(isotopes[1], fontsize=20, pad=7)
    ax_n_iso2.text(1.45, 0.92, isotopes[1], transform=axes[0][0].transAxes, 
                     fontsize=20, ha='left', va='top')

    all_potentials = set()

    isotope_axes = {
        isotopes[0]: (ax_n_iso1, ax_p_iso1),  # iso1 column
        isotopes[1]: (ax_n_iso2, ax_p_iso2),  # iso2 column
    }

    for isotope in isotopes:
        ax_n, ax_p = isotope_axes[isotope]
        entries = grouped_isotope[isotope]

        for i, entry in enumerate(entries):
            potential = entry["potential"]
            data = entry["data"]
            all_potentials.add(potential.replace("+", "plus").replace("-", "_").replace("*", "star"))

            #clean_pot = potential.replace("+", "plus").replace("-", "_").replace("*", "star")
            #potentials_found.add(clean_pot) 

            this_marker = markers[i % len(markers)]
            this_color = colors[i % len(colors)]

            ax_p.errorbar(data["r"], data["proton"], yerr=data["proton_err"], 
                    fmt=this_marker, color=this_color, ecolor=this_color, markersize=4, 
                    elinewidth=1, capsize=2, label=potential)

            ax_n.errorbar(data["r"], data["neutron"], yerr=data["neutron_err"], 
                    fmt=this_marker, color=this_color, ecolor=this_color, markersize=4, 
                    elinewidth=1, capsize=2, label=potential)

        clean_iso = re.sub(r'[\$\{\}\\\(\)\^]|mathrm', '', isotope)
       
    ax_n_iso1.legend(fontsize=16, frameon=False)

    pot_str = "_".join(sorted(all_potentials))
    iso_str = "_".join(re.sub(r'[\$\{\}\\\(\)\^]|mathrm', '', iso) for iso in isotopes)
    combined_filename = f"1body_density_{iso_str}_{pot_str}.pdf"
    plt.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0) 
    fig.savefig(combined_filename)
    plt.close(fig)
    print(f"  Saved 4-panel plot: {combined_filename}")
#-----------------------------------------------------------------------
if __name__ == "__main__":
    main()
#-----------------------------------------------------------------------