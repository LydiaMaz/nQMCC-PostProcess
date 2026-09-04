import re
import sys
import math
import matplotlib.pyplot as plt
#-----------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "figure.autolayout": True
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
    and finds all the times when the value expression is found, extracts it + error, and saves them.
    """
    
    r_vals0, r_vals1 = [], []
    nn0, nn0_err, np0, np0_err, pp0, pp0_err = [], [], [], [], [], []
    nn1, nn1_err, np1, np1_err, pp1, pp1_err = [], [], [], [], [], []

    current_spin = None
    in_block = False

    tot_avg = next((i for i in range(len(lines)-1, -1, -1) if "TOTAL AVERAGE" in lines[i]), None)
    if tot_avg is None:
        return None

    for i in range(tot_avg, len(lines)):
        line = lines[i]
        if "TWO-BODY R-SPACE DISTRIBUTIONS" in line:
            in_block = True
            continue
        
        if in_block:
            if "S=0" in line:
                current_spin = 0
                continue 
            elif "S=1" in line: 
                current_spin = 1 
                continue
            elif "PARTICLE COUNTS" in line: 
                break

        matches = VALUE_REGEX.findall(line)
        if len(matches) == 3 and current_spin is not None:
            parts = line.split()
            try:
                r = float(parts[0])
                if current_spin == 0:
                    r_vals0.append(r)
                    nn0.append(float(matches[0][0]))
                    nn0_err.append(float(matches[0][1]))
                    np0.append(float(matches[1][0]))
                    np0_err.append(float(matches[1][1]))
                    pp0.append(float(matches[2][0]))
                    pp0_err.append(float(matches[2][1]))
                elif current_spin == 1:
                    r_vals1.append(r)
                    nn1.append(float(matches[0][0]))
                    nn1_err.append(float(matches[0][1]))
                    np1.append(float(matches[1][0]))
                    np1_err.append(float(matches[1][1]))
                    pp1.append(float(matches[2][0]))
                    pp1_err.append(float(matches[2][1]))
            except ValueError:
                continue

    return {
        "r0": r_vals0, "r1": r_vals1,
        "nn0": nn0, "nn0_errs": nn0_err, "np0": np0, "np0_errs": np0_err, "pp0": pp0, "pp0_errs": pp0_err,
        "nn1": nn1, "nn1_errs": nn1_err, "np1": np1, "np1_errs": np1_err, "pp1": pp1, "pp1_errs": pp1_err,
    }
#-----------------------------------------------------------------------
def build_spin_figure(spin_label, isotopes, grouped_isotope, markers, colors, r_key, np_key, np_err_key, pp_key, pp_err_key):
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True, sharey='row')
    ax_np_iso1, ax_np_iso2 = axes[0]
    ax_pp_iso1, ax_pp_iso2 = axes[1]

    for ax in [ax_np_iso1, ax_np_iso2, ax_pp_iso1, ax_pp_iso2]:
        ax.tick_params(direction='in', top=True, right=True, length=10, labelsize=18)
        ax.set_xlim(0, 4)
    ax_pp_iso1.set_xlabel(r'$\mathrm{r\ (fm)}$', fontsize=20)
    ax_pp_iso2.set_xlabel(r'$\mathrm{r\ (fm)}$', fontsize=20)

    ax_np_iso1.set_ylabel(r'$\mathrm{\rho_{np}(r)}$', fontsize=20)
    ax_pp_iso1.set_ylabel(r'$\mathrm{\rho_{pp}(r)}$', fontsize=20)

    ax_np_iso1.set_title(isotopes[0], fontsize=26, pad=7, y=0.87, va="top")
    ax_np_iso2.set_title(isotopes[1], fontsize=26, pad=7, y=0.87, va="top")

    left_edge = ax_np_iso1.get_position().x0
    right_edge = ax_np_iso2.get_position().x1
    center_x = (left_edge + right_edge) / 2

    fig.suptitle(f"S = {spin_label}", fontsize=26, y=.98, x=0.53)

    isotope_axes = {
        isotopes[0]: (ax_np_iso1, ax_pp_iso1),
        isotopes[1]: (ax_np_iso2, ax_pp_iso2),
    }

    potentials_found = set()

    for isotope in isotopes:
        ax_np, ax_pp = isotope_axes[isotope]
        entries = grouped_isotope[isotope]

        for i, entry in enumerate(entries):
            potential = entry["potential"]
            data = entry["data"]
            potentials_found.add(potential.replace("+", "plus").replace("-", "_").replace("*", "star"))

            this_marker = markers[i % len(markers)]
            this_color = colors[i % len(colors)]

            ax_np.errorbar(data[r_key], data[np_key], yerr=data[np_err_key],
                    fmt=this_marker, color=this_color, ecolor=this_color, markersize=4,
                    elinewidth=1, capsize=2, label=potential)
            ax_pp.errorbar(data[r_key], data[pp_key], yerr=data[pp_err_key],
                    fmt=this_marker, color=this_color, ecolor=this_color, markersize=4,
                    elinewidth=1, capsize=2, label=potential)

    # y-axis ranges — keep consistent across isotopes, same fix as before
    #np_y_ticks = [0.1, 0.2, 0.3]
    #ax_np_iso1.set_yticks(np_y_ticks)
    #ax_np_iso1.set_ylim(0, 0.35)

    #pp_y_ticks = [0.02,0.04, 0.06]
    #ax_pp_iso1.set_yticks(pp_y_ticks)
    #ax_pp_iso1.set_ylim(0, 0.075)
    ax_np_iso1.locator_params(axis='y', nbins=6)
    ax_pp_iso1.locator_params(axis='y', nbins=6)

    # hide redundant right-column y labels
    ax_np_iso2.tick_params(labelleft=False)
    ax_pp_iso2.tick_params(labelleft=False)

    # hide seam-colliding x tick label on right column
    ax_np_iso2.xaxis.get_major_ticks()[0].label1.set_visible(False)
    ax_pp_iso2.xaxis.get_major_ticks()[0].label1.set_visible(False)
    ax_pp_iso1.xaxis.get_major_ticks()[4].label1.set_visible(False)

    ax_np_iso1.tick_params(right=False)
    ax_pp_iso1.tick_params(top=False, right=False)
    ax_pp_iso2.tick_params(top=False)

    ax_np_iso1.legend(fontsize=14, loc='upper right', frameon=False)

    plt.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)

    pot_str = "_".join(sorted(potentials_found))
    iso_str = "_".join(re.sub(r'[\$\{\}\\\(\)\^]|mathrm', '', iso) for iso in isotopes)
    filename = f"4panel_density_S{spin_label}_{iso_str}_{pot_str}.pdf"
    fig.savefig(filename, dpi=300)
    plt.close(fig)
    print(f"  Saved S={spin_label} 4-panel plot: {filename}")

#-----------------------------------------------------------------------
def main():
    if not FILE_PATHS:
        print("Error: Missing path(s). Please add file paths.")
        sys.exit(1)

    grouped_isotope = {} 
    markers = ['o', 's', '^', 'd', 'x', '*']
    colors = ['red', 'blue', 'green', 'purple', 'orange', 'pink']

    for path in FILE_PATHS: #file reading loop
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

        if isotope not in grouped_isotope:
            grouped_isotope[isotope] = []

        grouped_isotope[isotope].append({
            "potential": potential,
            "data": data
        })
    isotopes = list(grouped_isotope.keys())
    if len(isotopes) != 2:
        print(f"ERROR: Expected 2 isotopes, found {len(isotopes)}: {isotopes}")
        sys.exit(1)

    build_spin_figure("0", isotopes, grouped_isotope, markers, colors,
                       r_key="r0", np_key="np0", np_err_key="np0_errs",
                       pp_key="pp0", pp_err_key="pp0_errs")

    build_spin_figure("1", isotopes, grouped_isotope, markers, colors,
                       r_key="r1", np_key="np1", np_err_key="np1_errs",
                       pp_key="pp1", pp_err_key="pp1_errs")

#-----------------------------------------------------------------------
if __name__ == "__main__":
    main()
#-----------------------------------------------------------------------