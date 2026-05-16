import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import click
from PIL import Image
import io

# ==========================================
# === 🎨 CUSTOMIZATION HUB 🎨 ===
# Change these values to adjust the look of the plot!
# ==========================================

# 1. Image Size (Inches)
FIG_WIDTH_IN = 6.9  
FIG_HEIGHT_IN = 6.0 

# 2. Panel Size Ratios (Top to Bottom -> Epoch, Temp, Sea Level, PSMC)
# Adjusted to give Epoch, Temp, and Sea Level more vertical space, reducing PSMC
PANEL_RATIOS = [0.6, 1.2, 1.2, 2.5] 

# 3. Font Style and General Size
# Increased overall font size for better readability
FONT_SIZE = 10.5
FONT_FAMILY = 'Arial' 

# 4. Legend Position (X and Y coordinates relative to the whole figure)
LEGEND_POS_X = 0.5
LEGEND_POS_Y = 0.0008

# 5. Axis Title Positioning (Padding/distance from the axis tick marks)
LABEL_PAD_X_PSMC = 6
LABEL_PAD_Y_TEMP = 13
LABEL_PAD_Y_SEA = 0 
LABEL_PAD_Y_PSMC = 22

# 6. Y-Axis Limits for Top Panels
SEA_LEVEL_MIN = -140
SEA_LEVEL_MAX = 20

TEMP_MIN = -8 
TEMP_MAX = 4  

# 🚨 TEMPERATURE SHADING THRESHOLD (GLACIAL PERIODS) 🚨
# A drop of -2.0°C to -3.0°C is typically indicative of entering glacial conditions.
# The plot will draw vertical light blue bands through ALL panels where temperature is below this.
TEMP_COLD_THRESHOLD = -2.5 
COLD_SHADE_COLOR = '#EDF5FF' # Very light, seamless blue

# 7. Line Colors and Thicknesses (Main plots)
TEMP_LINE_COLOR = '#333333'
TEMP_LINE_WIDTH = 0.7

SEA_LEVEL_LINE_COLOR = '#333333'
SEA_LEVEL_LINE_WIDTH = 0.7

PSMC_LINE_WIDTH = 1.2 
BOOTSTRAP_LINE_WIDTH = 0.7
BOOTSTRAP_ALPHA = 0.08 

# 8. Vertical Event Markers (LGM)
MARKER_LINE_COLOR = 'black'
MARKER_LINE_WIDTH = 1.0
MARKER_LINE_STYLE = '--' 

MARKER_FONT_SIZE = 8.5 # Increased marker font size
MARKER_FONT_COLOR = 'black'
MARKER_FONT_WEIGHT = 'normal' 
MARKER_FONT_STYLE = 'normal'  

# 9. Epoch Bar Settings
EPOCH_FONT_SIZE = 7.5 # Increased epoch font size
EPOCH_LINE_COLOR = 'black'
EPOCH_LINE_WIDTH_THICK = 0.8  
EPOCH_LINE_WIDTH_THIN = 0.5   

# -> NUDGE TEXT POSITION: (Because X-axis is inverted, <1 moves text RIGHT, >1 moves text LEFT)
SHIFT_PLEISTOCENE = 1.0
SHIFT_HOLOCENE = 1.0

# Nudged the 'E' to the right (0.94) so it stays away from the left spine
SHIFT_EARLY_PLEISTOCENE = 0.98   
SHIFT_MIDDLE_PLEISTOCENE = 1.0
SHIFT_LATE_PLEISTOCENE = 1.0

SHIFT_EARLY_HOLOCENE = 1.0
SHIFT_MIDDLE_HOLOCENE = 1.0

# Export Resolution
DPI_EXPORT = 600    
# ==========================================

plt.rcParams.update({'font.size': FONT_SIZE, 'font.family': 'sans-serif'})
plt.rcParams['font.sans-serif'] = [FONT_FAMILY]

def load_timeseries_data(filepath):
    times = []
    values = []
    time_multiplier = 1000  

    with open(filepath, 'r') as f:
        lines = f.readlines()

    if not lines:
        return np.array([]), np.array([])

    header = lines[0].lower()
    if 'ma' in header:
        time_multiplier = 1_000_000
    elif 'kyr' in header or 'ka' in header:
        time_multiplier = 1_000
    elif 'yr' in header or 'year' in header:
        time_multiplier = 1

    for line in lines[1:]: 
        line = line.strip()
        if not line: 
            continue
        
        parts = line.split('\t')
        if len(parts) < 2:
            parts = line.split()
            
        if len(parts) >= 2:
            try:
                time_val = float(parts[0].replace(',', '.'))
                val = float(parts[1].replace(',', '.'))
                times.append(time_val * time_multiplier)
                values.append(val)
            except ValueError:
                continue 

    return np.array(times), np.array(values)


@click.command()
@click.option("--psmc_file_list", prompt=False, help="List of sample names and line colors (tab-separated)")
@click.option("--temperature_file", prompt=False, help="Temperature data file (tab-separated)")
@click.option("--sea_level_file", prompt=False, help="Sea level data file (tab-separated)")
@click.option("--bin_size", default=100)
@click.option("--mutation_rate", default=4.78e-9)
@click.option("--generation_time", default=20)
@click.option("--x_min", default=1e3, type=float)
@click.option("--x_max", default=1e6, type=float)
@click.option("--y_min", default=0, type=float)
@click.option("--y_max", default=60e3, type=float)
@click.option("--line_width", default=1.3, type=float) 
@click.option("--length", default=14, type=float)  
@click.option("--height", default=5, type=float)   
@click.option("--picture_dpi", default=300, type=int)  
@click.option("--span_min", default=0, type=float)
@click.option("--span_max", default=0, type=float)
@click.option("--span_color", default='purple')
@click.option("--span_alpha", default=0.3, type=float)
def plot_psmc_with_temperature_and_sea_level(psmc_file_list, temperature_file, sea_level_file, bin_size, mutation_rate, generation_time,
                                             x_min, x_max, y_min, y_max, line_width, length, height, picture_dpi, 
                                             span_min, span_max, span_color, span_alpha):
    # --- Load Data ---
    times_temp, temperatures = load_timeseries_data(temperature_file)
    times_sea_level, sea_levels = load_timeseries_data(sea_level_file)

    # --- Figure setup ---
    fig, (ax_epoch, ax_temp, ax_sea, ax_psmc) = plt.subplots(nrows=4, ncols=1, sharex=True, 
                                                 figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN),
                                                 gridspec_kw={'height_ratios': PANEL_RATIOS})
    
    fig.subplots_adjust(left=0.14, bottom=0.15, right=0.94, top=0.98, hspace=0.0)
    plt.style.use("classic")

    # ==========================================
    # 1. TOP PANE: Epochs 
    # ==========================================
    ax_epoch.set_ylim(0, 1)
    ax_epoch.set_yticks([]) 
    
    ax_epoch.spines['left'].set_visible(True)
    ax_epoch.spines['left'].set_linewidth(EPOCH_LINE_WIDTH_THICK)
    ax_epoch.spines['right'].set_visible(True)
    ax_epoch.spines['right'].set_linewidth(EPOCH_LINE_WIDTH_THICK)
    ax_epoch.spines['top'].set_linewidth(EPOCH_LINE_WIDTH_THICK)
    ax_epoch.spines['bottom'].set_linewidth(EPOCH_LINE_WIDTH_THICK)
    
    ax_epoch.axhline(y=0.5, color=EPOCH_LINE_COLOR, linewidth=EPOCH_LINE_WIDTH_THIN)
    
    # Official Geological Boundaries
    ep_late_holo = 11700     
    ep_early_mid = 774000    
    ep_mid_late = 129000     
    hol_early_mid = 8200     
    hol_mid_late = 4200      
    
    # Background Colors for the epoch box itself
    ax_epoch.axvspan(xmin=x_max, xmax=ep_late_holo, facecolor='#EADDC8', zorder=1) 
    ax_epoch.axvspan(xmin=ep_late_holo, xmax=x_min, facecolor='#F3EBDD', zorder=1) 

    ax_epoch.axvline(x=ep_late_holo, color=EPOCH_LINE_COLOR, linewidth=EPOCH_LINE_WIDTH_THICK, zorder=2)

    ax_epoch.axvline(x=ep_early_mid, ymin=0, ymax=0.5, color=EPOCH_LINE_COLOR, linewidth=EPOCH_LINE_WIDTH_THIN, zorder=2)
    ax_epoch.axvline(x=ep_mid_late, ymin=0, ymax=0.5, color=EPOCH_LINE_COLOR, linewidth=EPOCH_LINE_WIDTH_THIN, zorder=2)
    ax_epoch.axvline(x=hol_early_mid, ymin=0, ymax=0.5, color=EPOCH_LINE_COLOR, linewidth=EPOCH_LINE_WIDTH_THIN, zorder=2)

    def log_mid(x1, x2):
        return np.sqrt(max(x1, x_min) * min(x2, x_max))

    # --- TOP TIER TEXT ---
    ax_epoch.text(log_mid(x_max, ep_late_holo) * SHIFT_PLEISTOCENE, 0.75, 'Pleistocene', ha='center', va='center', fontsize=EPOCH_FONT_SIZE + 1, zorder=3)
    ax_epoch.text(log_mid(ep_late_holo, x_min) * SHIFT_HOLOCENE, 0.75, 'Holocene', ha='center', va='center', fontsize=EPOCH_FONT_SIZE + 1, zorder=3)

    # --- BOTTOM TIER TEXT ---
    ax_epoch.text(log_mid(x_max, ep_early_mid) * SHIFT_EARLY_PLEISTOCENE, 0.25, 'E', ha='center', va='center', fontsize=EPOCH_FONT_SIZE, zorder=3)
    ax_epoch.text(log_mid(ep_early_mid, ep_mid_late) * SHIFT_MIDDLE_PLEISTOCENE, 0.25, 'Middle', ha='center', va='center', fontsize=EPOCH_FONT_SIZE, zorder=3)
    ax_epoch.text(log_mid(ep_mid_late, ep_late_holo) * SHIFT_LATE_PLEISTOCENE, 0.25, 'Late', ha='center', va='center', fontsize=EPOCH_FONT_SIZE, zorder=3)

    # Holocene (E and M only)
    ax_epoch.text(log_mid(ep_late_holo, hol_early_mid) * SHIFT_EARLY_HOLOCENE, 0.25, 'E', ha='center', va='center', fontsize=EPOCH_FONT_SIZE, zorder=3)
    ax_epoch.text(log_mid(hol_early_mid, x_min) * SHIFT_MIDDLE_HOLOCENE, 0.25, 'M', ha='center', va='center', fontsize=EPOCH_FONT_SIZE, zorder=3)

    # ==========================================
    # 2. SECOND PANE: Temperature
    # ==========================================
    ax_temp.plot(times_temp, temperatures, color=TEMP_LINE_COLOR, linewidth=TEMP_LINE_WIDTH, zorder=5)
    ax_temp.set_ylim(TEMP_MIN, TEMP_MAX)
    
    ax_temp.set_yticks([3, 0, -3, -6])
    ax_temp.set_ylabel("Temperature\n(°C)", fontsize=FONT_SIZE, family='sans-serif', labelpad=LABEL_PAD_Y_TEMP) 

    # ==========================================
    # 3. THIRD PANE: Sea Level
    # ==========================================
    ax_sea.plot(times_sea_level, sea_levels, color=SEA_LEVEL_LINE_COLOR, linestyle='-', linewidth=SEA_LEVEL_LINE_WIDTH, label="Sea Level (m)")
    ax_sea.set_ylim(SEA_LEVEL_MIN, SEA_LEVEL_MAX)
    
    ax_sea.set_yticks([0, -40, -80, -120])
    ax_sea.set_ylabel("Sea Level\n(m)", fontsize=FONT_SIZE, family='sans-serif', labelpad=LABEL_PAD_Y_SEA)

    # ==========================================
    # 4. BOTTOM PANE: PSMC
    # ==========================================
    legend_handles = []
    with open(psmc_file_list, "r") as psmc_file:
        for line in psmc_file:
            line = line.rstrip()
            if not line: continue
            PSMC_RESULT, line_color, line_label = line.split("\t")
            
            runs = psmc_fun(PSMC_RESULT, bin_size, mutation_rate, generation_time)
            
            if runs:
                main_times, main_sizes = runs[0]
                ax_psmc.step(main_times, main_sizes, where='post', linestyle='-', color=line_color,
                             linewidth=PSMC_LINE_WIDTH, zorder=5)
                
                for b_times, b_sizes in runs[1:]:
                    ax_psmc.step(b_times, b_sizes, where='post', linestyle='-', color=line_color,
                                 linewidth=BOOTSTRAP_LINE_WIDTH, alpha=BOOTSTRAP_ALPHA, zorder=4)
                
                legend_copy, = ax_psmc.plot([], [], color=line_color, linewidth=2.5, label=line_label)
                legend_handles.append(legend_copy)

    ax_psmc.set_ylabel("Effective pop size\n(x10⁴)", fontsize=FONT_SIZE, family='sans-serif', labelpad=LABEL_PAD_Y_PSMC) 
    ax_psmc.set_xlabel("Years before present", fontsize=FONT_SIZE, family='sans-serif', labelpad=LABEL_PAD_X_PSMC)

    ax_psmc.set_yticks([0, 10000, 20000, 30000, 40000, 50000])

    # --- Axes formatting & VERTICAL COLD SHADING ---
    for ax in [ax_epoch, ax_temp, ax_sea, ax_psmc]:
        
        # 🚨 FIX: edgecolor='none' and linewidth=0 completely removes the black lines 🚨
        if ax != ax_epoch and len(temperatures) > 0:
            ax.fill_between(times_temp, 0, 1, where=(temperatures < TEMP_COLD_THRESHOLD),
                            transform=ax.get_xaxis_transform(), facecolor=COLD_SHADE_COLOR, 
                            edgecolor='none', linewidth=0, zorder=0)

        ax.tick_params(axis='both', labelsize=FONT_SIZE)
        
        if ax == ax_epoch:
            ax.tick_params(axis='y', which='both', left=False, right=False)
        else:
            ax.tick_params(axis='y', which='both', left=True, right=False, length=4) 
        
        if ax != ax_psmc:
            ax.tick_params(axis='x', which='both', bottom=False, top=False)
        
        if ax != ax_epoch:
            LGM_position = 26500
            ax.axvline(x=LGM_position, color=MARKER_LINE_COLOR, linestyle=MARKER_LINE_STYLE, linewidth=MARKER_LINE_WIDTH, zorder=2)
        
        if span_max != 0 or span_min != 0:
            ax.axvspan(xmin=span_min, xmax=span_max, facecolor=span_color, alpha=span_alpha, zorder=0)

    ax_psmc.tick_params(axis='x', which='major', bottom=True, top=False, length=4) 
    ax_psmc.tick_params(axis='x', which='minor', bottom=False, top=False)            
    ax_psmc.xaxis.set_minor_locator(plt.NullLocator())                               

    ax_psmc.text(LGM_position, y_min - 2000, 'LGM', ha='center', va='center', 
                 fontsize=MARKER_FONT_SIZE, color=MARKER_FONT_COLOR, 
                 fontweight=MARKER_FONT_WEIGHT, fontstyle=MARKER_FONT_STYLE) 

    ax_psmc.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2), useOffset=False)
    ax_psmc.yaxis.get_offset_text().set_visible(False)
    
    ax_psmc.set_xlim(x_min, x_max)
    ax_psmc.set_ylim(y_min, y_max)
    ax_psmc.set_xscale('log')
    ax_psmc.invert_xaxis()

    # --- Legend ---
    fig.legend(handles=legend_handles, loc='lower center', bbox_to_anchor=(LEGEND_POS_X, LEGEND_POS_Y),
                ncol=len(legend_handles), fontsize=FONT_SIZE, frameon=False)

    # === Export ===
    plt.savefig("psmc_temperature_sea_level_plot.pdf", dpi=DPI_EXPORT, bbox_inches="tight", pad_inches=0.3)
    buf = io.BytesIO()
    fig.savefig(buf, format='tiff', dpi=DPI_EXPORT, pil_kwargs={"compression": "tiff_lzw"})
    buf.seek(0)
    tiff_img = Image.open(buf)
    tiff_img.save("psmc_temperature_sea_level_plot.tiff", compression="tiff_lzw")
    buf.close()
    
    plt.show()

# --- PSMC parsing function ---
def psmc_fun(filename, size, mutation, generation_time):
    with open(filename, "r") as f:
        result = f.read()
    
    blocks = result.split('//\n')
    runs = []
    
    for i in range(len(blocks) - 1):
        current_block = blocks[i]
        
        is_final = False
        if i == len(blocks) - 2:
            is_final = True
        else:
            next_block = blocks[i+1].lstrip()
            if next_block.startswith("RD\t0") or next_block.startswith("RD 0") or next_block.startswith("CC"):
                is_final = True
                
        if is_final:
            time_windows = []
            estimated_lambdas = []
            theta = None
            
            for line in current_block.split('\n'):
                if line.startswith("RS"):
                    parts = line.split('\t')
                    if len(parts) >= 4:
                        time_windows.append(float(parts[2]))
                        estimated_lambdas.append(float(parts[3]))
                elif line.startswith("PA"):
                    try:
                        pa_str = line.split('PA\t')[-1]
                        pa_parts = pa_str.split(' ')
                        theta = float(pa_parts[1])
                    except:
                        parts = line.split()
                        theta = float(parts[2])
            
            if theta is not None and len(time_windows) > 0:
                N0 = theta / (4 * mutation) / size
                times = [generation_time * 2 * N0 * t for t in time_windows]
                sizes = [N0 * l for l in estimated_lambdas]
                
                raw_dict = dict(zip(times, sizes))
                false_result = sizes[-1]
                
                valid_times = []
                valid_sizes = []
                for t, s in zip(times, sizes):
                    if s != false_result:
                        valid_times.append(t)
                        valid_sizes.append(s)
                
                runs.append((valid_times, valid_sizes))
                
    return runs

if __name__ == '__main__':
    plot_psmc_with_temperature_and_sea_level()