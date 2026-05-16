import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from PIL import Image
import io

# === Journal Guidelines Settings ===
ONE_MM = 1 / 25.4
FIG_WIDTH_MM = 175      # double column width
FIG_HEIGHT_MM = 130
DPI = 600               # high-res for TIFF

FIG_WIDTH_IN = FIG_WIDTH_MM * ONE_MM
FIG_HEIGHT_IN = FIG_HEIGHT_MM * ONE_MM

# Enforce Arial font and line widths per journal
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
rcParams['lines.linewidth'] = 0.5

# === Load data ===
df = pd.read_csv('enrichment_data.txt', delimiter='\t')
df.columns = ['Genome', 'Term', 'Count', 'p-value']

# --- START OF NEW CODE TO FORMAT TERMS ---
def format_term(term):
    # Convert the entire string to lowercase first
    term = term.lower()
    
    # Capitalize the first letter of the first word only
    term = term.capitalize()
    
    # Handle 'dna' and 'rna' specifically within the string
    # We re-uppercase these acronyms wherever they might appear
    term = term.replace('dna', 'DNA').replace('rna', 'RNA')
    
    return term

df['Term'] = df['Term'].apply(format_term)
# --- END OF NEW CODE TO FORMAT TERMS ---

df['p-value'] = pd.to_numeric(df['p-value'].str.replace(',', '.'), errors='coerce')
df['p-value'] = 10 ** df['p-value']
df = df[df['Genome'] != 'S. quin.']
df = df.sort_values(by='p-value')

# Map terms to y-axis
terms_sorted = df['Term'].unique()
term_to_y = {term: i for i, term in enumerate(reversed(terms_sorted))}

# === Plot ===
fig, ax = plt.subplots(figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN))

# Define genome colors (MATCHING YOUR PSMC PLOT)
colors = {
    "CND": "darkblue",    # Matches your PSMC file
    "CNY": "purple",      # Matches your PSMC file
    "SA": "darkorange",   # Matches your PSMC file
    "USA": "darkgreen",   # Matches your PSMC file
    "AUS": "#FF1493",  # Deep Pink (Dalian - distinct from Purple)
}

for genome in df['Genome'].unique():
    subset = df[df['Genome'] == genome].copy()
    subset['y-axis'] = subset['Term'].map(term_to_y)
    ax.scatter(
        subset['p-value'],
        subset['y-axis'],
        s=subset['Count']*10,
        alpha=0.6,
        edgecolors="w",
        linewidth=0.5,
        color=colors.get(genome, 'gray'),
        zorder=3
    )

ax.margins(x=0.09)
ax.set_xlabel('log10(P-value)', fontsize=9)  # legend & axis fonts = 9 pt

# ---> CHANGED HERE: Added labelpad=15 to push the title left
ax.set_ylabel('Ontology terms', fontsize=10, labelpad=15)

ax.set_xscale('log')
ax.set_xlim(ax.get_xlim()[::-1])
ax.set_yticks(range(len(terms_sorted)))
ax.set_yticklabels(reversed(terms_sorted), fontsize=9)

# Remove minor ticks
ax.minorticks_off()
ax.tick_params(axis='x', which='minor', bottom=False)

ax.grid(True, linestyle='--', color='lightgray', alpha=0.7, linewidth=0.5)

# Legends (gene count)
size_legend = [3, 6, 12, 16]
legend_sizes = [s*10 for s in size_legend]
size_handles = [plt.scatter([], [], s=size, color='gray', alpha=0.6, edgecolors="w", linewidth=0.5)
                for size in legend_sizes]
legend2 = plt.legend(
    handles=size_handles,
    labels=[str(s) for s in size_legend],
    title="Gene count",
    loc="lower left",
    bbox_to_anchor=(1, -0.013),
    handletextpad=0.8,  # ---> CHANGED HERE: Decreased from 1.5 to shrink width
    labelspacing=0.6,
    borderpad=0.8,      # ---> CHANGED HERE: Decreased from 1.0 to shrink width
    frameon=True,
    fontsize=9,
    title_fontsize=9
)
plt.gca().add_artist(legend2)

# Legends (genomes)
color_handles = [plt.scatter([], [], s=100, color=colors[g], alpha=0.6, edgecolors="w", linewidth=0.5)
                 for g in colors]
legend1 = ax.legend(
    handles=color_handles,
    labels=list(colors.keys()),
    title="Genome",
    loc="lower left",
    bbox_to_anchor=(1, 0.25),
    handletextpad=0.8,  # ---> CHANGED HERE: Decreased from 1.0 to shrink width
    labelspacing=0.4,
    borderpad=0.8,      # ---> CHANGED HERE: Decreased from 1.0 to shrink width
    frameon=True,
    fontsize=9,
    title_fontsize=9
)
ax.add_artist(legend1)

plt.xticks(rotation=0, fontsize=9)
plt.yticks(fontsize=9)
plt.tight_layout(rect=[0, 0, 0.80, 1])

# === Save outputs ===
fig.savefig('ontology_vertical_plot.svg', bbox_inches='tight', pad_inches=0.2)
fig.savefig('ontology_vertical_plot.pdf', bbox_inches='tight', pad_inches=0.2)
buf = io.BytesIO()
fig.savefig(buf, format='tiff', dpi=DPI, pil_kwargs={"compression": "tiff_lzw"})
buf.seek(0)
tiff_img = Image.open(buf)
tiff_img.save('ontology_vertical_plot.tiff', compression='tiff_lzw')
buf.close()

plt.show()