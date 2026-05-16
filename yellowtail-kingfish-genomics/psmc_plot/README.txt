# PSMC Demographic History Visualization

## Description
This repository contains a Python 3 pipeline to visualize Pairwise Sequentially Markovian Coalescent (PSMC) demographic histories overlaid with paleoclimate data. 
The script generates a multi-panel figure incorporating geological epochs, historical temperature fluctuations, sea-level changes, and effective population size (Ne) trajectories.

## Prerequisites
Dependencies can be installed via standard package managers:
`pip install numpy matplotlib click Pillow`

## Input Data Requirements
You need three tab-separated text files to run this script:

### 1. PSMC Configuration File (`psmc_config.txt`)
A tab-separated file listing the path to the PSMC output, the desired line color, and the population/lineage label.
*Format:* `[filename.psmc] \t [color] \t [label]`
*Example:*
```text
CND.psmc    darkblue    CND
SA.psmc     darkorange  SA
USA.psmc    darkgreen   USA
CNY.psmc    purple      CNY

## 2. Paleotemperature Data (paleo_temperature.txt)
Time (kyr BP)    50%
1                -0.17
2                -0.09
3                0.03

## 3. Sea Level Data (paleo_sea_level.txt)
Age (ka)    PC1
0           8.49
1           7.63
2           4.01

## Usage
Execute the script via the command line. You must specify your generation time and mutation rate.Here is an example command utilizing a generation time of 3 years and a mutation rate of 2.5x10^-8: 
python plot_psmc_demography.py \
  --psmc_file_list psmc_config.txt \
  --temperature_file paleo_temperature.txt \
  --sea_level_file paleo_sea_level.txt \
  --generation_time 3 \
  --mutation_rate 2.5e-08 \
  --x_min 7e3 \
  --x_max 1e6 \
  --y_max 60e3 \
  --span_color yellow

## Outputs
The script generates a multi-panel figure exported in high resolution for publication:
psmc_temperature_sea_level_plot.pdf
psmc_temperature_sea_level_plot.tiff (600 DPI, LZW compressed)

## Customization
A # CUSTOMIZATION HUB is located at the top of the plot_psmc_demography.py script. Users can easily modify panel size ratios, font styles, vertical event markers (e.g., LGM), and export resolutions without altering the core plotting logic.

## References
Li, H., & Durbin, R. (2011). Inference of Human Population History From Whole Genome Sequence of A Single Individual. *Nature*, 475(7357), 493-496.