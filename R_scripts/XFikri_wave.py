#Fikir waves


################################################
#Peak period

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load CSV
file_path = r'C:\Users\gerar\OneDrive\1_Work\3_Results\11 Allee effects\Fikri_data\peak_period_canyon.csv'
df_peak = pd.read_csv(file_path)

# Convert Date column to datetime and filter for years > 2005
df_peak['date'] = pd.to_datetime(df_peak['date'], errors='coerce')
df_peak = df_peak[df_peak['date'].dt.year > 2005]

# Compute quantiles
quantiles_peak = df_peak.quantile([0.05, 0.50, 0.95])
quantiles_peak
# Plot
# Plot


###################################################
#Wave heights



# Load CSV
file_path1 = r'C:\Users\gerar\OneDrive\1_Work\3_Results\11 Allee effects\Fikri_data\wave_heights_canyon.csv'
df_wave = pd.read_csv(file_path1)

# Convert Date column to datetime and filter for years > 2005
df_wave['date'] = pd.to_datetime(df_wave['date'], errors='coerce')
df_wave = df_wave[df_wave['date'].dt.year > 2005]

# Compute quantiles
quantiles_wave = df_wave.quantile([0.05, 0.50, 0.95])
quantiles_wave
# Plot
# Plot



#############################
#Plots

#Period

fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(8, 10), sharex=True)

# Peak Period
for col in df_peak.columns:
    if col != 'date':  
        axes[0].scatter(df_peak['date'], df_peak[col], alpha=0.2, color='grey', s=15, label=f'Raw Data ({col})')

axes[0].axhline(quantiles_peak.loc[0.05].values, color='#4682B4', linestyle='dashed', label='5th percentile')
axes[0].axhline(quantiles_peak.loc[0.50].values, color='#D62728', linestyle='dashed', label='50th percentile')
axes[0].axhline(quantiles_peak.loc[0.95].values, color='#4682B4', linestyle='dashed', label='95th percentile')
axes[0].set_ylabel('Peak period (s)')
#axes[0].set_title('Peak Period Canyon')
axes[0].legend()

# Wave Heights
for col in df_wave.columns:
    if col != 'date':  
        axes[1].scatter(df_wave['date'], df_wave[col], alpha=0.2, color='grey', s=15, label=f'Raw Data ({col})')

axes[1].axhline(quantiles_wave.loc[0.05].values, color='#4682B4', linestyle='dashed', label='5th percentile')
axes[1].axhline(quantiles_wave.loc[0.50].values, color='#D62728', linestyle='dashed', label='50th percentile')
axes[1].axhline(quantiles_wave.loc[0.95].values, color='#4682B4', linestyle='dashed', label='95th percentile')
axes[1].set_ylabel('Significant wave heights (m)')
$axes[1].set_title('Significant Wave Heights')
axes[1].set_xlabel('Date')
fig.tight_layout()
plt.show()

