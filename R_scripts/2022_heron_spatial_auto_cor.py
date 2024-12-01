
#Moran's I Values (r):
#Moran's I is a measure of spatial autocorrelation. It quantifies how similar individuals are genetically, 
#relative to their geographic proximity. Pos =  pos corrlation, neg = neg correlation




import pandas as pd
import numpy as np
from scipy.spatial import distance_matrix
from sklearn.metrics.pairwise import euclidean_distances
import matplotlib.pyplot as plt

# Load the metadata
meta_platy = pd.read_csv("./data/meta_platy_ordered.csv")

# Load the SNP data
data_gl = pd.read_csv("./data/Report_DPlatyg23-7805_SNP_2 - Copy corrected.csv", skiprows=6)
data_gl1 = data_gl.columns[24:103]   #extract data

# Filter out eggs and larvae to keep only adults
meta_platy_adults = meta_platy[meta_platy["stage"] == "adults"].copy()
data_gl_adults = data_gl.loc[meta_platy_adults.index].copy()

# Preprocess the SNP data: convert non-numeric values to NaN and then to a numeric format
data_gl_adults.replace("-", np.nan, inplace=True)
data_gl_adults = data_gl_adults.apply(pd.to_numeric, errors='coerce')

# Calculate call rate for individuals (percentage of non-NaN values from columns 25 to 103)
individual_columns = data_gl_adults.columns[24:103]
ind_call_rate = data_gl_adults[individual_columns].notna().mean(axis=1)
# Filter individuals with call rate above the threshold (0.5)
ind_filtered = ind_call_rate[ind_call_rate >= 0.5].index
data_gl_adults = data_gl_adults.loc[ind_filtered]
meta_platy_adults = meta_platy_adults.loc[ind_filtered]

# Extract the call rate for loci
loc_call_rate = data_gl_adults.loc["CallRate"]

# Filter loci with call rate above the threshold (0.7)
loc_filtered = loc_call_rate[loc_call_rate >= 0.7].index
filtered_columns = data_gl_adults.columns[loc_filtered]

# Ensure coordinates match and are in the correct order
coordinates = meta_platy_adults[["latitude", "longitude"]].values

# Calculate the genetic distance matrix using Euclidean distances
#subset for 
genetic_dist_matrix = distance_matrix(data_gl_adults.fillna(0).values, data_gl_adults.fillna(0).values)

# Calculate the Euclidean distance matrix based on coordinates
euclidean_dist_matrix = euclidean_distances(coordinates, coordinates)

# Perform spatial autocorrelation analysis using Moran's I
def calculate_morans_i(gen_dist, geo_dist, bins=10):
    # Create distance bins
    bin_edges = np.linspace(geo_dist.min(), geo_dist.max(), bins + 1)
    bin_indices = np.digitize(geo_dist, bin_edges) - 1
    
    r_values = []
    for i in range(bins):
        # Select pairs in the current bin
        pairs = np.where(bin_indices == i)
        if len(pairs[0]) == 0:
            r_values.append(np.nan)
            continue
        # Calculate Moran's I
        gen_mean = gen_dist[pairs].mean()
        geo_mean = geo_dist[pairs].mean()
        num = np.sum((gen_dist[pairs] - gen_mean) * (geo_dist[pairs] - geo_mean))
        den = np.sum((gen_dist[pairs] - gen_mean)**2) * np.sum((geo_dist[pairs] - geo_mean)**2)
        r = num / np.sqrt(den)
        r_values.append(r)
    
    return pd.DataFrame({"dist": bin_edges[:-1], "r": r_values})

# Calculate Moran's I
spatial_autocor_results = calculate_morans_i(genetic_dist_matrix, euclidean_dist_matrix, bins=10)

# Display the results
print(spatial_autocor_results)

# Plot the results
plt.plot(spatial_autocor_results["dist"], spatial_autocor_results["r"], marker='o')
plt.xlabel("Distance Classes")
plt.ylabel("Autocorrelation Coefficient (r)")
plt.title("Spatial Autocorrelation Analysis")
plt.show()
