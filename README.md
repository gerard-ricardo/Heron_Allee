# 2022 Heron field experiment

The project aimed to assess small population effects (Allee) on coral reproduction. Corals along a 125 m section of reef slope were georeferenced, sequenced, and monitored for spawning. Coral fertilisation was determined in such a way that we knew that 'dam' of each embryo but hoped to use parentage assignment to identify 'sires'. 

## Workflow

Download the whole project and open the project file. In 'R_scripts' run '1_2022heron_seq_load_filt.R', which loads and filters the SNP data. At present, this creates a dart genlight object 'data_gl_filtered' and a genind object 'data_genind' for the whole experiment (adults  + larvae). A separate dart genlight object 'data_gl_filtered_adult' and genind 'data_genind_adult' is created for only adults for most of the pop gen analyses.


```bash
# Clone the repository
git clone https://github.com/yourusername/your-repository.git

# Navigate to the project directory
cd your-repository

# Install dependencies (if any)
