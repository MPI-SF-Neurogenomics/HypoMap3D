# 00 setup

First time? Congrats, you are starting in the right place.

Please do, however, read the repository README first, specifically the section on [HypoMap3D Development](../../README.md/#hypomap3d-development).

## Setting up your machine for running the HypoMap3D Pipeline

For this walkthough, meet John. You will need to replace John's filepaths below with your own directory of choice.

### Step 1: Cloning the HypoMap3D Repostiory

```bash
cd /scratch/John # replace with any directory of your choosing
git clone https://github.com/MPI-SF-Neurogenomics/HypoMap3D.git
```

This will create a new folder called 'HypoMap3D'. In John's case, `/scratch/John/HypoMap3D`.

### Step 2: Create a data directory

Create a data directory in the same directory as the HypoMap3D folder.
```bash
cd /scratch/John
mkdir HypoMap3D-data
```
So John's filesystem now looks like:
```bash
/scratch/John/
├── HypoMap3D/
└── HypoMap3D-data/
```

### Step 3: Pull requisite docker

We use containers to ensure fixed environments and a stable work environment.
For more details please refere to [containers.md](../../containers.md).

For this walkthough we will pull these images to the data directory but as long as the path is known, they can be put anywhere.

```bash
cd /scratch/John/HypoMap3D-data
mkdir containers

cd /scratch/John/HypoMap3D/atlas-construction/00_setup
bash 01_pull_containers.sh /scratch/John/HypoMap3D-data/containers
```

### Step 4: Updating the config.env

All paths defined through the pipeline scripts to follow are defined relative to the data directory defined in the `config.env`.

Firstly, navigate to the config.env:
```bash
cd /scratch/John/HypoMap3D
nano config.env # edit the config
```

**You will need to update the path to the newly created data directory, as well as paths to the containers**

In John's case, his config.env may look like this:
```bash
# Path to initialize the root directory for all data stored 
DATADIR="/scratch/John/HypoMap3D-data"

# Image paths
UV_PATH="/scratch/John/HypoMap3D-data/containers/uv_r_cuda_0.2.1.sif"
CELLRANGER_PATH="/scratch/John/HypoMap3D-data/containers/cellranger_9.0.0.sif"
STAR_PATH="/scratch/John/HypoMap3D-data/containers/STAR_2.7.11b.sif"
CELLBENDER_PATH="/scratch/John/HypoMap3D-data/containers/cellbender_image.sif"

# Project-wide defaults for matplotlib plotting functions
MPLSTYLE="/scratch/John/HypoMap3D/resources/figures.mplstyle"

# SLURM variables. 
# Required if you intend to use SLURM for a subset of our scripts.
# SLURM email -> get emails for ending/failing SLURM jobs
SLURM_EMAIL="John@scammer.anon"
CPU_PARTITION="cpu-node-1"
GPU_PARTITION="gpu-node-3"
```

### Step 5: Validate the config.env

A simple script to check that the `config.env` has been filled out properly can be used.

```bash
cd /scratch/John/HypoMap3D/atlas-construction/00_setup
bash 02_validate_config.sh
```

If these checks pass, your directory structure is ready to run the pipeline.

### Step 6: Restoring development environment

Information on restoring the development environment can be found on the repo README (see [here](../../README.md/#hypomap3d-development)).

A simplified script has been provided, however.

```bash
cd /scratch/John/HypoMap3D
bash atlas-construction/00_setup/03_restore_env.sh
```

At this stage, John (and you) are ready to run the full pipeline.

### Step 7: Adding environment as a Jupyter kernel (OPTIONAL)

You can add the newly created `.venv` as a jupyter kernel. This is entirely optional.

```bash
cd /scratch/John/HypoMap3D
bash atlas-construction/00_setup/04_add_kernel.sh
```

### John's Final Structure

Your directory structure should look something like this.
```bash
├── HypoMap3D
│   ├── .venv
│   ├── atlas-construction
│   ├── ... # rest of repository
│   └── config.env
└── HypoMap3D-data
    └── containers
```