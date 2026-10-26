# Spatial Objects

Understanding the spatial localization of cells within the hypothalamus is essential both for the classification of cell types within the hypothalamus but also for understanding the interaction between them. HypoMap3D integrates both spatial and single cell transcriptomics to further explore cellular heterogeneity in the region. Raw spatial transcriptomics are combined to form an **AnnData** object and combined cluster probabilties are integrated to form a smooth voxel grids as **Zarr** objects.

<!-- ======================================================================= -->
<!--                               Raw Spatial                               -->
<!-- ======================================================================= -->

## Spatial AnnData

Single cell resolution spatial transcriptomics from 4 major datasets were concatenated and harmonized. Cells within a defined bounding box around the hypothalamus were kept leaving ~1.4 mil cells of which ~400k reside within the hypothalamus itself. The label transfer probablities of cell-types *L1-5* are included as fuzzy sparse matrices in the `.obsm`.

### Structure

```bash
AnnData object with n_obs × n_vars = 1,401,313 × 56,939
    X: raw count matrix
    obs: 'Dataset', 'Sample', 'Author', 'Technology', 'Donor', 'Sex', 'RegionID', 'RegionAcronym', 'RegionName', 'total_counts', 'total_hvg_counts', 'total_genes', 'IsHY', 'L1_descr_max_p', 'L1_abc_like_max_p', 'L1_max_p', 'L2_descr_max_p', 'L2_abc_like_max_p', 'L2_max_p', 'L3_descr_max_p', 'L3_abc_like_max_p', 'L3_max_p', 'L4_descr_max_p', 'L4_abc_like_max_p', 'L4_max_p', 'L5_descr_max_p', 'L5_abc_like_max_p', 'L5_max_p'
    var: 'gene_stable_ids', 'IsHVG'
    uns: 'L1_ID', 'L1_abc_like', 'L1_descr', 'L2_ID', 'L2_abc_like', 'L2_descr', 'L3_ID', 'L3_abc_like', 'L3_descr', 'L4_ID', 'L4_abc_like', 'L4_descr', 'L5_ID', 'L5_abc_like', 'L5_descr', 'flags', 'pathmat'
    obsm: 'P_L1', 'P_L2', 'P_L3', 'P_L4', 'P_L5', 'X_embed', 'X_umap', 'X_umap_sc', 'coords_ccf_mm', 'coords_ccf_mm_wrap', 'coords_source', 'coords_voxels_25um', 'coords_voxels_25um_bounded'
    varm: 'IsCaptured'
```

### Access

```{caution}
Public data access to `hm3d_st_all.h5ad` is currently unavailable as HypoMap3D is not yet published.
```

### Using/Loading

The full anndata can be loaded using the `AnnData` python package or via the `ScanPy` wrapper

```python
import anndata as ad # via AnnData
scdata = ad.read_h5ad('hm3d_st_all.h5ad')

import scanpy as sc # via ScanPy
scdata = sc.read_h5ad('hm3d_st_all.h5ad')
```

### Subsetting to the Hypothalamus
You can subset to only hypothalamic cells using the `IsHY` mask in the `obs` dataframe.

```python
scdata_hy = scdata[scdata.obs['IsHY'], :]
```


<!-- ======================================================================= -->
<!--                              Voxel Grid(s)                              -->
<!-- ======================================================================= -->

## Voxel Grids Zarr

The smoothed spatial voxel grids are available in the Zarr file format.

### Structure

```bash
.
├── L1
│   ├── abc-like (6,) StringDType()
│   ├── descr (6,) StringDType()
│   ├── flags (6,) StringDType()
│   ├── ids (6,) StringDType()
│   ├── mass (6, 168, 144, 107) float32
│   ├── prob (6, 168, 144, 107) float32
│   └── rel_mass (6, 168, 144, 107) float32
├── L2
│   ├── abc-like (82,) StringDType()
│   ├── descr (82,) StringDType()
│   ├── flags (82,) StringDType()
│   ├── ids (82,) StringDType()
│   ├── mass (82, 168, 144, 107) float32
│   ├── prob (82, 168, 144, 107) float32
│   └── rel_mass (82, 168, 144, 107) float32
├── L3
│   ├── abc-like (274,) StringDType()
│   ├── descr (274,) StringDType()
│   ├── flags (274,) StringDType()
│   ├── ids (274,) StringDType()
│   ├── mass (274, 168, 144, 107) float32
│   ├── prob (274, 168, 144, 107) float32
│   └── rel_mass (274, 168, 144, 107) float32
├── L4
│   ├── abc-like (716,) StringDType()
│   ├── descr (716,) StringDType()
│   ├── flags (716,) StringDType()
│   ├── ids (716,) StringDType()
│   ├── mass (716, 168, 144, 107) float32
│   ├── prob (716, 168, 144, 107) float32
│   └── rel_mass (716, 168, 144, 107) float32
├── L5
│   ├── abc-like (932,) StringDType()
│   ├── descr (932,) StringDType()
│   ├── flags (932,) StringDType()
│   ├── ids (932,) StringDType()
│   ├── mass (932, 168, 144, 107) float32
│   ├── prob (932, 168, 144, 107) float32
│   └── rel_mass (932, 168, 144, 107) float32
├── bbox
│   ├── ccf_mm
│   │   ├── max (3,) float32
│   │   └── min (3,) float32
│   ├── shape (3,) float32
│   └── voxel_25_um
│       ├── max (3,) float32
│       └── min (3,) float32
├── counts (168, 144, 107) float32
├── evidence (168, 144, 107) float32
├── latent_embedding (45, 168, 144, 107) float32
└── region
    ├── acronyms (176,) StringDType()
    ├── ids (176,) int32
    ├── region_labels (168, 144, 107) uint32
    └── region_masks (176, 168, 144, 107) bool
```

### Access

```{caution}
Public data access to `hm3d_st_voxelized.zarr` is currently unavailable as HypoMap3D is not yet published.

```

### Using/Loading

The voxel grids can be loaded using the `Zarr` python package with **version >3.1**
```python
import zarr

voxel_store = zarr.open('hm3d_st_voxelized.zarr')

# Reading slice 50 of L1 smoothed probabilities for cluster `Gb1-3 HY Gaba`
clusters = voxel_store['L1'].attrs['descr']
idx = clusters.index('Gb1-3 HY Gaba')
view = voxel_store['L1'][idx, 50, :, :]
```

```{seealso}
For more information on handling the zarr voxel grids, please refer to
the {doc}`SpatialSlicePlotting Tutorial <../tutorials/SpatialSlicePlotting>`
```

## Cluster Region Assignment

The relative mass contribution of hypothalamic sub-nuclei to each cluster is provided as a CSV table.

### Structure
The first row of the table is as follows:

|celltype_id|celltype_descr|celltype_abc_like|flag      |TMv         |IF   |PT               |AAA               |IGL               |PVa      |PVH        |ACB      |MT |IMD|PVHd              |AD |LT           |INC          |ADP         |VL         |AHN        |SEZ          |PVi        |PVp         |MRN        |V3          |PVpo        |PVT         |RCH         |LGv         |RE          |RH           |LHA         |RL           |LM           |RN          |APN          |ARH        |LPO        |AV          |LSr        |NLOT1       |RT        |AVP         |LSv        |NLOT2       |AVPV        |SCH         |BAC          |BA           |MA          |SF          |MB          |SubG        |BMAa        |ASO         |SH           |BMAp         |SI        |SBPV      |BST        |PST         |MD         |PSTN        |SMT        |SNc         |SNr         |CA1         |SO           |MEA        |SPFm        |SPFp        |CA2          |MEPO         |CA3        |STN         |STR        |MPN        |MPO        |SUM         |CEAc         |TH          |CEAl         |CEAm       |MS         |CL           |TRS        |ND           |NDB        |CM           |SPA         |TU          |VAL         |DG-sg        |COAa       |COApl        |COApm       |CP       |VM          |VLPO        |VMH         |OLF         |CTXsp       |VPL         |VPLpc       |MMme         |VPM         |VPMpc       |VTA         |OT         |OV           |P           |PA         |PAA          |PAG          |ZI         |PAL        |FF          |DMH         |IIIn        |IIn          |PCN         |PD           |PF           |PH          |EPd         |PIR         |EW          |PMd         |FS         |PMv        |PO          |GPe        |POL         |GPi         |PP           |PR          |MGv         |MGm          |HPF          |AMd         |HY        |AMv         |IA           |PS          |IAD         |IAM          |TMd          |NLOT3       |ME        |DG-mo       |DG-po        |MA3         |Eth         |Xi          |PIL         |PoT          |IntG        |VMPO         |PeF         |HATA        |MMl         |MMm          |MMp          |MMd          |PN           |IPR          |IPC        |IPL          |IPI          |IPDM         |IPDL         |IPRL         |
|-----------|--------------|-----------------|----------|------------|-----|-----------------|------------------|------------------|---------|-----------|---------|---|---|------------------|---|-------------|-------------|------------|-----------|-----------|-------------|-----------|------------|-----------|------------|------------|------------|------------|------------|------------|-------------|------------|-------------|-------------|------------|-------------|-----------|-----------|------------|-----------|------------|----------|------------|-----------|------------|------------|------------|-------------|-------------|------------|------------|------------|------------|------------|------------|-------------|-------------|----------|----------|-----------|------------|-----------|------------|-----------|------------|------------|------------|-------------|-----------|------------|------------|-------------|-------------|-----------|------------|-----------|-----------|-----------|------------|-------------|------------|-------------|-----------|-----------|-------------|-----------|-------------|-----------|-------------|------------|------------|------------|-------------|-----------|-------------|------------|---------|------------|------------|------------|------------|------------|------------|------------|-------------|------------|------------|------------|-----------|-------------|------------|-----------|-------------|-------------|-----------|-----------|------------|------------|------------|-------------|------------|-------------|-------------|------------|------------|------------|------------|------------|-----------|-----------|------------|-----------|------------|------------|-------------|------------|------------|-------------|-------------|------------|----------|------------|-------------|------------|------------|-------------|-------------|------------|----------|------------|-------------|------------|------------|------------|------------|-------------|------------|-------------|------------|------------|------------|-------------|-------------|-------------|-------------|-------------|-----------|-------------|-------------|-------------|-------------|-------------|
|Gb1-2      |Gb1-2 GABA    |Gb1-2 GABA       |None      |0.0012162228|0.00015655666|0.0007011017     |0.008910024       |2.604616e-05      |0.0013842874|0.0011577017|0.06381319|0.0003416799|0.00043809053|0.0017023097      |2.639922e-05|0.00017739477|0.00029747118|0.0020205898|0.002436088|0.010450783|1.6771797e-06|0.001953325|0.0019771748|0.005507645|0.0042126277|0.0035268543|0.0011763936|0.0014934228|0.0014862448|0.0016202381|0.00033783494|0.0128647005|0.00011733183|0.00020410866|0.0015136964|0.00036336062|0.008227488|0.008641888|0.0012168279|0.024405263|0.0006090624|0.01752709|0.0014401844|0.017013999|0.0012381236|0.0027838773|0.0035406314|4.5887286e-05|0.00012445106|0.0071767657|0.0019897874|0.0097302785|0.0004271934|0.0048034517|0.0001381562|0.00022291896|0.00034220653|0.04667171|0.00322487|0.029999454|8.977511e-05|0.002004063|0.0006713049|0.001634424|0.0012884187|0.0066387546|0.0015726861|0.00028648277|0.028014915|0.0004048754|0.0003920776|0.00019321572|0.00079886295|0.009393225|0.0009428808|0.027355352|0.009530359|0.010141343|0.0005618859|0.00053159415|0.0037294286|0.00044019768|0.012636062|0.011418089|0.00017641549|0.000799447|0.00032995237|0.012780122|0.00081895955|0.0003811782|0.0064510047|0.0032412473|0.00028182066|0.004049578|2.0678754e-06|0.0047364566|0.2594797|0.0049051535|0.0008632637|0.0035660195|0.0054662316|0.0004969132|0.0054825847|0.0005060083|4.7758833e-05|0.0064058984|0.0008465315|0.0025367003|0.067800835|0.00014405268|1.139286e-05|0.008258109|0.00046846634|0.00083107647|0.024232045|0.012851092|0.0018063982|0.0057125525|4.396907e-06|0.00028940683|0.0007546474|0.00027128105|0.00058339076|0.0029196339|0.0004083942|0.0053130514|6.668573e-05|0.0001863813|0.009039877|0.001448437|0.0026516842|0.017183386|8.395138e-05|0.0025180348|0.00022454833|0.0012995326|0.0004234898|0.00069253385|0.00074663066|0.0013974423|0.03398079|0.0013081075|0.00044675352|0.0021949154|0.0004709228|0.00014847725|0.00040034397|0.0005949922|0.00135979|0.0009407896|7.7904166e-05|7.695236e-05|5.366119e-05|0.0003188306|0.0004784335|0.00046772635|0.0001581345|0.00092866743|0.0023081445|0.0011137992|0.0001318106|6.9092326e-05|1.4297381e-05|0.00013958917|4.3109252e-05|0.00039777902|3.44212e-05|5.3402524e-05|4.3644093e-05|4.4868156e-07|4.3928667e-05|0.00010212503|


### Access

```{caution}
Public data access to `hm3d_regionwise_relative_mass_wo_fibers.csv` is currently unavailable as HypoMap3D is not yet published.
```
