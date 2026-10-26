# Taxonomy Objects

The HypoMap3D cluster taxonomy spans five levels of varying granularity which can be summarised as follows:

```{mermaid}
graph TD
    A["Level 1 | Class"] --> B["Level 2 | Subclass"]
    B --> C["Level 3 | Supertype"]
    C --> D["Level 4 | Cell Type"]
    D --> E["Level 5 | Subtype"]
```

**Level 4** presents the highest level of granularity under our defined constraints of a cell type classification. **Level 5** provides further partitions with biological evidence but insufficient evidence for cell type classification. This could indicate the capture of finer cell states within each cell type. These fine partitions are included within our taxonomy to provide a more detailed view of cellular heterogeneity.

## Cluster Taxonomic Summary

A simple summary JSON file provides an overview of each cluster within the HypoMap3D taxonomy to include:

* Parent
* Children
* Sisters
* Path
* Alternative annotations (descriptive and allen brain-like)

### Example Entry
```json
"Gl2-35": {
    "descr": "Gl2-35 VMH Nr5a1.Acvr1c",
    "abc_like": "Gl2-35 VMH Nr5a1 Acvr1c Glut",
    "parent": "Gl1-1",
    "path": ["Gl1-1", "Gl2-35"],
    "children": ["Gl3-113", "Gl3-111", "Gl3-114", "Gl3-115", "Gl3-112"],
    "sisters": [],
    "n_sisters": 0,
    "n_children": 5
}
```

### Access

```{caution}
Public data access to `hm3d_clusthm3d_cluster_taxonomy.json` is currently unavailable as HypoMap3D is not yet published.
```