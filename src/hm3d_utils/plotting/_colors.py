from matplotlib.colors import LinearSegmentedColormap

cmap_marker = LinearSegmentedColormap.from_list(
    "marker_red", 
    [
        (0.0, "lightgrey"),  # value = 0
        (0.5, "crimson"),        # value = 0.5
        (1.0, "darkred")     # value = 1
    ]
)