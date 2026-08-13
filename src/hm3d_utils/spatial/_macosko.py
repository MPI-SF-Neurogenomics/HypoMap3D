import typing
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import skimage.io as ski_io
import os
from pathlib import Path
import nrrd
import re

def restrict_to_hy_boundary(adata, bbox_ccf):
    x_max, y_max, z_max = bbox_ccf.max 
    x_min, y_min, z_min = bbox_ccf.min

    keep_mask = (
        adata.obs.x_ccf <= x_max) & (
        adata.obs.x_ccf >= x_min) & (
        adata.obs.y_ccf <= y_max) & (
        adata.obs.y_ccf >= y_min) & (
        adata.obs.z_ccf <= z_max) & (
        adata.obs.z_ccf >= z_min
        )

    return adata[keep_mask, :]
 

def map_10um_atlas_annotations(image_array: np.ndarray, ccf_annotation_file_path: str) -> np.ndarray:
    """
    Maps a 3D image array to a 2D atlas annotation using a provided CCF annotation file.

    This function takes a 3D `image_array` (which represents the registration of a TIFF file) and maps 
    its voxel coordinates to the corresponding CCF (Common Coordinate Framework) annotation labels 
    from the specified `ccf_annotation_file_path`. The result is a 2D array where each element corresponds 
    to the CCF annotation at that position.

    Parameters:
    ----------
    image_array : np.ndarray
        A 3D numpy array representing the image data from a registered TIFF file.
        The dimensions of the array represent the x, y, and z coordinates in the brain's space.

    ccf_annotation_file_path : str
        Path to the CCF annotation file (typically in NRRD format). This file contains the anatomical labels 
        for the brain regions in the same space as `image_array`.

    Returns:
    -------
    np.ndarray
        A 2D numpy array where each element represents the CCF annotation label 
        for the corresponding pixel in the 2D projection of `image_array`. 
        If the voxel is out of bounds in the CCF annotation data, the value will be 0.

    Notes:
    -----
    - The function expects that the `image_array` is a registered 3D array, and it maps the coordinates 
      from the 3D space to a 2D annotation space.
    - The values from the `image_array` are scaled by a factor of 100 and converted to integers to match 
      the corresponding CCF annotation data's resolution.
    - If any coordinate is out of bounds of the CCF annotation array, the corresponding annotation is set to 0.
    """
    
    # Round and scale the image array values for mapping
    imarray_10um = np.round(image_array * 100).astype(int)

    # Read the CCF annotation data from the NRRD file
    ccf_anno, header = nrrd.read(ccf_annotation_file_path)

    # Initialize an empty matrix to store the 2D annotation result
    out = np.empty((imarray_10um.shape[0], imarray_10um.shape[1]))

    # Map the image array coordinates to the CCF annotation labels
    for i in range(imarray_10um.shape[0]):
        for j in range(imarray_10um.shape[1]):
            x = imarray_10um[i, j, 0]
            y = imarray_10um[i, j, 1]
            z = imarray_10um[i, j, 2]

            # Check if the coordinates are out of bounds in the CCF annotation
            if (x >= ccf_anno.shape[0]) or (y >= ccf_anno.shape[1]) or (z >= ccf_anno.shape[2]):
                out[i, j] = 0  # Set to 0 if out of bounds
            else:
                out[i, j] = ccf_anno[x, y, z]  # Map annotation value

    # Convert the output array to integers
    out = out.astype(int)
    
    return out


def load_pucks(puck_dir: str, registration_dir: str, only_load: list = []) -> ad.AnnData:
    """
    Loads and concatenates multiple puck datasets into a single AnnData object.

    This function loads puck data stored in `.h5ad` files from the specified `puck_dir`. By default, it loads 
    all pucks for which registration data is available in the `registration_dir`. However, you can override 
    this behavior by passing a list of specific puck numbers (integers) to the `only_load` parameter. 
    Each puck dataset is loaded, its metadata is adjusted, and the data is concatenated into one AnnData object.

    Parameters:
    ----------
    puck_dir : str
        The directory containing the `.h5ad` files for each puck. Each file should follow the naming convention
        `Puck_Num_<num>.h5ad`, where `<num>` is the puck number.
    
    registration_dir : str
        The directory containing the registration files for the pucks. The registration files should have 
        filenames that include the puck numbers.

    only_load : list, optional
        A list of integers representing specific puck numbers to load. If empty, the function will load all pucks 
        for which registration files are available in the `registration_dir`. Default is an empty list.

    Returns:
    -------
    ad.AnnData
        An AnnData object containing the concatenated data from all the loaded pucks. The AnnData object includes:
        - `obs`: The observation (rows) data from each puck, with an added "puck_num" column to denote which puck 
          each observation belongs to.
        - `var`: The feature (columns) data, with the variable names mapped based on the "features" in the `.h5ad` file.

    Notes:
    -----
    - The function assumes that the `.h5ad` files are named using the convention `Puck_Num_<num>.h5ad`, where `<num>`
      is a one- or two-digit integer (e.g., `Puck_Num_01.h5ad` or `Puck_Num_10.h5ad`).
    - If no puck numbers are provided in `only_load`, it will automatically detect all pucks available in the 
      `registration_dir` by extracting numbers from the filenames.
    - The puck number is added to the `obs` dataframe of each loaded AnnData object under the column "puck_num".
    - The variable names (`var_names`) are updated using the corresponding "features" from the `var` attribute of each 
      AnnData object.
    - The function concatenates all the pucks using the outer join on the index.
    """
    
    # get puck numbers to load
    if len(only_load) == 0:
        puck_nums = [int(re.findall("^\d{1,3}", string)[0]) for string in sorted(os.listdir(registration_dir))]
    else:
        puck_nums = only_load

    # load and append each puck dataset
    puck_ans = []
    for puck_num in puck_nums:
        print(f"Loading puck {puck_num}")

        # Format puck number with leading zero if necessary
        if puck_num < 10:
            puck_num = "0" + str(puck_num)

        # Load the AnnData object for the puck
        puck_an = ad.read_h5ad(puck_dir + f"/Puck_Num_{puck_num}.h5ad")
        print(f"n_obs x n_var = {len(puck_an.obs_names)} x {len(puck_an.var_names)}")

        # Add puck number to the obs dataframe
        puck_an.obs["puck_num"] = puck_num

        # Update variable names using the 'features' column in the var dataframe
        puck_an.var_names = puck_an.var_names.map(puck_an.var["features"])

        # Append the puck AnnData object to the list
        puck_ans.append(puck_an)

    # Concatenate all puck datasets
    print(f"Concatenating...")
    ss_an = ad.concat(puck_ans, join="outer")  # concatenate along the rows (observations)

    return ss_an


def curate_puck_an_obs(puck_an: ad.AnnData) -> ad.AnnData:
    """
    Curates the `.obs` DataFrame of an AnnData object by renaming columns and dropping unnecessary ones.

    This function modifies the `.obs` attribute of the provided AnnData object by:
    - Renaming specific columns to a more consistent or meaningful format.
    - Dropping columns that are deemed unnecessary for further analysis.

    The changes are made in place and the modified AnnData object is returned.

    Parameters:
    ----------
    puck_an : ad.AnnData
        An AnnData object containing spatial transcriptomics data, with the `.obs` attribute holding metadata 
        for the observations (e.g., coordinates and sample information).

    Returns:
    -------
    ad.AnnData
        The same AnnData object (`puck_an`), but with its `.obs` DataFrame curated by renaming and dropping 
        columns as described above.

    Notes:
    -----
    - The following columns are renamed:
        - "NisslTiffX" -> "x_nissl"
        - "NisslTiffY" -> "y_nissl"
    
    - The following columns are dropped:
        - "TopStruct", "DeepCCF", "CCF_acronym", "CCF_ID", "CCF_Name", 
        - "Raw_Slideseq_X", "Raw_Slideseq_Y", "nCount_Spatial.1", "nFeature_Spatial.1",
        - "CCF_LeftRight", "IsOutsideCCF", "Correponding_DissectateSegName", 
        - "Correponding_DissectateSegID", "BeadBarcode", "CCF_X", "CCF_Y", "CCF_Z"
    
    - This function is designed to streamline the `.obs` DataFrame, removing redundant or extraneous columns
      that may not be relevant for subsequent analyses or visualization.
    """
    
    obs = puck_an.obs

    # Rename specific columns for consistency or clarity
    obs.rename(columns={
        "NisslTiffX": "x_nissl",
        "NisslTiffY": "y_nissl",
    }, inplace=True)

    # Drop unnecessary columns
    obs.drop(columns=[
        'TopStruct', 
        'DeepCCF', 
        'CCF_acronym', 
        'CCF_ID', 
        'CCF_Name',
        'Raw_Slideseq_X', 
        'Raw_Slideseq_Y',
        'nCount_Spatial.1',
        'nFeature_Spatial.1',
        'CCF_LeftRight', 
        'IsOutsideCCF', 
        'Correponding_DissectateSegName', 
        'Correponding_DissectateSegID',
        'BeadBarcode',
        "CCF_X",
        "CCF_Y",
        "CCF_Z",
    ], inplace=True)

    return puck_an


def map_new_ccf_coords_to_macosko(ss_an: ad.AnnData, ccf_annotation_file: str, registration_dir: str) -> ad.AnnData:
    """
    Maps updated CCF coordinates and atlas annotations to the Macosko dataset, and filters out cells outside the brain.

    This function updates the spatial coordinates (`x_ccf`, `y_ccf`, `z_ccf`) and the corresponding Allen atlas annotations
    (`allen_annot`) for each cell in the given AnnData object (`ss_an`). The mapping is performed using new CCF registration 
    data (stored in `.png_Coords.tif` files) and an Allen CCF annotation file.

    The function also filters out cells that are outside the brain, as indicated by an annotation value of 0.

    Parameters:
    ----------
    ss_an : ad.AnnData
        The AnnData object containing the Macosko dataset. This object should have spatial metadata in `.obs`, 
        including the coordinates (`x_nissl`, `y_nissl`) of each cell. The updated coordinates and annotations 
        will be added to the `.obs` dataframe of this object.

    ccf_annotation_file : str
        Path to the CCF annotation file (e.g., in NRRD format). This file contains the brain atlas annotations 
        to be mapped to the cells.

    registration_dir : str
        Directory containing the updated CCF registration files (e.g., `.png_Coords.tif` files). These files 
        are used to map the cell coordinates to the CCF annotation space.

    Returns:
    -------
    ad.AnnData
        A new AnnData object, identical to `ss_an`, but with the following updates in its `.obs` dataframe:
        - `x_ccf`, `y_ccf`, `z_ccf`: Updated CCF coordinates for each cell.
        - `allen_annot`: The corresponding Allen CCF annotation for each cell.
        - `-y_ccf`, `-z_ccf`: Negative values of `y_ccf` and `z_ccf` for convenient plotting (flipping the coordinates).
        
        Additionally, cells outside the brain (where `allen_annot == 0`) are filtered out.

    Notes:
    -----
    - The function expects each puck to have a corresponding `.png_Coords.tif` file in the `registration_dir`. 
      If such a file is not found for a particular puck, no updates are made for that puck.
    - The function uses the `map_10um_atlas_annotations` function to map the updated coordinates to CCF annotations.
    - The `allen_annot` column is populated with the annotation from the Allen CCF that corresponds to the mapped coordinates.
    - Cells that do not map to any valid brain region (i.e., `allen_annot == 0`) are excluded from the final result.

    Example:
    --------
    updated_ss_an = map_new_ccf_coords_to_macosko(ss_an, "path/to/ccf_annotation.nrrd", "path/to/registration")
    """
    
    # Make a copy of the input AnnData object
    an = ss_an.copy()

    # Load the CCF annotation file
    print(f"Loading allen atlas annotation file...")
    ccf_anno, header = nrrd.read(ccf_annotation_file)

    # Initialize new columns in .obs for the updated coordinates and annotations
    an.obs["x_ccf"] = np.nan
    an.obs["y_ccf"] = np.nan
    an.obs["z_ccf"] = np.nan
    an.obs["allen_annot"] = np.nan
        
    for puck_num in an.obs.puck_num.unique():
        print(f"Checking data availability for puck {puck_num}")

        # Check if new registration is available
        registration_path = Path(registration_dir + f"/{puck_num}.png_Coords.tif")
        if registration_path.exists():
            print("New registration available! Updating...")

            # Load the updated CCF coordinates
            new_coords = ski_io.imread(registration_path)
            imarray = np.array(new_coords)

            # Map the coordinates to Allen CCF annotations
            mapped = map_10um_atlas_annotations(imarray, ccf_annotation_file,)

            # Update the coordinates and annotations for the current puck
            rows = an.obs.puck_num == puck_num
            an.obs.loc[rows, "x_ccf"] = imarray[an.obs.loc[rows, "y_nissl"], an.obs.loc[rows, "x_nissl"], 0]
            an.obs.loc[rows, "y_ccf"] = imarray[an.obs.loc[rows, "y_nissl"], an.obs.loc[rows, "x_nissl"], 1]
            an.obs.loc[rows, "z_ccf"] = imarray[an.obs.loc[rows, "y_nissl"], an.obs.loc[rows, "x_nissl"], 2]
            an.obs.loc[rows, "allen_annot"] = mapped[an.obs.loc[rows, "y_nissl"], an.obs.loc[rows, "x_nissl"]]

        else:
            print("No new registration available. No data updated.")

    # Add negative coordinates for convenient plotting
    an.obs["-y_ccf"] = -an.obs.loc[:, "y_ccf"]
    an.obs["-z_ccf"] = -an.obs.loc[:, "z_ccf"]

    # Filter cells that are outside the brain (allen_annot == 0)
    an = an[an.obs.allen_annot != 0]
            
    return an

def map_allen_annot(ss_an: ad.AnnData, query_file_path: str) -> ad.AnnData:
    """
    Maps Allen annotations to the cells in an AnnData object based on provided annotation data.

    This function takes an AnnData object (`ss_an`) and a query file (CSV format) containing Allen brain region annotations. 
    It maps the Allen region annotations from the query file to the corresponding cells in the AnnData object using the 
    `allen_annot` column of the AnnData's `.obs` dataframe. The mapped annotations are stored in a new column called 
    `allen_names` within the `.obs` dataframe of the AnnData object.

    Parameters:
    ----------
    ss_an : ad.AnnData
        The AnnData object containing the spatial transcriptomics data, with an `allen_annot` column in `.obs` 
        representing the Allen brain region identifiers for each cell.

    query_file_path : str
        Path to a CSV file containing the Allen brain region annotations. The file should have two columns:
        - "id": The identifiers of the Allen regions (matching the values in `allen_annot`).
        - "name": The names of the corresponding Allen regions.

    Returns:
    -------
    ad.AnnData
        The input AnnData object (`ss_an`), with an additional column `allen_names` in the `.obs` dataframe containing 
        the mapped Allen region names for each cell.

    Notes:
    -----
    - This function assumes that the `allen_annot` column in `.obs` of the AnnData object contains identifiers that 
      can be mapped to the `id` column in the query CSV.
    - The `query_file_path` CSV must have an "id" column with identifiers that match the values in `allen_annot` 
      and a "name" column with the corresponding region names.
    - The function modifies the AnnData object in place, adding the `allen_names` column to `.obs`.

    Example:
    --------
    updated_ss_an = map_allen_annot(ss_an, "path/to/query_file.csv")
    """
    
    # Load the query file containing Allen brain region annotations
    query = pd.read_csv(query_file_path)
    query = query.set_index("id")

    # Map the Allen region names to the AnnData object based on the allen_annot values
    ss_an.obs["allen_names"] = ss_an.obs.allen_annot.map(query["name"])

    return ss_an

def create_full_macosko(puck_dir: str, ccf_annotation_file: str, registration_dir: str,
                        query_file_path: str, only_load: list = [],
                        save_path: str|None = None,
                       save: bool = False,
                       bbox_ccf=None) -> ad.AnnData:
    """
    Creates a full Slide-seq object for the Macosko dataset with updated CCF coordinates, Allen annotations.
    Optionally restricts to the hypothalamus boundary, and saves the resulting AnnData object to a file.

    This wrapper function loads the Slide-seq puck data, curates the `.obs` dataframe, maps updated CCF coordinates 
    based on the provided CCF annotation file, adds Allen brain region annotations, and filters the dataset to 
    include only cells within the hypothalamus region. If the `save` parameter is set to `True`, the resulting 
    AnnData object will be saved to the specified `save_path`. If a BoundingBox3D object is passed to bbox_ccf, 
    the AnnData object will be filtered to that bounding box. The function returns an AnnData object representing the 
    full Macosko dataset with the necessary metadata and filtering.

    Parameters:
    ----------
    puck_dir : str
        Directory containing the puck data files. Each puck corresponds to a Slide-seq dataset.

    ccf_annotation_file : str
        Path to the CCF annotation file (e.g., in NRRD format). This file contains the brain atlas annotations 
        to be used for mapping CCF coordinates.

    registration_dir : str
        Directory containing the registration files (e.g., `.png_Coords.tif` files). These files are used to 
        map updated spatial coordinates to the CCF space.

    query_file_path : str
        Path to a CSV file containing Allen brain region annotations. This file should have two columns:
        - "id": The identifiers of the Allen regions (matching the values in `allen_annot`).
        - "name": The names of the corresponding Allen regions.

    only_load : list of int, optional, default=[]
        A list of specific puck indices to load. If empty, the function loads all pucks for which registration 
        files are available. If provided, only the specified puck indices will be loaded.

    save_path : str, optional, default=None
        The path to save the resulting AnnData object as an H5AD file. This parameter is required if the `save` 
        parameter is set to `True`. If not provided or `save` is `False`, the data will not be saved.

    save : bool, optional, default=False
        Whether to save the resulting AnnData object to the `save_path`. If `True`, the object will be saved as 
        an H5AD file to the specified location.

    bbox_ccf: None | BoundingBox3D object.
        If a BoundingBox3D object is passed, the data will be filtered using the restricht_to_hy_boundary function
        using the coordinates from the BoundingBox3D.

    Returns:
    -------
    ad.AnnData
        An AnnData object containing the Macosko dataset with the following updates:
        - Curated `.obs` DataFrame with cleaned metadata.
        - Updated CCF coordinates (`x_ccf`, `y_ccf`, `z_ccf`) and the corresponding Allen region annotations (`allen_names`).
        - Cells outside the brain (with `allen_annot == 0`) are filtered out.
        - Cells are further restricted to the hypothalamus boundary based on the `restrict_to_hy_boundary` function.

    Notes:
    -----
    - This function depends on the successful loading of the puck data, ccf_annotation_file, registration files, 
      and the Allen query file. If any of these files are missing or incorrect, the function may not work as expected.
    - The resulting AnnData object will have the `allen_names` column in `.obs` representing the mapped Allen brain 
      regions for each cell, based on the `allen_annot` values.
    - The `restrict_to_hy_boundary` function is applied to restrict the dataset to cells inside the hypothalamus region.
    - The `save` parameter must be `True` and `save_path` must be provided for the data to be saved.

    Example:
    --------
    macosko_anndata = create_full_macosko(puck_dir="path/to/pucks", 
                                          ccf_annotation_file="path/to/ccf_annotation.nrrd", 
                                          registration_dir="path/to/registrations", 
                                          query_file_path="path/to/query_file.csv",
                                          save=True, save_path="path/to/save_output.h5ad")
    """
    
    # Load and curate the Macosko data
    macosko = load_pucks(puck_dir=puck_dir, 
                       registration_dir=registration_dir,
                       only_load=only_load)

    macosko = curate_puck_an_obs(macosko)
    
    # Map updated CCF coordinates and Allen annotations
    macosko = map_new_ccf_coords_to_macosko(ss_an=macosko, 
                                            ccf_annotation_file=ccf_annotation_file, 
                                            registration_dir=registration_dir)

    macosko = map_allen_annot(ss_an=macosko, query_file_path=query_file_path)

    # Restrict the dataset to cells inside the hypothalamus boundary
    if bbox_ccf:
        print("restricting data to hypothalamus boundary")
        macosko = restrict_to_hy_boundary(macosko, bbox_ccf)

    if save==True:
        macosko.write_h5ad(save_path)

    return macosko

