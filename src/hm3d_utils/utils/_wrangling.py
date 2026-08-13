import pandas as pd 
import warnings 
import concurrent
import numpy as np


def restrict_dbscanmetrics_to_top_k(df: pd.DataFrame, metric: str, k: int = 3) -> pd.DataFrame:
    """
    function to restrict the desired metric of the dbscan metrics to the top k dbscan labels.
    all other labels will be set to the "rest" label for which the metric will be summed.

    returns df 
    """

    cluster_names_list = []
    otsu_thr_list = []
    label_list = []
    num_vox_list = []
    metric_list = []
    rank_list = []

    for cluster in df.cluster_name.unique():

        cluster_data = df.loc[df.cluster_name == cluster, :].sort_values(by=metric, ascending = False)

        # keep the top k
        top_k = cluster_data.iloc[:k, :]
        # create a column for the rank
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            top_k["size_rank"] = [str(i) for i in range(len(top_k))]
        
        cluster_names_list.extend(top_k.cluster_name.to_list())
        otsu_thr_list.extend(top_k.otsu_thr.to_list())
        label_list.extend(top_k.dbscan_label.to_list())
        num_vox_list.extend(top_k.num_vox.to_list())
        metric_list.extend(top_k[metric].to_list())
        rank_list.extend(top_k.size_rank.to_list())

        # calculate rests
        rest = cluster_data.iloc[k:, :]
        if len(rest) > 0:
            cluster_names_list.append(cluster)
            otsu_thr_list.append(rest.otsu_thr.iloc[0])
            label_list.append("rest")
            num_vox_list.append(rest.num_vox.sum())
            metric_list.append(rest[metric].sum())
            rank_list.append("rest")

    out = pd.DataFrame.from_dict({
        
        "cluster_name": cluster_names_list,
        "otsu_thr": otsu_thr_list,
        "dbscan_label": label_list,
        "num_vox": num_vox_list,
        metric: metric_list,
        "size_rank": rank_list
    })
    
    return out


def wrangle_dbscan(df):

    out=df.copy()

    # num clusters
    nums = out.groupby("cluster_name")["dbscan_label"].count().rename("num_postsmooth_dbscan_clusters")

    # top 2 and "rest" clusters 
    top_k = obtain_wide_top_k_df(out)

    
    return pd.concat([nums, top_k], axis=1)


def obtain_wide_top_k_df(df):
    
    top_k = restrict_dbscanmetrics_to_top_k(df, metric="percent_vox", k = 2)

    clusters = top_k.cluster_name.unique()

    # Use ThreadPoolExecutor to process clusters concurrently
    all_metrics = []
    with concurrent.futures.ThreadPoolExecutor() as executor:
        # Start the load opeartions by submitting tasks to the thread pool and collect the results
        future_to_cluster = {
            executor.submit(get_wide_top_k_for_cluster, top_k, cluster): cluster
            for cluster in clusters
        }

        # Retrieve the results as they complete
        for future in concurrent.futures.as_completed(future_to_cluster):
            cluster_metrics = future.result()
            all_metrics.append(cluster_metrics)

    # Create a DataFrame from the list of dictionaries
    df = pd.DataFrame(all_metrics)

    # extract values from lists
    for col in ["num_vox_biggest_dbscan_cluster", 
                "percent_vox_biggest_dbscan_cluster",
                "num_vox_second_biggest_dbscan_cluster",
                "percent_vox_second_biggest_dbscan_cluster",
                "num_vox_rest_dbscan_cluster",
                "percent_vox_rest_dbscan_cluster"]:

        df[col] = df[col].apply(lambda x: x[0] if len(x) > 0 else np.nan)
        
    return df.set_index("cluster_label")


def get_wide_top_k_for_cluster(df, cluster_name):
    data = df.loc[df.cluster_name == cluster_name,:]

    res = {
        "cluster_label": cluster_name,
        "num_vox_biggest_dbscan_cluster": data.loc[data.size_rank == "0", "num_vox"].values,
        "percent_vox_biggest_dbscan_cluster": data.loc[data.size_rank == "0", "percent_vox"].values,
        "num_vox_second_biggest_dbscan_cluster": data.loc[data.size_rank == "1", "num_vox"].values,
        "percent_vox_second_biggest_dbscan_cluster": data.loc[data.size_rank == "1", "percent_vox"].values,
        "num_vox_rest_dbscan_cluster": data.loc[data.size_rank == "rest", "num_vox"].values,
        "percent_vox_rest_dbscan_cluster": data.loc[data.size_rank == "rest", "percent_vox"].values,
        
    }
    return res
