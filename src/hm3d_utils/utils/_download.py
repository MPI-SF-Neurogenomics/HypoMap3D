from pathlib import Path
import requests
from tqdm.auto import tqdm


def download_file(url: str, output_path: str | Path, overwrite: bool = False):
    """
    Download a file from a URL.

    Parameters
    ----------
    url
        Direct download URL.
    output_path
        Where to save the file.
    overwrite
        Whether to overwrite an existing file.
    """
    output_path = Path(output_path)

    if output_path.exists() and not overwrite:
        print(f"{output_path.name} already exists.")
        return output_path

    output_path.parent.mkdir(parents=True, exist_ok=True)

    response = requests.get(url, stream=True)
    response.raise_for_status()

    total = int(response.headers.get("content-length", 0))

    with open(output_path, "wb") as f, tqdm(
        total=total,
        unit="B",
        unit_scale=True,
        desc=output_path.name,
    ) as pbar:
        for chunk in response.iter_content(chunk_size=1024 * 1024):
            if chunk:
                f.write(chunk)
                pbar.update(len(chunk))

    print(f"Saved to {output_path}")
    return output_path