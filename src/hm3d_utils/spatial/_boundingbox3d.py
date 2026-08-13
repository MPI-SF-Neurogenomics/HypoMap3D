import json
import numpy as np

# ============================================================================ #
#                               Exposed Class(es)                              #
# ============================================================================ #

class BoundingBox3D:
    """
    A class representing a 3D axis-aligned bounding box using NumPy arrays.

    Attributes:
        min (np.ndarray): Minimum bounds along x, y, z axes (shape: (3,)).
        max (np.ndarray): Maximum bounds along x, y, z axes (shape: (3,)).
        unit (str): Unit of measurement for the bounding box coordinates (e.g., 'mm').
    """

    def __init__(self,
                 min_coord: list[float] | tuple[float, float, float] | np.ndarray,
                 max_coord: list[float] | tuple[float, float, float] | np.ndarray,
                 unit: str = "unitless") -> None:
        """
        Initialize a 3D bounding box with min and max coordinates.

        Args:
            min_coord: Minimum coordinates (x, y, z).
            max_coord: Maximum coordinates (x, y, z).
            unit: Unit of measurement (default: 'unitless').
        """
        self.min = np.array(min_coord, dtype=float)
        self.max = np.array(max_coord, dtype=float)
        self.unit = unit

        if self.min.shape != (3,) or self.max.shape != (3,):
            raise ValueError("min_coord and max_coord must be 3-element vectors.")
        if np.any(self.min > self.max):
            raise ValueError("Each element of min_coord must be <= corresponding element of max_coord.")

    def save(self, filepath: str) -> None:
        """
        Save the bounding box coordinates and metadata to a JSON file.

        Args:
            filepath: Path to the JSON file.
        """
        data = {
            "min": self.min.tolist(),
            "max": self.max.tolist(),
            "unit": self.unit
        }
        with open(filepath, "w") as f:
            json.dump(data, f, indent=4)

    @classmethod
    def load(cls, filepath: str) -> "BoundingBox3D":
        """
        Load bounding box coordinates and metadata from a JSON file.

        Args:
            filepath: Path to the JSON file.

        Returns:
            A new BoundingBox3D instance with loaded coordinates and unit.
        """
        with open(filepath, "r") as f:
            data = json.load(f)
        return cls(min_coord=data["min"], max_coord=data["max"], unit=data.get("unit", "unitless"))

    def isinbounds(self, coordinates: np.ndarray) -> np.ndarray:
        """
        Check if points are inside the bounding box.

        Args:
            coordinates: Array of shape (N, 3) representing points.

        Returns:
            Boolean array of length N indicating points inside the box.
        """
        coordinates = np.asarray(coordinates, dtype=float)
        if coordinates.shape[1] != 3:
            raise ValueError("Input coordinates must have shape (N, 3).")

        return np.all((coordinates >= self.min) & (coordinates <= self.max), axis=1)

    def __repr__(self) -> str:
        """
        Return a human-readable string representation of the bounding box.
        """
        return (
            f"Bounding Box Ranges ({self.unit}):\n"
            f" x: {self.min[0]}<->{self.max[0]}\n"
            f" y: {self.min[1]}<->{self.max[1]}\n"
            f" z: {self.min[2]}<->{self.max[2]}"
        )
