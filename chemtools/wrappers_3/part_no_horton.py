import denspart


__all__ = ['DensPart']


class DensPart(object):
    """Density-based Atoms-in-Molecules Partitioning Class."""

    available_schemes = ["mbis"]

    def __init__(self, coordinates, numbers, density, grid, scheme="mbis", **kwargs):
        if scheme.lower() not in DensPart.available_schemes:
            raise NotImplementedError("MBIS is currently the only supported scheme.")
        self.coordinates = coordinates
        self.numbers = numbers
        self.density = density
        self.grid = grid
        self.part = denspart.mbis.partition(self.numbers,
                                            self.coordinates,
                                            self.grid,
                                            self.density
                                            **kwargs)

