# type: ignore

from dataclasses import dataclass, field
import numpy as np
import numpy.typing as npt
from typing import List, NoReturn


@dataclass
class Beam:
    fname: str
    data: npt.ArrayLike = field(init=False)
    phi_col: int = field(init=False)
    theta_col: int = field(init=False)
    volt_cols: List[int] = field(init=False)

    def __post_init__(self):

        f = open(self.fname, "r")
        header = f.readline().split(",")
        f.close()
        column_names = []
        units = []
        for col in header:
            column_name, unit = tuple(col.split("["))
            column_names.append("".join(column_name.lower().split(" ")))
            units.append(unit.lower().strip())

        self.data = np.loadtxt(self.fname, skiprows=1, delimiter=",")
        self.theta_col = np.where(np.array(column_names) == "theta")[0][0]
        self.phi_col = np.where(np.array(column_names) == "phi")[0][0]
        self.volt_cols = [2, 3, 4, 5, 6, 7]  # should be the last cols
        # check that volt col is the last col
        assert self.phi_col not in self.volt_cols  # not the same as phi
        assert self.theta_col not in self.volt_cols  # not the same as theta
        assert len(self.data[0]) == len(self.volt_cols) + 2

        # convert to radians
        for col in [self.theta_col, self.phi_col]:
            if "deg" in units[col]:
                self.data[:, col] = np.radians(self.data[:, col])
        for col in self.volt_cols:
            if "mv" in units[col]:
                self.data[:, col] /= 1e3  # convert mV to V

    def _delete_wrap(
        self, col: int = 0, period: float = 2 * np.pi
    ) -> np.ndarray:
        """
        In case one of the axis is sampled around the circle (eg includes both
        0 and 360), we delete the values for the second round around the circle
        """
        wrap_idx = np.argwhere(
            self.data[:, col] - self.data[:, col].min() >= period
        )[
            :, 0
        ]  # rows where the angle wraps around
        if len(wrap_idx):  # axis wraps around
            new_data = np.delete(
                self.data, wrap_idx, axis=0
            )  # delete the given rows
        assert new_data[:, col].max() - new_data[:, col].min() < period
        return new_data

    def _rev_indices(self) -> np.ndarray:
        """
        Convert from phi changing fastest to theta changing fastest
        (Fortran style to C style)
        """
        ph = self.data[:, self.phi_col]
        th = self.data[:, self.theta_col]
        new_data = np.empty_like(self.data)
        for i in range(new_data.shape[1]):  # loop over columns
            dcol = self.data[:, i].copy()
            new_data[:, i] = dcol.reshape(
                np.unique(ph).size, np.unique(th).size, order="F"
            ).flatten(order="C")

        return new_data

    def _shift_phi(self) -> NoReturn:
        """
        Change phi from going -pi->pi to 0->2pi
        """
        phi = self.data[:, self.phi_col]
        neg_phi = np.argwhere(phi < 0)[:, 0]
        # assert np.allclose(np.diff(neg_phi), np.ones_like(neg_phi)[:-1])
        new_data = self.data.copy()
        new_data[neg_phi, self.phi_col] += 2 * np.pi
        neg_phi_max = len(neg_phi)
        new_data = np.concatenate(
            (new_data[neg_phi_max:], new_data[:neg_phi_max])
        )
        assert (
            np.diff(new_data[:, self.phi_col]).min() >= 0
        )  # phi is still increasing
        assert (
            new_data[:, self.phi_col].min() >= 0
        )  # and is bigger than 0 everywhere
        return new_data

    def process_data(self) -> NoReturn:
        self.data = self._delete_wrap(col=self.theta_col, period=np.pi)
        self.data = self._delete_wrap(col=self.phi_col, period=2 * np.pi)
        self.data = self._rev_indices()
        self.data = self._shift_phi()

    def efield_to_power(self) -> np.ndarray:
        power = np.empty((len(self.data), 3))
        power[:, 0] = self.data[:, self.theta_col]
        power[:, 1] = self.data[:, self.phi_col]
        power[:, 2] = np.sqrt(
            np.sum(self.data[:, self.volt_cols] ** 2, axis=1)
        )

        return power

    def to_txt(self, outpath: str, to_power=True) -> NoReturn:
        if to_power:
            power = self.efield_to_power()
            np.savetxt(
                outpath,
                power,
                header="Theta [rad] Phi [rad] Abs(V) [V] \n\n",
                comments="",
            )
        else:
            raise NotImplementedError()
            # figure out what the format must be


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("in_fname", metavar="-if", type=str)
    parser.add_argument("out_fname", metavar="-of", type=str)
    args = parser.parse_args()
    beam = Beam(args.in_fname)
    beam.process_data()
    beam.to_txt(args.out_fname)
