from functools import reduce
import numpy as np


class Pauli:
    """
    """

    X_ = np.array([[0, 1], [1, 0]], dtype=np.complex128)
    Z_ = np.array([[1, 0], [0, -1]], dtype=np.complex128)

    def __init__(self, z, x):
        """
        """
        self.mat_ = np.eye(2, dtype=np.complex128)
        self.z_ = z % 2
        self.x_ = x % 2
        if self.z_ == 1:
            self.mat_ = self.mat_.dot(self.Z_)
        if self.x_ == 1:
            self.mat_ = self.mat_.dot(self.X_)
        if self.x_ == self.z_ and self.x_ == 1:
            self.mat_ *= -1 * 1.0j

    def __str__(self):
        """
        """
        return np.array2string(self.mat_)


if __name__ == "__main__":
    print("I\n", Pauli(0, 0))
    print("X\n", Pauli(0, 1))
    print("Z\n", Pauli(1, 0))
    print("Y\n", Pauli(1, 1))
