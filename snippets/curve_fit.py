# type: ignore
"""Curve fitting using scipy.optimize.curve_fit."""

from typing import Callable

import matplotlib.pyplot as plt
import numpy as np
from more_itertools import chunked
from scipy.optimize import curve_fit

data = np.array([0, 1, 3, 7, 3, 1, 0, 0.1, 0.3, 0.7, 0.3, 0.1, 0])

# Initial guess
p0 = [7, 3, 3]


# Define a model function
def gaussian(A: float, μ: float, σ: float) -> Callable[[float], float]:
    """Return a Gaussian function with the given parameters."""
    return lambda x: A * np.exp(-((x - μ) ** 2) / (2 * σ**2))


def f(x: float, *params: list[float]) -> float:
    """Return the sum of Gaussian functions with the given parameters."""
    if len(params) % 3:
        raise ValueError(
            f"Invalid number of parameters; should be divisible by 3: got {len(params)=}"
        )

    return np.sum((gaussian(A, μ, σ)(x) for A, μ, σ in chunked(params, 3)), axis=0)  # type:ignore


coeff, var_matrix = curve_fit(f, range(len(data)), data, p0=p0)

print(f"Coeff\n{coeff}")
print(f"Var Matrix\n{var_matrix}")

fig, ax = plt.subplots()
# Get the fitted curve
ax.plot(data)
ax.plot(f(data, *coeff))  # type:ignore

# Finally, lets get the fitting parameters, i.e. the mean and standard deviation:
print(f"Fitted mean = {coeff[1]}")
print(f"Fitted standard deviation {coeff[2]}")

plt.show()
