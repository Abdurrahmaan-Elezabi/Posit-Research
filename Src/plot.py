import pandas as pd
import matplotlib.pyplot as plt

filename = "./plots/test2.csv"

with open(filename, "r") as f:
    lines = [line.strip() for line in f if line.strip()]

for run_idx, line in enumerate(lines, start=1):
    residuals = [float(x) for x in line.split(",")]
    iterations = list(range(len(residuals)))
    plt.plot(iterations, residuals, label=f"Run {run_idx}")

plt.xlabel("Iteration")
plt.ylabel("Residual")
plt.title("Conjugate Gradient Convergence")
plt.yscale("log")
plt.legend()
plt.grid(True)
plt.show()