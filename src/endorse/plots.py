import matplotlib.pyplot as plt
import numpy as np

def plot_field(points, values, cut=(1,2), file=None):
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_aspect('equal', 'box')
    #ax.set_ylim(-50, 50)
    #ax.set_xlim(-50, 50)

    axis_label=('x','y','z')
    ax.set_xlabel(f"${axis_label[cut[0]]}$", fontsize=20)
    ax.set_ylabel(f"${axis_label[cut[1]]}$", fontsize=20)
    assert points.shape[0] == len(values)
    if len(values) > 1000:
        subset = np.random.randint(0,len(values),size=1000)
        points = points[subset, :]
        values = values[subset]
    sc = ax.scatter(points[:, cut[0]], points[:, cut[1]], c=values, s=1, cmap=plt.cm.viridis)

    axx = ax.twinx()
    isort = np.argsort(points[:,cut[0]])
    X = points[isort,cut[0]]
    Y = values[isort]
    axx.plot(X, Y)
    # levels = np.array([])
    #c = ax.contourf(X, Y, porosity, cmap=plt.cm.viridis)
    cb = fig.colorbar(sc)
    if file:
        fig.savefig(file)
    else:
        plt.show()
