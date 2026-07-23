import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def load_and_slice_data(filename, z_plane=0.0):
    """
    Reads the Fortran output file and extracts a 2D slice near z = z_plane
    by finding the closest Z-coordinate PER LEVEL.
    """
    cols = ['Level', 'i', 'k', 'l', 'Xpos', 'Ypos', 'Zpos', 'Numerical', 'Analytical', 'Rel_Error']

    # Read whitespace-delimited file, skipping comments (#)
    df = pd.read_csv(filename, sep=r'\s+', comment='#', header=None, names=cols)

    slice_list = []
    # Process each level independently to get its specific slice at Z ~ z_plane
    for level, group in df.groupby('Level'):
        unique_z = group['Zpos'].unique()
        closest_z = unique_z[np.argmin(np.abs(unique_z - z_plane))]
        level_slice = group[np.isclose(group['Zpos'], closest_z, atol=1e-5)].copy()
        print(f"Level {int(level)}: Extracted slice at Z = {closest_z:.5f} with {len(level_slice)} points.")
        slice_list.append(level_slice)

    # Combine all slices into a single multi-resolution 2D point cloud
    df_slice = pd.concat(slice_list, ignore_index=True)
    return df_slice, z_plane

def nearest_neighbor_2d_numpy(x_pts, y_pts, values, grid_x, grid_y, chunk_size=1000):
    """
    Computes 2D nearest-neighbor interpolation using pure NumPy without SciPy.
    Uses chunking to keep memory overhead low while remaining fully vectorized.
    """
    Ny, Nx = grid_x.shape
    pts = np.column_stack((x_pts, y_pts))                     # Shape: (N_pts, 2)
    grid_flat = np.column_stack((grid_x.ravel(), grid_y.ravel())) # Shape: (Ny*Nx, 2)

    result_flat = np.zeros(len(grid_flat))

    # Process query points in chunks to prevent memory spikes
    for i in range(0, len(grid_flat), chunk_size):
        chunk = grid_flat[i:i + chunk_size]                    # Shape: (C, 2)

        # Calculate squared 2D distance between chunk query points and all input points
        # chunk[:, None, :] shape: (C, 1, 2), pts[None, :, :] shape: (1, N_pts, 2)
        dists_sq = np.sum((chunk[:, None, :] - pts[None, :, :]) ** 2, axis=2)

        # Find index of the minimum distance in any direction
        nearest_idx = np.argmin(dists_sq, axis=1)
        result_flat[i:i + chunk_size] = values[nearest_idx]

    return result_flat.reshape(Ny, Nx)

def resample_nearest_neighbor(df_slice, grid_resolution=300):
    """
    Resamples scattered multi-resolution point data onto a uniform 2D grid.
    """
    x = df_slice['Xpos'].values
    y = df_slice['Ypos'].values
    num = df_slice['Numerical'].values
    ana = df_slice['Analytical'].values
    err = df_slice['Rel_Error'].values

    # Determine global domain bounds
    x_min, x_max = x.min(), x.max()
    y_min, y_max = y.min(), y.max()

    # Create uniform target 2D mesh
    grid_x, grid_y = np.meshgrid(
        np.linspace(x_min, x_max, grid_resolution),
        np.linspace(y_min, y_max, grid_resolution)
    )

    # Nearest-neighbor mapping using pure NumPy
    num_grid = nearest_neighbor_2d_numpy(x, y, num, grid_x, grid_y)
    ana_grid = nearest_neighbor_2d_numpy(x, y, ana, grid_x, grid_y)
    err_grid = nearest_neighbor_2d_numpy(x, y, err, grid_x, grid_y)

    return grid_x, grid_y, num_grid, ana_grid, err_grid

def plot_contours(grid_x, grid_y, num_grid, ana_grid, err_grid, z_val):
    """
    Generates a 1x3 figure with filled contours from the nearest-neighbor grid.
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharex=True, sharey=True)

    vmin = min(num_grid.min(), ana_grid.min())
    vmax = max(num_grid.max(), ana_grid.max())

    # 1. Numerical Contour
    c1 = axes[0].contourf(grid_x, grid_y, num_grid, levels=20, cmap='viridis', vmin=vmin, vmax=vmax)
    axes[0].set_title("Numerical Solution")
    fig.colorbar(c1, ax=axes[0])

    # 2. Analytical Contour
    c2 = axes[1].contourf(grid_x, grid_y, ana_grid, levels=20, cmap='viridis', vmin=vmin, vmax=vmax)
    axes[1].set_title("Analytical Solution")
    fig.colorbar(c2, ax=axes[1])

    # 3. Relative Error Contour
    c3 = axes[2].contourf(grid_x, grid_y, err_grid, levels=20, cmap='magma')
    axes[2].set_title("Relative Error")
    fig.colorbar(c3, ax=axes[2])

    for ax in axes:
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_aspect('equal')

    fig.suptitle(f"Nearest-Neighbor Contour Plots at Z ≈ {z_val:.4f}", fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f"contour_nearest_z_{z_val:.4f}.png")
    print(f"Contour plot saved as contour_nearest_z_{z_val:.4f}.png")

def plot_colorplots(grid_x, grid_y, num_grid, ana_grid, err_grid, df_slice, z_val):
    """
    Generates a 1x3 figure using pcolormesh (Nearest-Neighbor Heatmaps).
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharex=True, sharey=True)

    vmin = min(num_grid.min(), ana_grid.min())
    vmax = max(num_grid.max(), ana_grid.max())

    # 1. Numerical Heatmap
    c1 = axes[0].pcolormesh(grid_x, grid_y, num_grid, cmap='coolwarm', vmin=vmin, vmax=vmax, shading='nearest')
    axes[0].set_title("Numerical Solution")
    fig.colorbar(c1, ax=axes[0])

    # 2. Analytical Heatmap
    c2 = axes[1].pcolormesh(grid_x, grid_y, ana_grid, cmap='coolwarm', vmin=vmin, vmax=vmax, shading='nearest')
    axes[1].set_title("Analytical Solution")
    fig.colorbar(c2, ax=axes[1])

    # 3. Relative Error Heatmap
    c3 = axes[2].pcolormesh(grid_x, grid_y, err_grid, cmap='inferno', shading='nearest')
    axes[2].set_title("Relative Error")
    fig.colorbar(c3, ax=axes[2])

    # Overlay multi-resolution grid locations
    x_pts = df_slice['Xpos'].values
    y_pts = df_slice['Ypos'].values
    for ax in axes:
        ax.scatter(x_pts, y_pts, s=1, c='black', alpha=0.15)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_aspect('equal')

    fig.suptitle(f"Nearest-Neighbor Heatmaps at Z ≈ {z_val:.4f}", fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f"heatmap_nearest_z_{z_val:.4f}.png")
    print(f"Heatmap saved as heatmap_nearest_z_{z_val:.4f}.png")

if __name__ == "__main__":
    filename = "global_it_40.dat"  # Change to your output filename
    target_z_plane = 0.0          # Target Z slice plane
    grid_resolution = 300         # Canvas density

    try:
        # 1. Load multi-level slice data
        df_slice, z_val = load_and_slice_data(filename, z_plane=target_z_plane)

        # 2. Resample onto regular grid using pure NumPy 2D nearest neighbor
        grid_x, grid_y, num_grid, ana_grid, err_grid = resample_nearest_neighbor(
            df_slice, grid_resolution=grid_resolution
        )

        # 3. Plot Contours
        plot_contours(grid_x, grid_y, num_grid, ana_grid, err_grid, z_val)

        # 4. Plot Heatmaps
        plot_colorplots(grid_x, grid_y, num_grid, ana_grid, err_grid, df_slice, z_val)

    except FileNotFoundError:
        print(f"Error: Output file '{filename}' not found. Make sure the file exists in the directory.")
