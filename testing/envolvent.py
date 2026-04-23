from scipy.optimize import linprog
import numpy as np
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt


def closest_point_on_segment(point, segment_start, segment_end):
    segment_vector = segment_end - segment_start
    point_vector = point - segment_start

    s = np.dot(point_vector.flatten(), segment_vector.flatten()) / np.dot(segment_vector.flatten(), segment_vector.flatten())
    if s < 0:
        closest_point = segment_start
    elif s > 1:
        closest_point = segment_end
    else:
        closest_point = segment_start + s * segment_vector

    return closest_point


def distance_from_convex_hull(X, x_target):
    num_iter = 20
    eps = 10e-6
    # initialize
    norms = np.linalg.norm(X - x_target, axis = 0)
    index = np.argmin(norms)
    p = X[:, [index]]
    t = p

    for i in range(num_iter):
        prods = (X - t).T @ (x_target - t)
        index = np.argmax(prods)
        p = X[:, [index]]
        if np.linalg.norm(p - t) < eps:
            break
        t = closest_point_on_segment(x_target, t, p)

    return np.linalg.norm(x_target - t)


def optimal_point_to_add(X, target_points):
    num_target_points = target_points.shape[1]
    E = np.zeros([num_target_points, num_target_points])
    for i in range(num_target_points):
        for j in range(num_target_points):
            extended_points = np.hstack([X, target_points[:, [j]]])
            x_target = target_points[:, [i]]
            E[i, j] = distance_from_convex_hull(extended_points, x_target)

    max_values = np.linalg.norm(E, axis=0)
    best_index = np.argmin(max_values)

    return best_index

def optimal_point_to_remove(X):
    num_points = X.shape[1]
    E = np.zeros(num_points)
    for i in range(num_points):
        for j in range(num_points):
            x_target = X[:, [i]]
            X_reduced = np.delete(X, i, axis=1)
            E[i] = distance_from_convex_hull(X_reduced, x_target)

    # Find the min value
    best_index = np.argmin(E)

    return best_index


def remove_points(X, eps):
    num_points = X.shape[1]
    for i in range(num_points):
        x_target = X[:, [i]]
        X_reduced = np.delete(X, i, axis=1)
        if distance_from_convex_hull(X_reduced, x_target) < eps:
            X = X_reduced
    return X


def approximate_convex_hull(points, K, eps):
    num_points = points.shape[1]
    selection_mask = np.array([False] * num_points)
    selection_mask[0] = True
    max_dist = np.inf

    # Add until reached the desired number
    while np.sum(selection_mask) < K:
        current_points = points[:, selection_mask]
        target_points = points[:, ~selection_mask]
        index_mapping = [i for i in range(num_points) if not selection_mask[i]]
        index = optimal_point_to_add(current_points, target_points)
        index_to_insert = index_mapping[index]
        selection_mask[index_to_insert] = True

    # Remove and add to stabilize
    for i in range(K):
        # Remove
        current_points = points[:, selection_mask]
        index_mapping = [i for i in range(num_points) if selection_mask[i]]
        remove_index = optimal_point_to_remove(current_points)
        index_to_insert = index_mapping[remove_index]
        selection_mask[index_to_insert] = False

        # Add
        current_points = points[:, selection_mask]
        target_points = points[:, ~selection_mask]
        index_mapping = [i for i in range(num_points) if not selection_mask[i]]
        add_index = optimal_point_to_add(current_points, target_points)
        index_to_insert = index_mapping[add_index]
        selection_mask[index_to_insert] = True

        if remove_index==add_index:
            break


    return selection_mask

if __name__ == "__main__":
    dim = 3
    # points = np.random.randn(dim, 25)
    K = 4
    eps = 0.0001

    indices = np.random.choice(np.arange(dim), size=2, replace=False)
    indices.sort()

    selection_mask = approximate_convex_hull(points, K, eps)
    selected_points = points[:, selection_mask]

    hull = ConvexHull(selected_points[indices, :].T)

    plt.scatter(points[indices[0],:], points[indices[1], :], c='b')
    plt.fill(selected_points[indices[0], hull.vertices], selected_points[indices[1], hull.vertices], 'r', alpha=0.5)

    plt.show()


