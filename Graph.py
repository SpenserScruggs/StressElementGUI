import matplotlib.pyplot as plt
import numpy as np
from math import cos, sin


def stressGraph(Fx, Fy, Fz, Vxy, Vxz, Vyz):

    # Stress Tensor
    A = [[Fx, Vxy, Vxz], 
         [Vxy, Fy, Vyz], 
         [Vxz, Vyz, Fz]]

    # calculating principal stresses
    eigval, eigvect = np.linalg.eigh(A)
    eigvect = np.round(eigvect, 5)
    eigval = np.round (eigval, 5)

    print("\nEigen Values: ")
    print(eigval, "\n")

    print("Eigenvectors (x, y, z): ")
    print(eigvect[:, 0])
    print(eigvect[:, 1])
    print(eigvect[:, 2])

    # generating points for the cube
    points = np.array([[1.0, 1.0, 1.0], 
                       [1.0, -1.0, 1.0], 
                       [-1.0, -1.0, 1.0], 
                       [-1.0, 1.0, 1.0], 
                       [1.0, 1.0, 1.0], 
                       [1.0, 1.0, -1.0], 
                       [1.0, -1.0, -1.0], 
                       [1.0, -1.0, 1.0],
                       [1.0, -1.0, -1.0],
                       [-1.0, -1.0, -1.0], 
                       [-1.0, -1.0, 1.0],
                       [-1.0, -1.0, -1.0],
                       [-1.0, 1.0, -1.0],
                       [-1.0, 1.0, 1.0],
                       [-1.0, 1.0, -1.0], 
                       [1.0, 1.0, -1.0]])

    original_points = points.copy()



    # z axis rotation
    theta = - np.arctan2(eigvect[1, 0], eigvect[0, 0])

    z_rot = np.matrix([[cos(theta), -sin(theta), 0], 
                       [sin(theta), cos(theta),  0], 
                       [0,      0,       1]])

    # xy rotation
    Rxy = np.arccos(eigvect[2, 0])

    magnitude = np.sqrt(eigvect[1, 0] ** 2 + eigvect[0, 0] ** 2)
    if magnitude == 0:
            magnitude = 1
    rx = round(eigvect[1, 0] / magnitude, 5)
    ry = - round(eigvect[0, 0] / magnitude, 5)
    rz = 0
    if rx == 0 and ry == 0:
        rx = 1
        ry = 0
        rz = 0
        Rxy_rot = 0

    Rxy_rot = np.matrix([[cos(Rxy) + (rx ** 2) * (1 - cos(Rxy)), rx * ry * (1 - cos(Rxy)) - rz * sin(Rxy), rx * rz * (1 - cos(Rxy)) + ry * sin(Rxy)],
                         [ry * rx * (1 - cos(Rxy)) + rz * sin(Rxy), cos(Rxy) + (ry ** 2) * (1 - cos(Rxy)), ry * rz * (1 - cos(Rxy)) - rx * sin(Rxy)],
                         [rz * rx * (1 - cos(Rxy)) - ry * sin(Rxy), rz * ry * (1 - cos(Rxy)) + rx * sin(Rxy), cos(Rxy) + (rz ** 2) * (1 - cos(Rxy))]])


    # xyz rotation
    Uxyz = - np.arccos(np.dot(eigvect[:, 1], [rx, ry, 0]) / (np.linalg.norm(eigvect[:, 1]) * np.linalg.norm([rx, ry, 0])))

    ux = eigvect[0, 0]
    uy = eigvect[1, 0]
    uz = eigvect[2, 0]

    Uxyz_rot = np.matrix([[cos(Uxyz) + (ux ** 2) * (1 - cos(Uxyz)), ux * uy * (1 - cos(Uxyz)) - uz * sin(Uxyz), ux * uz * (1 - cos(Uxyz)) + uy * sin(Uxyz)],
                         [uy * ux * (1 - cos(Uxyz)) + uz * sin(Uxyz), cos(Uxyz) + (uy ** 2) * (1 - cos(Uxyz)), uy * uz * (1 - cos(Uxyz)) - ux * sin(Uxyz)],
                         [uz * ux * (1 - cos(Uxyz)) - uy * sin(Uxyz), uz * uy * (1 - cos(Uxyz)) + ux * sin(Uxyz), cos(Uxyz) + (uz ** 2) * (1 - cos(Uxyz))]])

    # rotating points
    for i in range(len(points)):
        points[i] = np.dot(points[i], z_rot)
        points[i] = np.dot(points[i], Rxy_rot)
        points[i] = np.dot(points[i], Uxyz_rot)
        

    # assigning P1 P2 and P3
    stress = [0, 0, 0]
    for i in range(len(eigval)):
        if eigval[i] == max(eigval):
            stress[i] = "P1"
        elif eigval[i] == min(eigval):
            stress[i] = "P3"
        else:
            stress[i] = "P2"

    # drawing the graph
    ax = plt.axes(projection="3d")

    ax.plot(original_points[:, 0] * np.sqrt(3), original_points[:, 1] * np.sqrt(3), original_points[:, 2] * np.sqrt(3), color="red", linestyle='dashed')
    ax.plot(points[:, 0], points[:, 1], points[:, 2], color="blue")

    ax.quiver(eigvect[0, 0], eigvect[1, 0], eigvect[2, 0], eigvect[0, 0] * 3, eigvect[1, 0] * 3, eigvect[2, 0] * 3, color="blue")
    ax.quiver(eigvect[0, 1], eigvect[1, 1], eigvect[2, 1], eigvect[0, 1] * 3, eigvect[1, 1] * 3, eigvect[2, 1] * 3, color="blue")
    ax.quiver(eigvect[0, 2], eigvect[1, 2], eigvect[2, 2], eigvect[0, 2] * 3, eigvect[1, 2] * 3, eigvect[2, 2] * 3, color="blue")

    ax.text(eigvect[0, 0] * 4, eigvect[1, 0] * 4, eigvect[2, 0] * 4, stress[0] + ": " + str(eigval[0]))
    ax.text(eigvect[0, 1] * 4, eigvect[1, 1] * 4, eigvect[2, 1] * 4, stress[1] + ": " + str(eigval[1]))
    ax.text(eigvect[0, 2] * 4, eigvect[1, 2] * 4, eigvect[2, 2] * 4, stress[2] + ": " + str(eigval[2]))

    ax.quiver(np.sqrt(3), 0, 0, 3, 0, 0, color="red", linestyle="dashed")
    ax.quiver(0, np.sqrt(3), 0, 0, 3, 0, color="red", linestyle="dashed")
    ax.quiver(0, 0, np.sqrt(3), 0, 0, 3,  color="red", linestyle="dashed")

    ax.text(4.5, 0, 0, r'$\sigma$x: ' + str(Fx))
    ax.text(0, 4.5, 0, r'$\sigma$y: ' + str(Fy))
    ax.text(0, 0, 4.5, r'$\sigma$z: ' + str(Fz))

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    ax.set_zlim3d(-2.5, 2.5)
    ax.set_xlim3d(-3, 3)
    ax.set_ylim3d(-3, 3)

    ax.text2D(-0.15, 1.1, r'$\sigma$ X: ' + str(Fx), transform=ax.transAxes)
    ax.text2D(-0.15, 1.05, r'$\sigma$ Y: ' + str(Fy), transform=ax.transAxes)
    ax.text2D(-0.15, 1, r'$\sigma$ Z: ' + str(Fz), transform=ax.transAxes)
    ax.text2D(-0.15, 0.95, r'$\tau$ xy: ' + str(Vxy), transform=ax.transAxes)
    ax.text2D(-0.15, 0.9, r'$\tau$ xz: ' + str(Vxz), transform=ax.transAxes)
    ax.text2D(-0.15, 0.85, r'$\tau$ yz: ' + str(Vyz), transform=ax.transAxes)

    ax.text2D(0.95, 1.1, stress[0] + ": " + str(round(eigval[0], 3)), transform=ax.transAxes)
    ax.text2D(0.95, 1.05, stress[1] + ": " + str(round(eigval[1], 3)), transform=ax.transAxes)
    ax.text2D(0.95, 1, stress[2] + ": " + str(round(eigval[2], 3)), transform=ax.transAxes)

    plt.show()