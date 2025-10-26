#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

#include <Eigen/Eigen>

#include <iostream>

int main(int argc, char *argv[])
{
    polyscope::init();

    Eigen::MatrixX3f vertices{
        {0, 0, 0},
        {1, 0, 0},
        {0.5f, std::sqrt(3.f) * 0.5f, 0},
    };

    Eigen::MatrixX3i faces{
        {0, 1, 2},
    };

    Eigen::Matrix3f colors{
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1},
    };

    polyscope::registerSurfaceMesh("triangle", vertices, faces)
        ->addVertexColorQuantity("vColors", colors)
        ->setEnabled(true);

    polyscope::show();

    return 0;
}
