#include <vector>
#include <string>
#include <filesystem>
#include "readGoZFile.h"
#include "calculateSymmetricPlane.h"
#include "igl/writeOBJ.h"
#include "Eigen/Geometry"

void main(int argc, char *argv[])
{
    if (argc >= 2)
    {
        std::filesystem::path inputGoZFileName(argv[1]);
        Mesh<double, int> meshIn;

        FromZ::readGoZFile(inputGoZFileName.string(), meshIn.meshName, meshIn.V, meshIn.F, meshIn.UV, meshIn.VC, meshIn.M, meshIn.G);

        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> V;
        Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> F;
        igl::list_to_matrix(meshIn.V, V);
        igl::polygon_mesh_to_triangle_mesh(meshIn.F, F);

        igl::writeOBJ("in.obj", V, F);

        Eigen::Matrix<double, 1, Eigen::Dynamic> eulerXZY;
        Eigen::Matrix<double, 1, Eigen::Dynamic> translateXYZ;
        calculateSymmetricPlane(meshIn, eulerXZY, translateXYZ);

    }
    else
    {
        std::cerr << "Please drag-and-drop GoZ file to this .exe file." << std::endl;
    }

    return;
}