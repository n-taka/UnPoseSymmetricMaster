#include <vector>
#include <string>
#include <filesystem>
#include "readGoZFile.h"
#include "calculateRotation.h"
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
        calculateRotation(meshIn, eulerXZY, translateXYZ);

        std::cout << eulerXZY << std::endl;
        std::cout << translateXYZ << std::endl;
        V *= Eigen::AngleAxisd(eulerXZY(0), Eigen::Vector3d::UnitX()).toRotationMatrix();
        igl::writeOBJ("rot1.obj", V, F);
        V *= Eigen::AngleAxisd(eulerXZY(1), Eigen::Vector3d::UnitZ()).toRotationMatrix();
        igl::writeOBJ("rot2.obj", V, F);
        V *= Eigen::AngleAxisd(eulerXZY(2), Eigen::Vector3d::UnitY()).toRotationMatrix();
        igl::writeOBJ("rot3.obj", V, F);
    }
    else
    {
        std::cerr << "Please drag-and-drop GoZ file to this .exe file." << std::endl;
    }

    return;
}