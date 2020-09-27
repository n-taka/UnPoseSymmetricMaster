#include <vector>
#include <string>
#include <filesystem>
#include "readGoZFile.h"
#include "UnPoseSymmetricMaster.h"
#include "igl/writeOBJ.h"
#include "Eigen/Geometry"

void main(int argc, char *argv[])
{
    if (argc >= 2)
    {
        std::filesystem::path inputGoZFileName(argv[1]);
        char fileName[1024];
        char outputBuffer[1024];
        char dummyBuffer[1024];
        char **dummyBuffer1;
        union
        {
            char c[sizeof(float)];
            float f;
        } loader;
        loader.f = 1.0f;
        memcpy(outputBuffer + sizeof(float) * 3, loader.c, sizeof(float));
        loader.f = 2.0f;
        memcpy(outputBuffer + sizeof(float) * 4, loader.c, sizeof(float));
        loader.f = -4.5f;
        memcpy(outputBuffer + sizeof(float) * 5, loader.c, sizeof(float));

        sprintf(fileName, "%s", inputGoZFileName.string().c_str());
        initialize(fileName, 0.0f, outputBuffer, 0, dummyBuffer, 0, dummyBuffer1);
        getCachedRotation(fileName, 0.0f, outputBuffer, 0, dummyBuffer, 0, dummyBuffer1);
    }
    else
    {
        std::cerr << "Please drag-and-drop GoZ file to this .exe file." << std::endl;
    }

    return;
}