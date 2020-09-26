#include "UnPoseSymmetricMaster.h"

#include <vector>
#include <string>
#include <filesystem>
#include "readGoZFile.h"
#include "calculateEulerAnglesForSymmetrize.h"
#include "writeGoZFile.h"

////
// implementation
////
#if defined(_WIN32) || defined(_WIN64)
#define DLLEXPORT __declspec(dllexport)
#else
#define DLLEXPORT __attribute__((visibility("default")))
#endif

extern "C" DLLEXPORT float version(char *someText, double optValue, char *outputBuffer, int optBuffer1Size, char *pOptBuffer2, int optBuffer2Size, char **zData)
{
	return 1.0f;
}

extern "C" DLLEXPORT float computeRotation(char *someText, double optValue, char *outputBuffer, int optBuffer1Size, char *pOptBuffer2, int optBuffer2Size, char **zData)
{
	// parse parameters
	std::filesystem::path inputGoZFileName(someText);
	union
	{
		char c[sizeof(float)];
		float f;
	} loader;

	{
		Mesh<double, int> mesh;
		FromZ::readGoZFile(inputGoZFileName.string(), mesh.meshName, mesh.V, mesh.F, mesh.UV, mesh.VC, mesh.M, mesh.G);

		Eigen::Matrix<double, 1, Eigen::Dynamic> eulerXYZ;
		Eigen::Matrix<double, 1, Eigen::Dynamic> translateXYZ;
		if (calculateEulerAnglesForSymmetrize(mesh, eulerXYZ, translateXYZ))
		{
			eulerXYZ *= 180.0;
			eulerXYZ /= M_PI;
			// I don't fully understand, but eulerAngles in Eigen is left-handed?
			eulerXYZ(2) *= -1;
			// std::ofstream log("log.txt");
			// log << eulerXYZ << std::endl;
			// log.close();
			loader.f = static_cast<float>(eulerXYZ(0));
			memcpy(outputBuffer, loader.c, sizeof(float));
			loader.f = static_cast<float>(eulerXYZ(1));
			memcpy(outputBuffer + sizeof(float), loader.c, sizeof(float));
			loader.f = static_cast<float>(eulerXYZ(2));
			memcpy(outputBuffer + sizeof(float) * 2, loader.c, sizeof(float));
			loader.f = static_cast<float>(translateXYZ(0));
			memcpy(outputBuffer + sizeof(float) * 3, loader.c, sizeof(float));
			loader.f = static_cast<float>(translateXYZ(1));
			memcpy(outputBuffer + sizeof(float) * 4, loader.c, sizeof(float));
			loader.f = static_cast<float>(translateXYZ(2));
			memcpy(outputBuffer + sizeof(float) * 5, loader.c, sizeof(float));
			return 1.0f;
		}
		else
		{
			return 0.0f;
		}
	}
}
