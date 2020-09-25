#include "UnPoseSymmetricMaster.h"

#include <vector>
#include <string>
#include <filesystem>
#include "readGoZFile.h"
#include "calculateSymmetricPlane.h"
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
	memcpy(loader.c, outputBuffer, sizeof(float));
	double height = static_cast<float>(loader.f);
	memcpy(loader.c, outputBuffer + sizeof(float), sizeof(float));
	double minClearance = static_cast<float>(loader.f);
	memcpy(loader.c, outputBuffer + sizeof(float) * 2, sizeof(float));
	double maxClearance = static_cast<float>(loader.f);
	double scale = -1.0;

	{
		Mesh<double, int> mesh;
		FromZ::readGoZFile(inputGoZFileName.string(), mesh.meshName, mesh.V, mesh.F, mesh.UV, mesh.VC, mesh.M, mesh.G);
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> V;
		igl::list_to_matrix(mesh.V, V);
		double maxY, minY;
		maxY = V.col(1).maxCoeff();
		minY = V.col(1).minCoeff();
		if (maxY <= minY)
		{
			return 0.0f;
		}
		scale = height / (maxY - minY);
	}

	std::filesystem::path inputGoZFileDir = inputGoZFileName.parent_path();
	std::vector<Mesh<double, int>> meshIn, meshOut;

	for (const std::filesystem::directory_entry &x : std::filesystem::directory_iterator(inputGoZFileDir))
	{
		if (x.path().extension() == ".GoZ")
		{
			Mesh<double, int> mesh;
			FromZ::readGoZFile(x.path().string(), mesh.meshName, mesh.V, mesh.F, mesh.UV, mesh.VC, mesh.M, mesh.G);
			meshIn.push_back(mesh);
		}
	}

	calculateClearance(meshIn, minClearance, maxClearance, scale, meshOut);

	int meshIdx = 0;
	for (const std::filesystem::directory_entry &x : std::filesystem::directory_iterator(inputGoZFileDir))
	{
		if (x.path().extension() == ".GoZ")
		{
			std::filesystem::path output = x.path().parent_path();
			output /= x.path().stem();
			output += ".GoZ";

			const Mesh<double, int> &mesh = meshOut.at(meshIdx++);

			FromZ::writeGoZFile(output.string(), mesh.meshName, mesh.V, mesh.F, mesh.UV, mesh.VC, mesh.M, mesh.G);
		}
	}

	return 1.0f;
}
