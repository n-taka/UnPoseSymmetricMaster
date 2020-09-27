#include "UnPoseSymmetricMaster.h"

#include <vector>
#include <string>
#include <filesystem>
#include "cpprest/json.h"
#include "readGoZFile.h"
#include "calculateEulerAnglesForSymmetrize.h"

#include "Eigen/Core"
#include "Eigen/Geometry"

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

extern "C" DLLEXPORT float initialize(char *someText, double optValue, char *outputBuffer, int optBuffer1Size, char *pOptBuffer2, int optBuffer2Size, char **zData)
{
	std::filesystem::path cache = std::filesystem::temp_directory_path();
	cache.append("UnPoseSymmetric.json");

	if (std::filesystem::exists(cache))
	{
		// erase older cache (for accidental collision in subtool name)
		std::filesystem::remove_all(cache);
	}

	return 1.0f;
}

extern "C" DLLEXPORT float isCached(char *someText, double optValue, char *outputBuffer, int optBuffer1Size, char *pOptBuffer2, int optBuffer2Size, char **zData)
{
	// find and parse tmp file
	web::json::value json;
	std::filesystem::path cache = std::filesystem::temp_directory_path();
	cache.append("UnPoseSymmetric.json");

	if (std::filesystem::exists(cache))
	{
		std::ifstream ifs(cache);
		json = web::json::value::parse(ifs);
		ifs.close();
	}

// parse parameters
#if defined(_WIN32) || defined(_WIN64)
	std::filesystem::path inputGoZFileName(someText);
#else
	// for macOS, zsc adds meaningless prefix
	std::string tmpStr(someText);
	std::filesystem::path inputGoZFileName(tmpStr.substr(2));
#endif
	utility::string_t key;
	{
		std::string fileName = inputGoZFileName.filename().string();
		key = utility::string_t(fileName.begin(), fileName.end());
	}

	if (json.has_field(key))
	{
		return 1.0f;
	}
	else
	{
		return 0.0f;
	}
}

extern "C" DLLEXPORT float getCachedRotation(char *someText, double optValue, char *outputBuffer, int optBuffer1Size, char *pOptBuffer2, int optBuffer2Size, char **zData)
{
	// find and parse tmp file
	web::json::value json;
	std::filesystem::path cache = std::filesystem::temp_directory_path();
	cache.append("UnPoseSymmetric.json");

	if (std::filesystem::exists(cache))
	{
		std::ifstream ifs(cache);
		json = web::json::value::parse(ifs);
		ifs.close();
	}

// parse parameters
#if defined(_WIN32) || defined(_WIN64)
	std::filesystem::path inputGoZFileName(someText);
#else
	// for macOS, zsc adds meaningless prefix
	std::string tmpStr(someText);
	std::filesystem::path inputGoZFileName(tmpStr.substr(2));
#endif
	utility::string_t key;
	{
		std::string fileName = inputGoZFileName.filename().string();
		key = utility::string_t(fileName.begin(), fileName.end());
	}

	Eigen::Matrix<double, 3, 1> eulerZYX, invEulerZYX, eulerZYX_Zb, invEulerZYX_Zb;
	Eigen::Matrix<double, 3, 1> translateXYZ;
	if (json.has_field(key))
	{
		// hit
		eulerZYX << json[key][_XPLATSTR("rotZ")].as_double(), json[key][_XPLATSTR("rotY")].as_double(), json[key][_XPLATSTR("rotX")].as_double();

		invEulerZYX << json[key][_XPLATSTR("invRotZ")].as_double(), json[key][_XPLATSTR("invRotY")].as_double(), json[key][_XPLATSTR("invRotX")].as_double();

		translateXYZ << json[key][_XPLATSTR("posX")].as_double(), json[key][_XPLATSTR("posY")].as_double(), json[key][_XPLATSTR("posZ")].as_double();
	}
	else
	{
		// no hit
		// calculate symmetrization
		Mesh<double, int> mesh;
		FromZ::readGoZFile(inputGoZFileName.string(), mesh.meshName, mesh.V, mesh.F, mesh.UV, mesh.VC, mesh.M, mesh.G);

		// remove temporary file
		std::filesystem::remove_all(inputGoZFileName);

		Eigen::Matrix<double, 3, 3> rotMatrix;
		if (calculateEulerAnglesForSymmetrize(mesh, rotMatrix))
		{
			eulerZYX = rotMatrix.eulerAngles(2, 1, 0);
			invEulerZYX.setZero();

			translateXYZ.setZero();
		}
		else
		{
			return 0.0f;
		}
	}

	union
	{
		char c[sizeof(float)];
		float f;
	} loader;
	web::json::value update;
	{
		float posX, posY, posZ, extraRotX;
		memcpy(loader.c, outputBuffer + sizeof(float) * 0, sizeof(float));
		extraRotX = loader.f;
		memcpy(loader.c, outputBuffer + sizeof(float) * 1, sizeof(float));
		posX = loader.f;
		memcpy(loader.c, outputBuffer + sizeof(float) * 2, sizeof(float));
		posY = loader.f;
		memcpy(loader.c, outputBuffer + sizeof(float) * 3, sizeof(float));
		posZ = loader.f;

		update[_XPLATSTR("posX")] = web::json::value(posX);
		update[_XPLATSTR("posY")] = web::json::value(posY);
		update[_XPLATSTR("posZ")] = web::json::value(posZ);
		// convert [degree] -> [rad]
		extraRotX *= (M_PI / 180.0);
		invEulerZYX(2) -= extraRotX;
	}

	eulerZYX_Zb = eulerZYX;
	eulerZYX_Zb *= 180.0;
	eulerZYX_Zb /= M_PI;
	// I don't fully understand, but eulerAngles in Eigen is left-handed?
	eulerZYX_Zb(0) *= -1;
	invEulerZYX_Zb = invEulerZYX;
	invEulerZYX_Zb *= 180.0;
	invEulerZYX_Zb /= M_PI;
	// I don't fully understand, but eulerAngles in Eigen is left-handed?
	invEulerZYX_Zb(0) *= -1;

	// pose -> symmetry
	loader.f = static_cast<float>(eulerZYX_Zb(0));
	memcpy(outputBuffer + sizeof(float) * 0, loader.c, sizeof(float));
	loader.f = static_cast<float>(eulerZYX_Zb(1));
	memcpy(outputBuffer + sizeof(float) * 1, loader.c, sizeof(float));
	loader.f = static_cast<float>(eulerZYX_Zb(2));
	memcpy(outputBuffer + sizeof(float) * 2, loader.c, sizeof(float));

	// symmetry -> pose
	loader.f = static_cast<float>(invEulerZYX_Zb(0));
	memcpy(outputBuffer + sizeof(float) * 3, loader.c, sizeof(float));
	loader.f = static_cast<float>(invEulerZYX_Zb(1));
	memcpy(outputBuffer + sizeof(float) * 4, loader.c, sizeof(float));
	loader.f = static_cast<float>(invEulerZYX_Zb(2));
	memcpy(outputBuffer + sizeof(float) * 5, loader.c, sizeof(float));

	loader.f = static_cast<float>(translateXYZ(0));
	memcpy(outputBuffer + sizeof(float) * 6, loader.c, sizeof(float));
	loader.f = static_cast<float>(translateXYZ(1));
	memcpy(outputBuffer + sizeof(float) * 7, loader.c, sizeof(float));
	loader.f = static_cast<float>(translateXYZ(2));
	memcpy(outputBuffer + sizeof(float) * 8, loader.c, sizeof(float));

	update[_XPLATSTR("rotZ")] = web::json::value(-invEulerZYX(0));
	update[_XPLATSTR("rotY")] = web::json::value(-invEulerZYX(1));
	update[_XPLATSTR("rotX")] = web::json::value(-invEulerZYX(2));
	update[_XPLATSTR("invRotZ")] = web::json::value(-eulerZYX(0));
	update[_XPLATSTR("invRotY")] = web::json::value(-eulerZYX(1));
	update[_XPLATSTR("invRotX")] = web::json::value(-eulerZYX(2));
	// write json to cache file.
	json[key] = update;
	std::ofstream ofs(cache);
	json.serialize(ofs);
	ofs.close();

	return 1.0f;
}
