#pragma once
#include <vector>
#include "Eigen/Eigen"

struct Point
{
	Eigen::Vector3d v;
	Eigen::Vector3d vn;
	Eigen::Vector2d vt;
};
struct Mesh
{
	int v[3];
	int vn[3];
	int vt[3];
	// int order[3];
};

struct TriangleMesh
{
	//std::vector<Point> Points;
	std::vector<Eigen::Vector3d> v;
	std::vector<Eigen::Vector3d> vn;
	std::vector<Eigen::Vector2d> vt;
	std::vector<Mesh> Meshes;

	void LoadFile(std::string file_name);
};