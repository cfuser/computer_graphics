#pragma once
#include<Eigen/eigen>
#include "scan_line_z-buffer.h"
struct Quadtree
{
	int width;
	int height;
	double *z_depth;
	Eigen::Vector3d *color;
	int n;
	Quadtree(int width, int height)
	{
		this->width = width;
		this->height = height;
		n = width * height;
		z_depth = new double[n * 6];
		color = new Eigen::Vector3d[n * 6];

		for (int i = 0; i < n * 6; i++)
		{
			z_depth[i] = -DBL_MAX;
			color[i] = { 0, 0, 0 };
		}
	}

	void insideTriangle(double x, double y, nodeActivedPolygon Polygon, int count)
	{

	}
	void query(int width_L, int width_R, int height_L, int height_R, int x)
	{

	}
	void update(int width_L, int width_R, int height_L, int height_R, int x)
	{
		double z_min = DBL_MAX;
		int count = 0;
		z_min = insideTriangle(width_L, height_L, triangle, &count);
		z_min = 0;
		if (z_min < z_depth[x]) return;
		int width_mid = (width_L + width_R) / 2;
		int height_mid = (height_L + height_R) / 2;

		update(width_L, width_mid, height_L, height_mid, 4 * x + 1);
		update(width_mid, width_R, height_L, height_mid, 4 * x + 2);
		update(width_L, width_mid, height_mid, height_R, 4 * x + 3);
		update(width_mid, width_R, height_mid, height_R, 4 * x + 4);
		z_depth[i] = min(z_depth[4 * x + 1], z_depth[4 * x + 2]);
		z_depth[i] = min(z_depth[4 * x + 3], z_depth[i]);
		z_depth[i] = min(z_depth[4 * x + 4], z_depth[i]);

	}
};