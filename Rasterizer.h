#pragma once
#include <vector>
#include "Eigen/Eigen"
#include "scan_line_z-buffer.h"

//std::tuple<double, double, double> computeBarycentric2D(double x, double y, const Eigen::Vector4d* v) {
std::vector<double> computeBarycentric2D(double x, double y, const Eigen::Vector3d* v) {
	double c1 = (x*(v[1].y() - v[2].y()) + (v[2].x() - v[1].x())*y + v[1].x()*v[2].y() - v[2].x()*v[1].y()) / (v[0].x()*(v[1].y() - v[2].y()) + (v[2].x() - v[1].x())*v[0].y() + v[1].x()*v[2].y() - v[2].x()*v[1].y());
	double c2 = (x*(v[2].y() - v[0].y()) + (v[0].x() - v[2].x())*y + v[2].x()*v[0].y() - v[0].x()*v[2].y()) / (v[1].x()*(v[2].y() - v[0].y()) + (v[0].x() - v[2].x())*v[1].y() + v[2].x()*v[0].y() - v[0].x()*v[2].y());
	double c3 = (x*(v[0].y() - v[1].y()) + (v[1].x() - v[0].x())*y + v[0].x()*v[1].y() - v[1].x()*v[0].y()) / (v[2].x()*(v[0].y() - v[1].y()) + (v[1].x() - v[0].x())*v[2].y() + v[0].x()*v[1].y() - v[1].x()*v[0].y());
	return std::vector<double>{c1, c2, c3};
}

Eigen::Vector2d interpolate(double alpha, double beta, double gamma, const Eigen::Vector2d& vert1, const Eigen::Vector2d& vert2, const Eigen::Vector2d& vert3, double weight)
{
	auto u = (alpha * vert1[0] + beta * vert2[0] + gamma * vert3[0]);
	auto v = (alpha * vert1[1] + beta * vert2[1] + gamma * vert3[1]);

	u /= weight;
	v /= weight;

	return Eigen::Vector2d(u, v);
}

Eigen::Vector3d interpolate(double alpha, double beta, double gamma, const Eigen::Vector3d& vert1, const Eigen::Vector3d& vert2, const Eigen::Vector3d& vert3, double weight)
{
	return (alpha * vert1 + beta * vert2 + gamma * vert3) / weight;

	auto u = (alpha * vert1[0] + beta * vert2[0] + gamma * vert3[0]);
	auto v = (alpha * vert1[1] + beta * vert2[1] + gamma * vert3[1]);
	auto w = (alpha * vert1[2] + beta * vert2[2] + gamma * vert3[2]);

	u /= weight;
	v /= weight;
	w /= weight;

	return Eigen::Vector3d(u, v, w);
}

struct Rasterizer
{
	int width;
	int height;
	std::vector<Eigen::Vector3d> Color_buffer;
	std::vector<double> z_buffer;
	
	Eigen::Matrix4d model;
	Eigen::Matrix4d view;
	Eigen::Matrix4d projection;
	cv::Mat tex;
	int tex_width;
	int tex_height;

	Rasterizer(int width, int height)
	{
		this->width = width;
		this->height = height;
		Color_buffer.clear();
		Color_buffer.resize(width * height);
		z_buffer.clear();
		z_buffer.resize(width * height, -DBL_MAX);
	}
	int get_index(int x, int y)
	{
		return (height - y - 1) * width + x;
	}
	void update_Color(int x, int y, Eigen::Vector3d color)
	{
		Color_buffer[get_index(x, y)] = color;
	}
	void Draw(nodeActivedPolygon* polygon, double y)
	{
		double L = polygon->edge->xl;
		double R = polygon->edge->xr;
		if (polygon->edge->xl > polygon->edge->xr)
		{
			printf("error: xl > xr\n");
			return;
		}
		double zx = polygon->edge->zl;
		zx = zx + (int(polygon->edge->xl) - polygon->edge->xl) * polygon->edge->dzx;
		zx = zx + 0.5 * polygon->edge->dzx;
		zx = zx - 0.5 * polygon->edge->dzy;
		for (int x = polygon->edge->xl; x <= polygon->edge->xr; x++)
		{
			int index = get_index(x, y);
			if (zx > z_buffer[index])
			{
				//Point* Triangle = polygon->Triangle_node->point;
				nodeClassifiedPolygon* Triangle = polygon->Triangle_node;
				//std::vector<Eigen::Vector3d> v{ Triangle[0].v , Triangle[1].v, Triangle[2].v};
				Eigen::Vector3d v[3] = { Triangle->point[0].v , Triangle->point[1].v, Triangle->point[2].v };
				std::vector<double> res = computeBarycentric2D(x + 0.5, y + 0.5, v);
				double alpha = res[0], beta = res[1], gamma = res[2];
				Eigen::Vector2d interpolated_texcoords = interpolate(alpha, beta, gamma, Triangle->tex_coor[0], Triangle->tex_coor[1], Triangle->tex_coor[2], 1);
				//Eigen::Vector3d interpolated_texcolors = interpolate(alpha, beta, gamma, Triangle->tex_color[0], Triangle->tex_color[1], Triangle->tex_color[2], 1);
				//std::cout << interpolated_texcoords[0] << " " << interpolated_texcoords[1] << std::endl;

				Eigen::Vector2d pos = interpolated_texcoords;
				pos[0] = pos[0] - floor(pos[0]);
				pos[1] = pos[1] - floor(pos[1]);

				pos[0] = pos[0] * tex_height;
				pos[1] = (1 - pos[1]) * tex_width;

				cv::Vec3b color = tex.at<cv::Vec3b>(int(pos[1]), int(pos[0]));
				Eigen::Vector3d interpolated_texcolors = { double(color[0]), double(color[1]), double(color[2]) };
				z_buffer[index] = zx;
				//update_Color(x, y, polygon->color);
				update_Color(x, y, interpolated_texcolors);
				//update_Color(x, y, Eigen::Vector3d(255, 0, 0));
			}
			zx = zx + polygon->edge->dzx;
		}
	}

	void set_model(const Eigen::Matrix4d& m)
	{
		model = m;
	}
	void set_view(const Eigen::Matrix4d& v)
	{
		view = v;
	}
	void set_projection(const Eigen::Matrix4d& p)
	{
		projection = p;
	}

};
