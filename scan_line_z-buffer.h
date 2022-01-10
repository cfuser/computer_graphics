#pragma once
#pragma once
#include "TriangleMesh.h"

struct nodeClassifiedEdge {
	double x;
	double dx;
	int dy;
	int id;
	nodeClassifiedEdge *next;
	int top_ID, bot_ID;
	void setEdge(Point a, Point b, int idx)
	{
		x = a.v[0];// std::min(a.v[0], b.v[0]);
		dx = -1 / ((b.v[1] - a.v[1]) / (b.v[0] - a.v[0]));
		//dx = -(b.v[0] - a.v[0]) / (b.v[1] - a.v[1]);
		dy = ceil(a.v[1]) - floor(b.v[1]);
		//dy = std::abs(a.v[1] - b.v[1]);
		id = idx;
		next = NULL;
	}
	nodeClassifiedEdge(Point a, Point b, int idx, int top_ID, int bot_ID)
	{
		x = a.v[0];// std::min(a.v[0], b.v[0]);
		dx = -1 / ((b.v[1] - a.v[1]) / (b.v[0] - a.v[0]));
		//dx = -(b.v[0] - a.v[0]) / (b.v[1] - a.v[1]);
		dy = ceil(a.v[1]) - floor(b.v[1]);
		// dy = std::abs(a.v[1] - b.v[1]);
		id = idx;
		next = NULL;
		this->top_ID = top_ID;
		this->bot_ID = bot_ID;

	}
};

struct nodeClassifiedPolygon {
	/*
	attritubes
	a, b, c, d is the coefficient of face,
	id is the identify of face,
	dy is the number of scan line,
	color is the color of face,
	next is the point to next face,
	edge is the vector of edge,

	functions
	nodeClassifiedPolygon() to initialize

	*/
	double a, b, c, d;
	int id;
	int dy;
	Eigen::Vector3d color;
	Point point[3];
	nodeClassifiedPolygon *next;
	std::vector<nodeClassifiedEdge*> edge;
	Eigen::Vector2d tex_coor[3];
	Eigen::Vector3d tex_color[3];
	// std::list<nodeClassifiedPolygon> list;
	nodeClassifiedPolygon() {}
	nodeClassifiedPolygon(Point point[3], int idx)
	{
		if (point[0].v[1] < point[1].v[1]) std::swap(point[0], point[1]);
		if (point[0].v[1] < point[2].v[1]) std::swap(point[0], point[2]);
		if (point[1].v[1] < point[2].v[1]) std::swap(point[1], point[2]);
		//if (point[1].v[0] > point[2].v[0]) std::swap(point[1], point[2]);

		this->point[0] = point[0];
		this->point[1] = point[1];
		this->point[2] = point[2];
		Eigen::Vector3d vec01 = point[1].v - point[0].v;
		Eigen::Vector3d vec02 = point[2].v - point[0].v;
		Eigen::Vector3d vn = vec01.cross(vec02);
		vn.normalize();
		a = vn[0];
		b = vn[1];
		c = vn[2];
		d = -point[0].v.dot(vn);
		//printf("%d %lf %lf %lf %lf\n", idx, a, b, c, d);
		double on_plane[3];
		on_plane[0] = a * point[0].v[0] + b * point[0].v[1] + c * point[0].v[2] + d;
		on_plane[1] = a * point[1].v[0] + b * point[1].v[1] + c * point[1].v[2] + d;
		on_plane[2] = a * point[2].v[0] + b * point[2].v[1] + c * point[2].v[2] + d;
		
		if (abs(on_plane[0]) > 1e-6 || abs(on_plane[0]) > 1e-6 || abs(on_plane[2]) > 1e-6)
		{
			printf("%d %lf %lf %lf %lf\n", idx, a, b, c, d);
			printf("%lf %lf %lf\n", on_plane[0], on_plane[1], on_plane[2]);
			system("pause");
		}
		if (point[0].v[1] < point[1].v[1] || point[0].v[1] < point[2].v[1])
			system("pause");
		id = idx;
		double y_max = ceil(point[0].v[1]);// ceil(std::max(std::max(point[0].v[1], point[1].v[1]), point[2].v[1]));
		double y_min = floor(std::min(point[1].v[1], point[2].v[1]));// floor(std::max(std::min(point[0].v[1], point[1].v[1]), point[2].v[1]));
		dy = y_max - y_min;

		//if (point[0].v[1] != point[1].v[1])
		if (ceil(point[0].v[1]) > ceil(point[1].v[1]))
		edge.push_back(new nodeClassifiedEdge(point[0], point[1], idx, 0, 1));

		//if (point[0].v[1] != point[2].v[1])
		if (ceil(point[0].v[1]) > ceil(point[2].v[1]))
			edge.push_back(new nodeClassifiedEdge(point[0], point[2], idx, 0, 2));

		//if (point[1].v[1] != point[2].v[1])
		if (ceil(point[1].v[1]) > ceil(point[2].v[1]))
			edge.push_back(new nodeClassifiedEdge(point[1], point[2], idx, 1, 2));
		
		if (edge.size() < 2) return;
		if (edge[0]->x > edge[1]->x || edge[0]->x == edge[1]->x && edge[0]->dx > edge[1]->dx)
		{
			auto temp = edge[0];
			edge[0] = edge[1];
			edge[1] = temp;
			//std::swap(edge[0], edge[1]);
		}
		next = NULL;
	}
};

struct nodeActivedEdge {
	double xl, xr, dxl, dyl, dxr, dyr;
	double zl;
	double dzx, dzy;
	int id;
	nodeActivedEdge *next;

	nodeClassifiedPolygon *Triangle_node;

	nodeActivedEdge(nodeClassifiedPolygon* nodePolygon, double y)
	{
		xl = nodePolygon->edge[0]->x;
		dxl = nodePolygon->edge[0]->dx;
		dyl = nodePolygon->edge[0]->dy;

		xr = nodePolygon->edge[1]->x;
		dxr = nodePolygon->edge[1]->dx;
		dyr = nodePolygon->edge[1]->dy;
		if (xl > xr)
			system("pause");
		Triangle_node = nodePolygon;
		dzx = -Triangle_node->a / Triangle_node->c;
		dzy = Triangle_node->b / Triangle_node->c;
		id = Triangle_node->id;
		zl = (-Triangle_node->d - Triangle_node->a * xl - Triangle_node->b * y) / Triangle_node->c;
		next = NULL;
	}

	void Print()
	{
		printf("xl, dxl, dyl = %lf %lf %lf\n", xl, dxl, dyl);
		printf("xr, dxr, dyr = %lf %lf %lf\n", xr, dxr, dyr);
		printf("zl, dzx, dzy = %lf %lf %lf\n", zl, dzx, dzy);
		printf("%d\n", id);
	}
};

struct nodeActivedPolygon {
	double a, b, c, d;
	int id;
	int dy;
	Eigen::Vector3d color;
	nodeActivedPolygon *next;
	nodeClassifiedPolygon *Triangle_node;
	//std::vector<nodeClassifiedEdge*> edge;
	nodeActivedEdge *edge;

	nodeActivedPolygon(nodeClassifiedPolygon* nodePolygon, double y)
	{
		a = nodePolygon->a;
		b = nodePolygon->b;
		c = nodePolygon->c;
		d = nodePolygon->d;
		id = nodePolygon->id;
		dy = nodePolygon->dy;
		color = nodePolygon->color;
		edge = new nodeActivedEdge(nodePolygon, y);
		Triangle_node = nodePolygon;
		next = NULL;
	}

	void update()
	{
		dy = dy - 1;
		if (dy == 0) return ;
		edge->dyl = edge->dyl - 1;
		edge->dyr = edge->dyr - 1;
		
		if (edge->dyl == 0)
		{
			edge->xl = Triangle_node->edge[2]->x;
			edge->dxl = Triangle_node->edge[2]->dx;
			// edge->dyl = Triangle_node->edge[2]->dy - 1;
			// edge->zl = (-Triangle_node->d - Triangle_node->a * edge->xl - Triangle_node->b * edge->y) / Triangle_node->c;
		}
		int bot_ID_l = edge->Triangle_node->edge[0]->bot_ID;
		double depth_l = 1;
		if (edge->dyl == 1)
		{
			depth_l = edge->Triangle_node->point[bot_ID_l].v[1];
			depth_l = ceil(depth_l) - depth_l;
		}
		//depth_l = 1;
		edge->xl = edge->xl + edge->dxl * depth_l;
		
		if (edge->dyr == 0)
		{
			edge->xr = Triangle_node->edge[2]->x;
			edge->dxr = Triangle_node->edge[2]->dx;
			edge->dyr = Triangle_node->edge[2]->dy - 1;
		}
		int bot_ID_r = edge->Triangle_node->edge[0]->bot_ID;
		double depth_r = 1;
		if (edge->dyr == 1)
		{
			depth_r = edge->Triangle_node->point[bot_ID_r].v[1];
			depth_r = ceil(depth_r) - depth_r;
		}
		//depth_r = 1;
		edge->xr = edge->xr + edge->dxr * depth_r;
		
		if (!true)
		if (edge->xl > edge->xr)
		{
			printf("%d\n", id);
			std::cout << dy << " " << edge->dyl << " " << edge->dyr << std::endl;
			std::cout << Triangle_node->edge[0]->x << " " << Triangle_node->edge[0]->dx << " " << Triangle_node->edge[0]->dy << std::endl;
			std::cout << Triangle_node->edge[1]->x << " " << Triangle_node->edge[1]->dx << " " << Triangle_node->edge[1]->dy << std::endl;

			std::cout << Triangle_node->point[0].v << "\n" << Triangle_node->point[1].v << "\n" << Triangle_node->point[2].v << std::endl;
			system("pause");
		}
		
		edge->zl = edge->zl + edge->dzx * edge->dxl + edge->dzy;

	}
	void Print()
	{
		printf("a, b, c, d = %lf %lf %lf %lf\n", a, b, c, d);
		printf("id = %d\n", id);
		printf("dy = %d\n", dy);
		this->edge->Print();
	}
};
