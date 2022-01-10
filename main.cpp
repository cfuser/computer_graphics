#include <bits/stdc++.h>
#include <vector>
#include <opencv2/opencv.hpp>
#include "scan_line_z-buffer.h"
#include "TriangleMesh.h"
#include "Rasterizer.h"

#define MY_PI 3.1415926

Eigen::Matrix4d get_view_matrix(Eigen::Vector3d eye_pos)
{
	Eigen::Matrix4d view = Eigen::Matrix4d::Identity();

	Eigen::Matrix4d translate;
	translate << 1, 0, 0, -eye_pos[0],
		0, 1, 0, -eye_pos[1],
		0, 0, 1, -eye_pos[2],
		0, 0, 0, 1;

	view = translate * view;

	return view;
}

Eigen::Matrix4d get_model_matrix(float angle)
{
	Eigen::Matrix4d rotation;
	angle = angle * MY_PI / 180.f;
	rotation << cos(angle), 0, sin(angle), 0,
		0, 1, 0, 0,
		-sin(angle), 0, cos(angle), 0,
		0, 0, 0, 1;

	Eigen::Matrix4d scale;
	scale << 2.5, 0, 0, 0,
		0, 2.5, 0, 0,
		0, 0, 2.5, 0,
		0, 0, 0, 1;

	Eigen::Matrix4d translate;
	translate << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1;

	return translate * rotation * scale;
}

Eigen::Matrix4d get_projection_matrix(float eye_fov, float aspect_ratio, float zNear, float zFar)
{
	// TODO: Use the same projection matrix from the previous assignments
	Eigen::Matrix4d projection;
	Eigen::Matrix4d Persp_to_Ortho = Eigen::Matrix4d::Zero();
	//std::swap(zNear, zFar);
	//zNear *= -1; zFar *= -1;
	Persp_to_Ortho(0, 0) = zNear; Persp_to_Ortho(1, 1) = zNear;
	Persp_to_Ortho(2, 2) = zNear + zFar; Persp_to_Ortho(2, 3) = -zNear * zFar;
	Persp_to_Ortho(3, 2) = 1;

	Eigen::Matrix4d Ortho = Eigen::Matrix4d::Zero();
	Eigen::Matrix4d Ortho_Translate = Eigen::Matrix4d::Identity();
	Eigen::Matrix4d Ortho_scale = Eigen::Matrix4d::Identity();

	double t = abs(zNear) * tan(eye_fov / 180 * MY_PI / 2);
	double r = t * aspect_ratio;

	Ortho_Translate(0, 0) = 1 / r; Ortho_Translate(1, 1) = 1 / t; Ortho_Translate(2, 2) = 2 / (zNear - zFar);

	Ortho_scale(0, 3) = -0; Ortho_scale(1, 3) = -0; Ortho_scale(2, 3) = -(zNear + zFar) / 2;

	Ortho = Ortho_Translate * Ortho_scale;
	projection = Ortho * Persp_to_Ortho;

	Eigen::Matrix4d Reverse = Eigen::Matrix4d::Identity();
	Reverse(0, 0) = Reverse(1, 1) = -1;
	projection = projection * Reverse;

	return projection;
}


int main()
{
	Eigen::Vector3d tri[3];
	tri[0] << 1, 0, 0;
	tri[1] << 0, 1, 0;
	tri[2] << 0, 0, 1;

	Eigen::Vector3d test;
	test << 0.2, 0.4, 0.3;
	std::vector<double> res = computeBarycentric2D(0.2, 0.4, tri);
	std::cout << res[0] << " " << res[1] << " " << res[2] << std::endl;
	std::cout << std::numeric_limits<float>::infinity() << std::endl;
	int Height = 1024;
	int Width = 1024;
	std::vector<nodeClassifiedPolygon*> ClassifiedPolygonTable;
	ClassifiedPolygonTable.resize(Height + 1);
	std::vector<nodeClassifiedPolygon> Triangle_List;
	TriangleMesh triMesh;
	std::string file_name;
	//file_name = "D:/subject/graduate/computer graphics/example-scenes-cg21/veach-mis/veach-mis.obj";
	//file_name = "D:/subject/graduate/computer graphics/example-scenes-cg21/veach-mis/test.obj";
	file_name = "D:/subject/graduate/computer graphics/project/PA2/Z-buffer/Z-buffer/model/spot/spot_triangulated_good.obj";
	//scanf("%s", file_name);

	triMesh.LoadFile(file_name);
	std::cout << triMesh.Meshes.size() << std::endl;
	
	Rasterizer rasterizer(Width, Height);

	double angle = 140.0;
	rasterizer.set_model(get_model_matrix(angle));
	Eigen::Vector3d eye_pos = { 0, 0, 10 };
	rasterizer.set_view(get_view_matrix(eye_pos));
	rasterizer.set_projection(get_projection_matrix(45.0, 1, 0.1, 50));
	Eigen::Matrix4d mvp = rasterizer.projection * rasterizer.view * rasterizer.model;

	double f1 = (50 - 0.1) / 2.0;
	double f2 = (50 + 0.1) / 2.0;

	FILE *fp = fopen("position_transformed.txt", "w");

	cv::Mat stand(Width, Height, CV_8UC3);
	cv::Mat tex_image = cv::imread("D:/subject/graduate/computer graphics/project/PA2/Z-buffer/Z-buffer/model/spot/spot_texture.png");
	rasterizer.tex = tex_image;
	rasterizer.tex_width = tex_image.cols;
	rasterizer.tex_height = tex_image.rows;

	std::cout << tex_image.cols << " " << tex_image.rows << std::endl;
	//cv::imshow("texture", tex_image);
	//cv::waitKey(0);

	int tex_width = tex_image.cols;
	int tex_height = tex_image.rows;

	double depth_max = -DBL_MAX;
	double depth_min = DBL_MAX;
	srand(time(0));
	for (int idx = 0; idx < triMesh.Meshes.size(); idx++)
	{
		auto triangle = triMesh.Meshes[idx];
		Point point[3];
		for (int i = 0; i < 3; i++)
		{
			Eigen::Vector3d pos = triMesh.v[triangle.v[i] - 1];
			Eigen::Vector4d temp = Eigen::Vector4d(pos.x() , pos.y(), pos.z(), 1.0);

			//Eigen::Vector4d temp = mvp * Eigen::Vector4d(triMesh.v[triangle.v[i] - 1], 1.0);
			temp = mvp * temp;
			temp.x() = temp.x() / temp.w();
			temp.y() = temp.y() / temp.w();
			temp.z() = temp.z() / temp.w();

			temp.x() = 0.5 * Width * (temp.x() + 1.0);
			temp.y() = 0.5 * Height * (temp.y() + 1.0);
			temp.z() = temp.z() * f1 + f2;
			temp.z() = -temp.z();
			depth_max = std::max(depth_max, temp.z());
			depth_min = std::min(depth_min, temp.z());
			stand.at<cv::Vec3b>(Height - int(temp.y()), int(temp.x())) = cv::Vec3b(255, 0, 0);
			//std::cout << idx << " " << triangle.v[i] << " " << triangle.vn[i] << " " << triangle.vt[i] << std::endl;

			point[i].v = temp.head<3>();
			fprintf(fp, "%lf %lf %lf\n", point[i].v[0], point[i].v[1], point[i].v[2]);
			//std::cout << point[i].v << std::endl;
			
			//point[i].v = triMesh.v[triangle.v[i] -1];
			
			//point[i].v[0] = (point[i].v[0] * 10 + Height) / 2;
			//point[i].v[1] = (point[i].v[1] * 10 + Width) / 2;
			point[i].vn = triMesh.vn[triangle.vn[i] - 1];
			point[i].vt = triMesh.vt[triangle.vt[i] - 1];
		}
		
		nodeClassifiedPolygon *Polygon = new nodeClassifiedPolygon(point, idx);
		

		if (Polygon->edge.size() < 2) continue;
		for (int c = 0; c < 3; c++)
			Polygon->color[c] = rand() % 255;

		for (int c = 0; c < 3; c++)
		{
			Polygon->tex_coor[c] = triMesh.vt[triangle.vt[c] - 1];
			Eigen::Vector2d pos = Polygon->tex_coor[c];
			pos[0] = pos[0] - floor(pos[0]);
			pos[1] = pos[1] - floor(pos[1]);

			pos[0] = pos[0] * tex_width;
			pos[1] = (1 - pos[1]) * tex_height;
			//std::cout << pos << std::endl;
			//system("pause");
			auto color = tex_image.at<cv::Vec3b>(int(pos[1]), int(pos[0]));
			Polygon->color[0] = color[0];
			Polygon->color[1] = color[1];
			Polygon->color[2] = color[2];
			Polygon->tex_color[c] = { double(color[0]), double(color[1]), double(color[2]) };
		}

		//Polygon->color = Eigen::Vector3d(255.0, 0, 0);
		
		int y_max = ceil(std::fmax(std::fmax(point[0].v[1], point[1].v[1]), point[2].v[1]));
		Polygon->next = ClassifiedPolygonTable[y_max];
		ClassifiedPolygonTable[y_max] = Polygon;
		
		//ClassifiedPolygonTable[y_max].addNewPolygon(Polygon);
	}
	std::cout << "Load done" << std::endl;
	std::cout << "depth_max: " << depth_max << std::endl;
	std::cout << "depth_min: " << depth_min << std::endl;
	/*
	int tot = 0;
	for (int y = Height - 1; y >= 0; y--)
	{
		int num = 0;
		auto nodePolygon = ClassifiedPolygonTable[y];
		while (nodePolygon)
		{
			num = num + 1;
			nodePolygon = nodePolygon->next;
		}
		printf("%d %d\n", y, num);
		tot = tot + num;
	}
	printf("%d\n", tot);
	system("pause");
	return 0;
	*/
	


	int y = Height - 1;
	std::list<nodeActivedPolygon*> ActivedPolygonList;
	//std::list<nodeActivedEdge*> ActivedEdgeList;

	for (; y >= 0; y--)
	{
		// printf("%d\n", y);
		//ActivedPolygon.addNewPolygon(ClassifiedPolygonTable[y]);
		for (auto nodePolygon = ClassifiedPolygonTable[y]; nodePolygon != NULL; nodePolygon = nodePolygon->next)
			if (nodePolygon->c != 0)
			{
				ActivedPolygonList.push_back(new nodeActivedPolygon(nodePolygon, y));
				// nodeActivedPolygon *iter = *(ActivedPolygonList.end());
				// ActivedEdgeList.push_back(new nodeActivedEdge(nodePolygon->edge[0], nodePolygon->edge[1], nodePolygon, y));
				// ActivedEdgeList.push_back(new nodeActivedEdge(nodePolygon->edge[0], nodePolygon->edge[2], nodePolygon, y));
			}
		
		for (auto iter = ActivedPolygonList.begin(); iter != ActivedPolygonList.end();)
		{
			auto temp_Polygon = *iter;
			//Draw(temp_Polygon);
			//temp_Polygon->Draw(y);
			rasterizer.Draw(temp_Polygon, y);
			temp_Polygon->update();
			if (temp_Polygon->dy == 0)
				iter = ActivedPolygonList.erase(iter);
			else
				iter++;
		}
		
		//Draw();
	}
	
	/*
	for (auto iter = ActivedPolygonList.begin(); iter != ActivedPolygonList.end(); iter++)
	{
		(*iter)->Print();
	}
	*/

	//cv::Mat image(Width, Height, CV_8UC3, rasterizer.Color_buffer.data());
	cv::Mat image(Width, Height, CV_8UC3);
	
	for (int i = 0; i < Height; i++)
		for (int j = 0; j < Width; j++)
			for (int c = 0; c < 3; c++)
				image.at<cv::Vec3b>(Height - j - 1, i)[c] = rasterizer.Color_buffer[rasterizer.get_index(i, j)][c];
	
	cv::imshow("stand", stand);
	cv::imshow("image", image);
	cv::waitKey(0);
	cv::destroyAllWindows();
	system("pause");
	return 0;
}
