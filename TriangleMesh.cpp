#include <bits/stdc++.h>
#include "TriangleMesh.h"

void TriangleMesh::LoadFile(std::string file_name)
{
	FILE *fp = fopen(file_name.c_str(), "r");
	std::ifstream ifile(file_name);
	//std::cout << fp << std::endl;
	//std::cout << file_name.c_str() << std::endl;
	std::string tag;
	tag.resize(20);
	char ch;
	while (ifile >> tag)
	{
		//std::cout << tag << std::endl;
		if (tag == std::string("v"))
		{
			Eigen::Vector3d temp;
			// fscanf(fp, "%lf %lf %lf", &temp[0], &temp[1], &temp[2]);
			ifile >> temp[0] >> temp[1] >> temp[2];
			v.push_back(temp);
		}
		else if (tag == "vn")
		{
			Eigen::Vector3d temp;
			//fscanf(fp, "%lf %lf %lf", &temp[0], &temp[1], &temp[2]);
			ifile >> temp[0] >> temp[1] >> temp[2];
			vn.push_back(temp);
		}
		else if (tag == "vt")
		{
			Eigen::Vector2d temp;
			//fscanf(fp, "%lf %lf", &temp[0], &temp[1]);
			ifile >> temp[0] >> temp[1];
			
			/*
			if (temp[0] < 0 || temp[1] < 0)
			{
				std::cout << "in load\n" << temp << std::endl;
				system("pause");
			}
			*/
			vt.push_back(temp);
		}
		else if (tag == "f")
		{
			Mesh temp;
			//fscanf(fp, "%d/%d/%d", &temp.v[0], &temp.vt[0], &temp.vn[0]);
			//fscanf(fp, "%d/%d/%d", &temp.v[1], &temp.vt[1], &temp.vn[1]);
			//fscanf(fp, "%d/%d/%d", &temp.v[2], &temp.vt[2], &temp.vn[2]);
			// temp.order[0] = 0; temp.order[1] = 1; temp.order[2] = 2;
			ifile >> temp.v[0] >> ch >> temp.vt[0] >> ch >> temp.vn[0];
			ifile >> temp.v[1] >> ch >> temp.vt[1] >> ch >> temp.vn[1];
			ifile >> temp.v[2] >> ch >> temp.vt[2] >> ch >> temp.vn[2];

			Meshes.push_back(temp);
		}
	}
}