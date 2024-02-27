// SlabCutOpt.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "SlabCutOpt.h"
#include "tri_tri_intersect.cpp"
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <algorithm> //std::sort
#include <stdio.h>
#include <omp.h> 
#include <ctime>
int main()
{
    //model parameters
    
	int start_s = clock();

	
    
    FILE *fp_log;
    errno_t err;
    err = fopen_s(&fp_log, "results.log", "w");
    fprintf(fp_log, "Optimization results:\n");
    fclose(fp_log);
    if (readParameter("SlabCutOpt.par"))
    {
        printf("Problem during opening of the parameter file\n");

    }
    //
	write_parameters();
    initialize();
    if (read_PLY_FileList) 
    {
        readPLY_FileList();
        for (int i = 0; i < ply_file_name.size(); i++)
        {
            FILE *f_in = fopen(ply_file_name[i].c_str(), "rt");
            read_PLY_file(f_in);
            fclose(f_in);
        }

    }
    else 
    {
        FILE *f_in = fopen("test_fractures.ply", "rt");
        read_PLY_file(f_in);
        fclose(f_in);
    }
    
    
    if (read_bound)
    {
        FILE *f_bounds = fopen("bounds.dat", "rt");
        read_BOUNDS_file(f_bounds);

    }
    if (read_block_dimension)
    {
        
        read_blocks_dimensions();
    }
    else 
    {
        Block_dimensions tmp_dim;
        tmp_dim.dim_x = dim_block_x;
        tmp_dim.dim_y = dim_block_y;
        tmp_dim.dim_z = dim_block_z;
        blocks_dimensions_vector.push_back(tmp_dim);
    }
    //CreateTriangles3();
    int id_test = 0;
	if (angle_type == 1) 
	{
		//conversion from DEG to RAD:
		tetha_max = tetha_max / 180.0*pigreco;
		tetha_step = tetha_step / 180.0*pigreco;
		phi_max = phi_max / 180.0*pigreco;
		phi_step = phi_step / 180.0*pigreco;
		csi_max = csi_max / 180.0*pigreco;
		csi_step = csi_step / 180.0*pigreco;

	}

    err = fopen_s(&fp_log, "results.log", "a");
    fprintf(fp_log, "Number of dimensions to test: %d\n", blocks_dimensions_vector.size());
    for (int i = 0; i < blocks_dimensions_vector.size(); i++)
    {
        fprintf(fp_log, "[%d] dim_block_x= %f dim_block_y= %f dim_block_z= %f\n",i, blocks_dimensions_vector[i].dim_x, blocks_dimensions_vector[i].dim_y, blocks_dimensions_vector[i].dim_z);
    }
    fprintf(fp_log, "n_triangles=%d\n", mesh_to_approx.size());
    fprintf(fp_log, "\n");
    fprintf(fp_log, "N_iter i_x_zone i_y_zone i_k_zone dim_block_x dim_block_y dim_block_z Tetha Phi Csi dx dy dz n_block_inside n_block_no_intersect\n");
    fclose(fp_log);

    //////////////////////////////////////////////////////////////////////
    /////////Main cicle///////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    
	//solutions: The vector holding the whole space of solutions
	std::vector<solution_struc> solutions;
    printf("Preparing solutions...\n");
	std::vector<std::vector<std::vector<std::vector<best_results_struc>>>>optimum_solutions;
	
	optimum_solutions.resize(n_x_division);
	for (int i = 0; i < n_x_division; i++)
	{
		optimum_solutions[i].resize(n_y_division);
		for (int j = 0; j < n_y_division; j++)
		{
			optimum_solutions[i][j].resize(n_z_division);
			for (int k = 0; k < n_z_division; k++)
			{
				optimum_solutions[i][j][k].resize(blocks_dimensions_vector.size());
			}
		}
			
	}

	
    for (int i = 0; i < blocks_dimensions_vector.size(); i++) {
        dim_block_x = blocks_dimensions_vector[i].dim_x;
        dim_block_y = blocks_dimensions_vector[i].dim_y;
        dim_block_z = blocks_dimensions_vector[i].dim_z;
		for (int i_n_x_division = 0; i_n_x_division < n_x_division; i_n_x_division++)
		{
			for (int i_n_y_division = 0; i_n_y_division < n_y_division; i_n_y_division++)
			{
				for (int i_n_z_division = 0; i_n_z_division < n_z_division; i_n_z_division++)
				{
					for (double tetha = 0.0; tetha < tetha_max; tetha = tetha + tetha_step) {
						for (double phi = 0.0; phi < phi_max; phi = phi + phi_step) {
							for (double csi = 0.0; csi < csi_max; csi = csi + csi_step) {
								dx_max = dim_block_x;
								dy_max = dim_block_y;
								dz_max = dim_block_z;
								for (double dx = -dx_max / 2.0; dx < dx_max / 2.0; dx = dx + dx_step) {
									for (double dy = -dy_max / 2.0; dy < dy_max / 2.0; dy = dy + dy_step) {
										for (double dz = -dz_max / 2.0; dz < dz_max / 2.0; dz = dz + dz_step) {
											solution_struc tmp_solution;
											tmp_solution.csi = csi;
											tmp_solution.dim_block_x = dim_block_x;
											tmp_solution.dim_block_y = dim_block_y;
											tmp_solution.dim_block_z = dim_block_z;
											tmp_solution.dx = dx;
											tmp_solution.dy = dy;
											tmp_solution.dz = dz;
											tmp_solution.phi = phi;
											tmp_solution.tetha = tetha;
											tmp_solution.n_dimension = i;
											tmp_solution.xmin = xmin + i_n_x_division*(xmax - xmin) / (double)n_x_division;
											tmp_solution.xmax = xmin + (i_n_x_division+1)*(xmax - xmin) / (double)n_x_division;
											tmp_solution.ymin = ymin + i_n_y_division*(ymax - ymin) / (double)n_y_division;
											tmp_solution.ymax = ymin + (i_n_y_division + 1)*(ymax - ymin) / (double)n_y_division;
											tmp_solution.zmin = zmin + i_n_z_division*(zmax - zmin) / (double)n_z_division;
											tmp_solution.zmax = zmin + (i_n_z_division + 1)*(zmax - zmin) / (double)n_z_division;
											tmp_solution.i_domain = i_n_x_division;
											tmp_solution.j_domain = i_n_y_division;
											tmp_solution.k_domain = i_n_z_division;
											
											tmp_solution.X0 = (tmp_solution.xmax + tmp_solution.xmin) / 2.0;
											tmp_solution.Y0 = (tmp_solution.ymax + tmp_solution.ymin) / 2.0;
											tmp_solution.Z0 = (tmp_solution.zmax + tmp_solution.zmin) / 2.0;
											solutions.push_back(tmp_solution);
										}
									}
								}
							}
						}
					}
				}
			}
        }
    }
    //////////////////////////////////////////////////////////////////////
    printf("Start computation...\n");
	printf("Number of solutions is=%d\n",solutions.size());
	int count = 0;
	int reported_count = 0;
	int step_size = (int)solutions.size() / 10;

#pragma omp parallel for

    for (int i = 0; i < solutions.size(); i++) 
    {
        int n_proc = omp_get_thread_num();
		int n_threads = omp_get_num_threads();
        std::vector<std::vector<std::vector<Blocks>>>vectBlocks;
		CreateBlocks5(solutions[i], vectBlocks);
		MoveBlocks3(solutions[i], vectBlocks);
		int n_block_inside2 = CheckBlocks3(solutions[i], vectBlocks);
		int n_block_no_intersect2 = block_intersection3(solutions[i], vectBlocks);
		solutions[i].n_block_inside = n_block_inside2;
		solutions[i].n_block_no_intersect = n_block_no_intersect2;
        //printf("Calculated solution %d of %d with processor %d block_inside=%d block_no_intersect=%d\n", i, solutions.size(), n_proc, solutions[i].n_block_inside, solutions[i].n_block_no_intersect);
        char buffer[200];
        sprintf_s(buffer, "Blocks_vect%d.vtu", i);
		if (write_vtu==1)printVTU_Blocks3(buffer, solutions[i], vectBlocks);
		freemem3(solutions[i], vectBlocks);
#pragma omp atomic
		count++;

		if (omp_get_thread_num() == 0 && count >= step_size+reported_count)
		{
			double percentage = (double)count / (double)solutions.size()*100.0;
			double elaspsedtime_s = (clock() - start_s)/1000.0;
			double elaspsedtime_days = elaspsedtime_s / (24 * 60 * 60);
			
			double remaing_time_s = elaspsedtime_s / count*(double)(solutions.size()-count);
			double remaing_time_days=remaing_time_s/ (24 * 60 * 60);
			printf("Iter: %d of %d (%.2f%%) Elapsed t: %.1f(s) %.2f(days).End in: %.1f(s) %.2f(days) using %d threads\n", count, solutions.size(), percentage, elaspsedtime_s, elaspsedtime_days, remaing_time_s, remaing_time_days, n_threads);
			reported_count = count;
		}




    }
    ////
	////
	//MODIFICARE SCRITTURA PER LE SINGOLE ZONE
	////
	////
	err = fopen_s(&fp_log, "results.log", "a");
	int old_dim = -1;
	for (int i = 0; i < solutions.size(); i++)
    {
		
		if (solutions[i].n_dimension != old_dim)
		{
			fprintf(fp_log, "\nResults_for_slab_[%d] (%f %f %f)\n", solutions[i].n_dimension, solutions[i].dim_block_x, solutions[i].dim_block_y, solutions[i].dim_block_z);
			old_dim = solutions[i].n_dimension;
		}
		int i_zone = solutions[i].i_domain;
		int j_zone = solutions[i].j_domain;
		int k_zone = solutions[i].k_domain;

		if (optimum_solutions[i_zone][j_zone][k_zone][solutions[i].n_dimension].n_block_no_intersect < solutions[i].n_block_no_intersect)
		{
			optimum_solutions[i_zone][j_zone][k_zone][solutions[i].n_dimension].n_block_no_intersect = solutions[i].n_block_no_intersect;
			optimum_solutions[i_zone][j_zone][k_zone][solutions[i].n_dimension].n_block_inside = solutions[i].n_block_inside;
			optimum_solutions[i_zone][j_zone][k_zone][solutions[i].n_dimension].solution_number=i;
		}
		
		fprintf(fp_log, " %d %d %d %d %f %f %f %f %f %f %f %f %f %d %d\n", i, solutions[i].i_domain, solutions[i].j_domain, solutions[i].k_domain, solutions[i].dim_block_x, solutions[i].dim_block_y, solutions[i].dim_block_z, solutions[i].tetha, solutions[i].phi, solutions[i].csi, solutions[i].dx, solutions[i].dy, solutions[i].dz, solutions[i].n_block_inside, solutions[i].n_block_no_intersect);
    
	}
	for (int i_n_x_division = 0; i_n_x_division < n_x_division; i_n_x_division++)
	{
		for (int i_n_y_division = 0; i_n_y_division < n_y_division; i_n_y_division++)
		{
			for (int i_n_z_division = 0; i_n_z_division < n_z_division; i_n_z_division++)
			{

				fprintf(fp_log, "zone i=%d,j=%d,k=%d\n", i_n_x_division, i_n_y_division, i_n_z_division);
				for (int i = 0; i < optimum_solutions[i_n_x_division][i_n_y_division][i_n_z_division].size(); i++)
				{
					fprintf(fp_log, "Optimum solution for slab_[%d] is solution[%d],inside_blocks=%d, no_intersected_block=%d\n", i, optimum_solutions[i_n_x_division][i_n_y_division][i_n_z_division][i].solution_number, optimum_solutions[i_n_x_division][i_n_y_division][i_n_z_division][i].n_block_inside,optimum_solutions[i_n_x_division][i_n_y_division][i_n_z_division][i].n_block_no_intersect);
					if (write_vtu == 2)
					{
						int i_opt = optimum_solutions[i_n_x_division][i_n_y_division][i_n_z_division][i].solution_number;
						char buffer[200];
						sprintf_s(buffer, "Blocks_opt_i_%d_j_%d_k_%d_%d_%d.vtu", i_n_x_division,i_n_y_division,i_n_z_division, i, i_opt);
						std::vector<std::vector<std::vector<Blocks>>>vectBlocks;
						CreateBlocks5(solutions[i_opt], vectBlocks);
						MoveBlocks3(solutions[i_opt], vectBlocks);
						int n_block_inside2 = CheckBlocks3(solutions[i_opt], vectBlocks);
						int n_block_no_intersect2 = block_intersection3(solutions[i_opt], vectBlocks);
						printVTU_Blocks3(buffer, solutions[i_opt], vectBlocks);
					}
				}
			}
		}
	}
	int stop_s= clock();
	int time =( stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;
	fprintf(fp_log, "Elapsed time %d (ms)\n", time);
	fclose(fp_log);
    return 0;
}
int write_parameters() 
{
	FILE *fp_log;
	errno_t err;
	err = fopen_s(&fp_log, "results.log", "w");
	fprintf(fp_log, "*************************************************\n");
	fprintf(fp_log, "*         DICAM - University of Bologna         *\n");
	fprintf(fp_log, "*          SlabCutOpt - ver.1.0 (2018)          *\n");
	fprintf(fp_log, "*************************************************\n");
	fprintf(fp_log, "*A computer code for Slabs Cutting Optimization.*\n");
	fprintf(fp_log, "*************************************************\n");
	fprintf(fp_log, "*Authors: Stefano Bondua, Mohamed Elkarmoty      *\n");
	fprintf(fp_log, "*************************************************\n");
	fprintf(fp_log, "*For any information please contact:            *\n");
	fprintf(fp_log, "*                 stefano.bondua@unibo.it       *\n");
	fprintf(fp_log, "*************************************************\n");
	fprintf(fp_log, "*Parameters used                                *\n");
	fprintf(fp_log, "*************************************************\n");
	fprintf(fp_log, "x_max=%f\n", xmax);
	fprintf(fp_log, "x_min=%f\n", xmin);
	fprintf(fp_log, "y_max=%f\n", ymax);
	fprintf(fp_log, "y_min=%f\n", ymin);
	fprintf(fp_log, "z_max=%f\n", zmax);
	fprintf(fp_log, "z_min=%f\n", zmin);
	fprintf(fp_log, "tetha_step=%f\n", tetha_step);
	fprintf(fp_log, "tetha_max=%f\n", tetha_max);
	fprintf(fp_log, "phi_max=%f\n", phi_max);
	fprintf(fp_log, "phi_step=%f\n", phi_step);
	fprintf(fp_log, "csi_max=%f\n", csi_max);
	fprintf(fp_log, "csi_step=%f\n", csi_step);
	fprintf(fp_log, "n_x_division=%d\n", n_x_division);
	fprintf(fp_log, "n_y_division=%d\n", n_y_division);
	fprintf(fp_log, "n_z_division=%d\n", n_z_division);
	fprintf(fp_log, "dx_step=%f\n", dx_step);
	fprintf(fp_log, "dy_step=%f\n", dy_step);
	fprintf(fp_log, "dz_step=%f\n", dz_step);
	fprintf(fp_log, "dx_max=%f\n", dx_max);
	fprintf(fp_log, "dy_max=%f\n", dy_max);
	fprintf(fp_log, "dz_max=%f\n", dz_max);
	fprintf(fp_log, "dim_block_x=%f\n", dim_block_x);
	fprintf(fp_log, "dim_block_y=%f\n", dim_block_y);
	fprintf(fp_log, "dim_block_z=%f\n", dim_block_z);
	fprintf(fp_log, "read_block_dimension=%d\n", read_block_dimension);
	fprintf(fp_log, "read_bound=%d\n", read_bound);
	fprintf(fp_log, "BiDimensional=%d\n", BiDimensional);
	fprintf(fp_log, "write_vtu=%d\n", write_vtu);
	fprintf(fp_log, "read_PLY_FileList=%d\n", read_PLY_FileList);
	fprintf(fp_log, "n_x_division=%d\n", n_x_division);
	fprintf(fp_log, "n_y_division=%d\n", n_y_division);
	fprintf(fp_log, "n_z_division=%d\n", n_z_division);
	fprintf(fp_log, "rotation_method=%d\n", rotation_method);
	fprintf(fp_log, "cut_saw_thickness=%f\n", cut_saw_thickness);
	fprintf(fp_log, "angle_type=%d\n", angle_type);
	fprintf(fp_log, "bound_tolerance=%.17f\n", bound_tolerance);
	fprintf(fp_log, "end=end\n");
	fprintf(fp_log, "*************************************************\n");
	fprintf(fp_log, "*               END PARAMETERS                  *\n");
	fprintf(fp_log, "*************************************************\n");
	fclose(fp_log);
	return 0;
}
int initialize()
{
    cube_rect_faces = new int*[6];
    for (int i = 0; i < 6; i++)cube_rect_faces[i] = new int[4];
    //bottom
    cube_rect_faces[0][0] = 0;
    cube_rect_faces[0][1] = 3;
    cube_rect_faces[0][2] = 2;
    cube_rect_faces[0][3] = 1;
    //right
    cube_rect_faces[1][0] = 2;
    cube_rect_faces[1][1] = 6;
    cube_rect_faces[1][2] = 5;
    cube_rect_faces[1][3] = 1;
    //top
    cube_rect_faces[2][0] = 4;
    cube_rect_faces[2][1] = 5;
    cube_rect_faces[2][2] = 6;
    cube_rect_faces[2][3] = 7;
    //left
    cube_rect_faces[3][0] = 0;
    cube_rect_faces[3][1] = 4;
    cube_rect_faces[3][2] = 7;
    cube_rect_faces[3][3] = 3;
    //back
    cube_rect_faces[4][0] = 2;
    cube_rect_faces[4][1] = 3;
    cube_rect_faces[4][2] = 7;
    cube_rect_faces[4][3] = 6;
    //front
    cube_rect_faces[5][0] = 0;
    cube_rect_faces[5][1] = 1;
    cube_rect_faces[5][2] = 5;
    cube_rect_faces[5][3] = 4;
    
    cube_tri_faces = new int*[12];
    for (int i = 0; i < 12; i++)cube_tri_faces[i] = new int[3];
//bottom
    cube_tri_faces[0][0] = 0;
    cube_tri_faces[0][1] = 3;
    cube_tri_faces[0][2] = 1;
    
    
    cube_tri_faces[1][0] = 1;
    cube_tri_faces[1][1] = 3;
    cube_tri_faces[1][2] = 2;
    
//right	
    cube_tri_faces[2][0] = 5;
    cube_tri_faces[2][1] = 1;
    cube_tri_faces[2][2] = 2;
    
    cube_tri_faces[3][0] = 2;
    cube_tri_faces[3][1] = 6;
    cube_tri_faces[3][2] = 5;

//top
    cube_tri_faces[4][0] = 4;
    cube_tri_faces[4][1] = 5;
    cube_tri_faces[4][2] = 7;
    
    
    cube_tri_faces[5][0] = 6;
    cube_tri_faces[5][1] = 7;
    cube_tri_faces[5][2] = 5;
//back	
    cube_tri_faces[6][0] = 3;
    cube_tri_faces[6][1] = 7;
    cube_tri_faces[6][2] = 6;

    cube_tri_faces[7][0] = 6;
    cube_tri_faces[7][1] = 2;
    cube_tri_faces[7][2] = 3;
//left
    cube_tri_faces[8][0] = 4;
    cube_tri_faces[8][1] = 7;
    cube_tri_faces[8][2] = 3;

    cube_tri_faces[9][0] = 0;
    cube_tri_faces[9][1] = 4;
    cube_tri_faces[9][2] = 3;
//front
    cube_tri_faces[10][0] = 0;
    cube_tri_faces[10][1] = 5;
    cube_tri_faces[10][2] = 4;

    cube_tri_faces[11][0] = 5;
    cube_tri_faces[11][1] = 0;
    cube_tri_faces[11][2] = 1;






    FILE *fp_bounds;
    errno_t err;
    err = fopen_s(&fp_bounds, "bounds.ply", "w");
    fprintf(fp_bounds, "ply\nformat ascii 1.0\nComment Exported by voro2mesh\nelement vertex 8\nproperty float x\nproperty float y\nproperty float z\nelement face 6\nproperty list uchar int vertex_index\nend_header\n");
    fprintf(fp_bounds, "%f %f %f\n", xmin, ymin, zmin);
    fprintf(fp_bounds, "%f %f %f\n", xmax, ymin, zmin);
    fprintf(fp_bounds, "%f %f %f\n", xmax, ymax, zmin);
    fprintf(fp_bounds, "%f %f %f\n", xmin, ymax, zmin);
    fprintf(fp_bounds, "%f %f %f\n", xmin, ymin, zmax);
    fprintf(fp_bounds, "%f %f %f\n", xmax, ymin, zmax);
    fprintf(fp_bounds, "%f %f %f\n", xmax, ymax, zmax);
    fprintf(fp_bounds, "%f %f %f\n", xmin, ymax, zmax);
    for (int i = 0; i < 6; i++)
    {
        fprintf(fp_bounds, "4 ");
        for (int j = 0; j < 4; j++)fprintf(fp_bounds, "%d ", cube_rect_faces[i][j]);
        fprintf(fp_bounds, "\n");
    }
    fclose(fp_bounds);

    return 0;
}
int CreateTriangles3() {
    const int n_triangle = 2;//this only for test. We need to read somewhere triangles data.
    n_triangles = n_triangle;

    FILE *fp_mesh;
    errno_t err;
    err = fopen_s(&fp_mesh, "test_surface.ply", "w");
    
    fprintf(fp_mesh, "ply\nformat ascii 1.0\nComment Exported by voro2mesh\nelement vertex 6\nproperty float x\nproperty float y\nproperty float z\nelement face 2\nproperty list uchar int vertex_index\nend_header\n");


    mesh_to_approx.resize(n_triangle);

    
    mesh_to_approx[0].p[0] = 15.0;
    mesh_to_approx[0].p[1] =  0.0;
    mesh_to_approx[0].p[2] =  -20.0;

    mesh_to_approx[0].q[0] = 15.0;
    mesh_to_approx[0].q[1] = 10.0;
    mesh_to_approx[0].q[2] =  -20.0;

    mesh_to_approx[0].r[0] = 25.0;
    mesh_to_approx[0].r[1] = 10.0;
    mesh_to_approx[0].r[2] = 40.0;

    mesh_to_approx[1].p[0] = 25.0;
    mesh_to_approx[1].p[1] = 10.0;
    mesh_to_approx[1].p[2] = 40.0;

    mesh_to_approx[1].q[0] = 25.0;
    mesh_to_approx[1].q[1] = 0.0;
    mesh_to_approx[1].q[2] = 40.0;

    mesh_to_approx[1].r[0] = 15.0;
    mesh_to_approx[1].r[1] = 0.0;
    mesh_to_approx[1].r[2] = -20.0;

    for (int i = 0; i < 2; i++)
    {
        fprintf(fp_mesh, "%f10 %f10 %f10\n", mesh_to_approx[i].p[0], mesh_to_approx[i].p[1], mesh_to_approx[i].p[2]);
        fprintf(fp_mesh, "%f10 %f10 %f10\n", mesh_to_approx[i].q[0], mesh_to_approx[i].q[1], mesh_to_approx[i].q[2]);
        fprintf(fp_mesh, "%f10 %f10 %f10\n", mesh_to_approx[i].r[0], mesh_to_approx[i].r[1], mesh_to_approx[i].r[2]);
    }
    for (int i = 0; i < 2; i++)
    {
        fprintf(fp_mesh, "3 %d %d %d\n", i * 3 + 0, i * 3 + 1, i * 3 + 2);
    }

    fclose(fp_mesh);
    return 0;
}

//int printVTU_Blocks(char *nomefile) 
//{ 
//	//printing the output blocks configuration
//	FILE *fp;
//	errno_t err;
//	err = fopen_s(&fp, nomefile, "w");
//	int total_particles = nx*ny*nz;
//	int total_vertex = total_particles*8;
//	double xmin_l = +1e35, xmax_l = -1e35, ymin_l = +1e35, ymax_l = -1e35, zmin_l = +1e35, zmax_l = -1e35;
//	for (int i = 0; i < nx; i++)
//	{
//		for (int j = 0; j < ny; j++)
//		{
//			for (int k = 0; k < nz; k++)
//			{
//				for (int i_p = 0; i_p < 8; i_p++)
//				{
//					xmin_l = min(myBlocks[i][j][k].P1[i_p][0], xmin_l);
//					xmax_l = max(myBlocks[i][j][k].P1[i_p][0], xmax_l);
//					ymin_l = min(myBlocks[i][j][k].P1[i_p][1], ymin_l);
//					ymax_l = max(myBlocks[i][j][k].P1[i_p][1], ymax_l);
//					zmin_l = min(myBlocks[i][j][k].P1[i_p][2], zmin_l);
//					zmax_l = max(myBlocks[i][j][k].P1[i_p][2], zmax_l);
//				}
//			}
//		}
//	}
//	double rangemax_points = sqrt((xmax_l - xmin_l)*(xmax_l - xmin_l) + (ymax_l - ymin_l)*(ymax_l - ymin_l) + (zmax_l - zmin_l)*(zmax_l - zmin_l));
//	fprintf(fp, "<?xml version=\"1.0\"?>\n");
//	fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n");
//	fprintf(fp, "  <UnstructuredGrid>\n");
//	fprintf(fp, "    <Piece NumberOfPoints=\"%lu\" NumberOfCells=\"%lu\">\n", total_vertex, total_particles);
//	fprintf(fp, "      <PointData>\n");
//	fprintf(fp, "      </PointData>\n");
//	fprintf(fp,"      <CellData Scalars=\"myTYPES\">\n");
//	fprintf(fp,"        <DataArray type=\"Float32\" Name=\"RockTypes\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\">\n");
//	for (int i = 0; i < nx; i++)
//	{
//		for (int j = 0; j < ny; j++)
//		{
//			for (int k = 0; k < nz; k++)
//			{
//				
//				fprintf(fp, "        %d\n", myBlocks[i][j][k].rock_type);
//				
//			}
//		}
//	}
//	
//	fprintf(fp,"        </DataArray>\n");
//	fprintf(fp,"      </CellData>\n");
//	fprintf(fp, "      <Points>\n");
//	fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"%f\">\n", rangemax_points);
//	int n_int = 0;
//	for (int i = 0; i < nx; i++)
//	{
//		for (int j = 0; j < ny; j++)
//		{
//			for (int k = 0; k < nz; k++)
//			{
//				for (int i_p = 0; i_p < 8; i_p++)
//				{
//					fprintf(fp, "          ");
//					for (int i_c = 0; i_c < 3; i_c++)
//					{
//						fprintf(fp, "%f ", myBlocks[i][j][k].P1[i_p][i_c]);
//					}
//					fprintf(fp, "\n");
//				}
//			}
//		}
//	}
//	fprintf(fp, "\n");
//	fprintf(fp, "        </DataArray>\n");
//	fprintf(fp, "      </Points>\n");
//	fprintf(fp, "      <Cells>\n");
//	fprintf(fp, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"%lu\">\n", total_vertex);
//	fprintf(fp, "          ");
//	n_int = 0;
//	for (long int i = 0; i<total_vertex; i++)
//	{
//		fprintf(fp, "%lu ", i);
//		if (n_int == 5)
//		{
//			n_int = 0;
//			fprintf(fp, "\n          ");
//		}
//		else
//		{
//			n_int++;
//		}
//	}
//	fprintf(fp, "\n");
//	fprintf(fp, "        </DataArray>\n");
//	fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\" RangeMin=\"%d\" RangeMax=\"%d\">\n", 8, 8*total_particles);
//	fprintf(fp, "           ");
//	n_int = 0;
//	for (unsigned int i = 0; i<total_particles; i++)
//	{
//		fprintf(fp, "%d ", 8+i*8);
//		if (n_int == 5)
//		{
//			n_int = 0;
//			fprintf(fp, "\n          ");
//		}
//		else
//		{
//			n_int++;
//		}
//	}
//	fprintf(fp, "\n");
//	fprintf(fp, "        </DataArray>\n");
//	fprintf(fp, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\" RangeMin=\"42\" RangeMax=\"42\">\n");
//	fprintf(fp, "          ");
//	n_int = 0;
//	for (long int i = 0; i<total_particles; i++)
//	{
//		fprintf(fp, "42 ");//42=Polyhedron
//		if (n_int == 5)
//		{
//			n_int = 0;
//			fprintf(fp, "\n          ");
//		}
//		else
//		{
//			n_int++;
//		}
//	}
//	fprintf(fp, "\n");
//	fprintf(fp, "        </DataArray>\n");
//	int maxrange_faces = total_particles * 8-1;
//	fprintf(fp, "        <DataArray type=\"Int32\" Name=\"faces\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"%d\">\n", maxrange_faces);
//	n_int = 0;
//	
//	for (int i = 0; i < nx; i++)
//	{
//		for (int j = 0; j < ny; j++)
//		{
//			for (int k = 0; k < nz; k++)
//			{
//				fprintf(fp, "          6\n");
//				for (int i_faces = 0; i_faces < 6;i_faces++)
//				{
//					fprintf(fp, "          4");
//					for (int i_vertex = 0; i_vertex < 4;i_vertex++)
//					{
//						fprintf(fp, " %d", cube_rect_faces[i_faces][i_vertex]+n_int);
//					}
//					fprintf(fp, "\n");
//				}
//				n_int = n_int + 8;
//			}
//		}
//	}
//	/*for (unsigned int i = 0; i<faces.size(); i++)
//	{
//		fprintf(fp, "%d ", faces[i]);
//		if (n_int == 5)
//		{
//			n_int = 0;
//			fprintf(fp, "\n          ");
//		}
//		else
//		{
//			n_int++;
//		}
//	}*/
//	fprintf(fp, "\n");
//	fprintf(fp, "        </DataArray>\n");
//	int RangeMin_faceoffsets = 31;
//	int RangeMax_faceoffsets = RangeMin_faceoffsets*total_particles;
//	fprintf(fp, "        <DataArray type=\"Int32\" Name=\"faceoffsets\" format=\"ascii\" RangeMin=\"%d\" RangeMax=\"%d\">\n", RangeMin_faceoffsets, RangeMax_faceoffsets);
//	n_int = 0;
//	fprintf(fp, "          ");
//	for (unsigned int i = 0; i<total_particles; i++)
//	{
//		fprintf(fp, "%d ", (i+1)*RangeMin_faceoffsets);
//		if (n_int == 5)
//		{
//			n_int = 0;
//			fprintf(fp, "\n          ");
//		}
//		else
//		{
//			n_int++;
//		}
//	}
//	fprintf(fp, "\n");
//
//	fprintf(fp, "        </DataArray>\n");
//	fprintf(fp, "      </Cells>\n");
//	fprintf(fp, "    </Piece>\n");
//	fprintf(fp, "  </UnstructuredGrid>\n");
//	fprintf(fp, "</VTKFile>\n");
//
//	fclose(fp);
//	return 0;
//
//}
int printVTU_Blocks2(char *nomefile,solution_struc &solution, Blocks ***&blocks)
{
    //printing the output blocks configuration
    FILE *fp;
    errno_t err;
    err = fopen_s(&fp, nomefile, "w");
    int total_particles = solution.nx*solution.ny*solution.nz;
    int total_vertex = total_particles * 8;
    double xmin_l = +1e35, xmax_l = -1e35, ymin_l = +1e35, ymax_l = -1e35, zmin_l = +1e35, zmax_l = -1e35;
    for (int i = 0; i < solution.nx; i++)
    {
        for (int j = 0; j < solution.ny; j++)
        {
            for (int k = 0; k < solution.nz; k++)
            {
                for (int i_p = 0; i_p < 8; i_p++)
                {
                    xmin_l = min(blocks[i][j][k].Center[0] + solution.P1[i_p][0], xmin_l);
                    xmax_l = max(blocks[i][j][k].Center[0] + solution.P1[i_p][0], xmax_l);
                    ymin_l = min(blocks[i][j][k].Center[1] + solution.P1[i_p][1], ymin_l);
                    ymax_l = max(blocks[i][j][k].Center[1] + solution.P1[i_p][1], ymax_l);
                    zmin_l = min(blocks[i][j][k].Center[2] + solution.P1[i_p][2], zmin_l);
                    zmax_l = max(blocks[i][j][k].Center[2] + solution.P1[i_p][2], zmax_l);
                }
            }
        }
    }
    double rangemax_points = sqrt((xmax_l - xmin_l)*(xmax_l - xmin_l) + (ymax_l - ymin_l)*(ymax_l - ymin_l) + (zmax_l - zmin_l)*(zmax_l - zmin_l));
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n");
    fprintf(fp, "  <UnstructuredGrid>\n");
    fprintf(fp, "    <Piece NumberOfPoints=\"%lu\" NumberOfCells=\"%lu\">\n", total_vertex, total_particles);
    fprintf(fp, "      <PointData>\n");
    fprintf(fp, "      </PointData>\n");
    fprintf(fp, "      <CellData Scalars=\"myTYPES\">\n");
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"RockTypes\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\">\n");
    for (int i = 0; i < solution.nx; i++)
    {
        for (int j = 0; j < solution.ny; j++)
        {
            for (int k = 0; k < solution.nz; k++)
            {

                fprintf(fp, "        %d\n", blocks[i][j][k].rock_type);

            }
        }
    }

    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </CellData>\n");
    fprintf(fp, "      <Points>\n");
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"%f\">\n", rangemax_points);
    int n_int = 0;
    for (int i = 0; i < solution.nx; i++)
    {
        for (int j = 0; j < solution.ny; j++)
        {
            for (int k = 0; k < solution.nz; k++)
            {
                for (int i_p = 0; i_p < 8; i_p++)
                {
                    fprintf(fp, "          ");
                    for (int i_c = 0; i_c < 3; i_c++)
                    {
                        fprintf(fp, "%f ", blocks[i][j][k].Center[i_c] + solution.P1[i_p][i_c]);
                    }
                    fprintf(fp, "\n");
                }
            }
        }
    }
    fprintf(fp, "\n");
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </Points>\n");
    fprintf(fp, "      <Cells>\n");
    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"%lu\">\n", total_vertex);
    fprintf(fp, "          ");
    n_int = 0;
    for (long int i = 0; i<total_vertex; i++)
    {
        fprintf(fp, "%lu ", i);
        if (n_int == 5)
        {
            n_int = 0;
            fprintf(fp, "\n          ");
        }
        else
        {
            n_int++;
        }
    }
    fprintf(fp, "\n");
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\" RangeMin=\"%d\" RangeMax=\"%d\">\n", 8, 8 * total_particles);
    fprintf(fp, "           ");
    n_int = 0;
    for (unsigned int i = 0; i<total_particles; i++)
    {
        fprintf(fp, "%d ", 8 + i * 8);
        if (n_int == 5)
        {
            n_int = 0;
            fprintf(fp, "\n          ");
        }
        else
        {
            n_int++;
        }
    }
    fprintf(fp, "\n");
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\" RangeMin=\"42\" RangeMax=\"42\">\n");
    fprintf(fp, "          ");
    n_int = 0;
    for (long int i = 0; i<total_particles; i++)
    {
        fprintf(fp, "42 ");//42=Polyhedron
        if (n_int == 5)
        {
            n_int = 0;
            fprintf(fp, "\n          ");
        }
        else
        {
            n_int++;
        }
    }
    fprintf(fp, "\n");
    fprintf(fp, "        </DataArray>\n");
    int maxrange_faces = total_particles * 8 - 1;
    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"faces\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"%d\">\n", maxrange_faces);
    n_int = 0;

    for (int i = 0; i < solution.nx; i++)
    {
        for (int j = 0; j < solution.ny; j++)
        {
            for (int k = 0; k < solution.nz; k++)
            {
                fprintf(fp, "          6\n");
                for (int i_faces = 0; i_faces < 6; i_faces++)
                {
                    fprintf(fp, "          4");
                    for (int i_vertex = 0; i_vertex < 4; i_vertex++)
                    {
                        fprintf(fp, " %d", cube_rect_faces[i_faces][i_vertex] + n_int);
                    }
                    fprintf(fp, "\n");
                }
                n_int = n_int + 8;
            }
        }
    }
    
    fprintf(fp, "\n");
    fprintf(fp, "        </DataArray>\n");
    int RangeMin_faceoffsets = 31;
    int RangeMax_faceoffsets = RangeMin_faceoffsets*total_particles;
    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"faceoffsets\" format=\"ascii\" RangeMin=\"%d\" RangeMax=\"%d\">\n", RangeMin_faceoffsets, RangeMax_faceoffsets);
    n_int = 0;
    fprintf(fp, "          ");
    for (unsigned int i = 0; i<total_particles; i++)
    {
        fprintf(fp, "%d ", (i + 1)*RangeMin_faceoffsets);
        if (n_int == 5)
        {
            n_int = 0;
            fprintf(fp, "\n          ");
        }
        else
        {
            n_int++;
        }
    }
    fprintf(fp, "\n");

    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </Cells>\n");
    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "  </UnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");

    fclose(fp);
    return 0;

}
int printVTU_Blocks3(char *nomefile, solution_struc &solution, std::vector<std::vector<std::vector<Blocks>>>&blocks)
{
	//printing the output blocks configuration
	FILE *fp;
	errno_t err;
	err = fopen_s(&fp, nomefile, "w");
	int total_particles = solution.nx*solution.ny*solution.nz;
	int total_vertex = total_particles * 8;
	double xmin_l = +1e35, xmax_l = -1e35, ymin_l = +1e35, ymax_l = -1e35, zmin_l = +1e35, zmax_l = -1e35;
	for (int i = 0; i < solution.nx; i++)
	{
		for (int j = 0; j < solution.ny; j++)
		{
			for (int k = 0; k < solution.nz; k++)
			{
				for (int i_p = 0; i_p < 8; i_p++)
				{
					xmin_l = min(blocks[i][j][k].Center[0] + solution.P1[i_p][0], xmin_l);
					xmax_l = max(blocks[i][j][k].Center[0] + solution.P1[i_p][0], xmax_l);
					ymin_l = min(blocks[i][j][k].Center[1] + solution.P1[i_p][1], ymin_l);
					ymax_l = max(blocks[i][j][k].Center[1] + solution.P1[i_p][1], ymax_l);
					zmin_l = min(blocks[i][j][k].Center[2] + solution.P1[i_p][2], zmin_l);
					zmax_l = max(blocks[i][j][k].Center[2] + solution.P1[i_p][2], zmax_l);
				}
			}
		}
	}
	double rangemax_points = sqrt((xmax_l - xmin_l)*(xmax_l - xmin_l) + (ymax_l - ymin_l)*(ymax_l - ymin_l) + (zmax_l - zmin_l)*(zmax_l - zmin_l));
	fprintf(fp, "<?xml version=\"1.0\"?>\n");
	fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n");
	fprintf(fp, "  <UnstructuredGrid>\n");
	fprintf(fp, "    <Piece NumberOfPoints=\"%lu\" NumberOfCells=\"%lu\">\n", total_vertex, total_particles);
	fprintf(fp, "      <PointData>\n");
	fprintf(fp, "      </PointData>\n");
	fprintf(fp, "      <CellData Scalars=\"myTYPES\">\n");
	fprintf(fp, "        <DataArray type=\"Float32\" Name=\"RockTypes\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\">\n");
	for (int i = 0; i < solution.nx; i++)
	{
		for (int j = 0; j < solution.ny; j++)
		{
			for (int k = 0; k < solution.nz; k++)
			{

				fprintf(fp, "        %d\n", blocks[i][j][k].rock_type);

			}
		}
	}

	fprintf(fp, "        </DataArray>\n");
	fprintf(fp, "      </CellData>\n");
	fprintf(fp, "      <Points>\n");
	fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"%f\">\n", rangemax_points);
	int n_int = 0;
	for (int i = 0; i < solution.nx; i++)
	{
		for (int j = 0; j < solution.ny; j++)
		{
			for (int k = 0; k < solution.nz; k++)
			{
				for (int i_p = 0; i_p < 8; i_p++)
				{
					fprintf(fp, "          ");
					for (int i_c = 0; i_c < 3; i_c++)
					{
						fprintf(fp, "%f ", blocks[i][j][k].Center[i_c]+solution.P1[i_p][i_c]);
					}
					fprintf(fp, "\n");
				}
			}
		}
	}
	fprintf(fp, "\n");
	fprintf(fp, "        </DataArray>\n");
	fprintf(fp, "      </Points>\n");
	fprintf(fp, "      <Cells>\n");
	fprintf(fp, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"%lu\">\n", total_vertex);
	fprintf(fp, "          ");
	n_int = 0;
	for (long int i = 0; i<total_vertex; i++)
	{
		fprintf(fp, "%lu ", i);
		if (n_int == 5)
		{
			n_int = 0;
			fprintf(fp, "\n          ");
		}
		else
		{
			n_int++;
		}
	}
	fprintf(fp, "\n");
	fprintf(fp, "        </DataArray>\n");
	fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\" RangeMin=\"%d\" RangeMax=\"%d\">\n", 8, 8 * total_particles);
	fprintf(fp, "           ");
	n_int = 0;
	for (unsigned int i = 0; i<total_particles; i++)
	{
		fprintf(fp, "%d ", 8 + i * 8);
		if (n_int == 5)
		{
			n_int = 0;
			fprintf(fp, "\n          ");
		}
		else
		{
			n_int++;
		}
	}
	fprintf(fp, "\n");
	fprintf(fp, "        </DataArray>\n");
	fprintf(fp, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\" RangeMin=\"42\" RangeMax=\"42\">\n");
	fprintf(fp, "          ");
	n_int = 0;
	for (long int i = 0; i<total_particles; i++)
	{
		fprintf(fp, "42 ");//42=Polyhedron
		if (n_int == 5)
		{
			n_int = 0;
			fprintf(fp, "\n          ");
		}
		else
		{
			n_int++;
		}
	}
	fprintf(fp, "\n");
	fprintf(fp, "        </DataArray>\n");
	int maxrange_faces = total_particles * 8 - 1;
	fprintf(fp, "        <DataArray type=\"Int32\" Name=\"faces\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"%d\">\n", maxrange_faces);
	n_int = 0;

	for (int i = 0; i < solution.nx; i++)
	{
		for (int j = 0; j < solution.ny; j++)
		{
			for (int k = 0; k < solution.nz; k++)
			{
				fprintf(fp, "          6\n");
				for (int i_faces = 0; i_faces < 6; i_faces++)
				{
					fprintf(fp, "          4");
					for (int i_vertex = 0; i_vertex < 4; i_vertex++)
					{
						fprintf(fp, " %d", cube_rect_faces[i_faces][i_vertex] + n_int);
					}
					fprintf(fp, "\n");
				}
				n_int = n_int + 8;
			}
		}
	}
	
	fprintf(fp, "\n");
	fprintf(fp, "        </DataArray>\n");
	int RangeMin_faceoffsets = 31;
	int RangeMax_faceoffsets = RangeMin_faceoffsets*total_particles;
	fprintf(fp, "        <DataArray type=\"Int32\" Name=\"faceoffsets\" format=\"ascii\" RangeMin=\"%d\" RangeMax=\"%d\">\n", RangeMin_faceoffsets, RangeMax_faceoffsets);
	n_int = 0;
	fprintf(fp, "          ");
	for (unsigned int i = 0; i<total_particles; i++)
	{
		fprintf(fp, "%d ", (i + 1)*RangeMin_faceoffsets);
		if (n_int == 5)
		{
			n_int = 0;
			fprintf(fp, "\n          ");
		}
		else
		{
			n_int++;
		}
	}
	fprintf(fp, "\n");

	fprintf(fp, "        </DataArray>\n");
	fprintf(fp, "      </Cells>\n");
	fprintf(fp, "    </Piece>\n");
	fprintf(fp, "  </UnstructuredGrid>\n");
	fprintf(fp, "</VTKFile>\n");

	fclose(fp);
	return 0;

}
int freemem2(solution_struc &solution,Blocks ***&blocks)
{
    for (int i = 0; i < solution.nx; i++)
    {
        for (int j = 0; j < solution.ny; j++)
        {
            for (int k = 0; k < solution.nz; k++)
            {
//               for (int i_p = 0; i_p < 8; i_p++)
//               {
//                    delete[] blocks[i][j][k].P1[i_p];
//               }
//              delete[] blocks[i][j][k].P1;
                delete[] blocks[i][j][k].Center;
            }
            delete[] blocks[i][j];
        }
        delete[] blocks[i];
    }

    delete[] blocks;
    return 0;
}
int freemem3(solution_struc &solution, std::vector<std::vector<std::vector<Blocks>>>&blocks)
{
	for (int i = 0; i < solution.nx; i++)
	{
		for (int j = 0; j < solution.ny; j++)
		{
			for (int k = 0; k < solution.nz; k++)
			{
//				for (int i_p = 0; i_p < 8; i_p++)
//				{
//					delete[] blocks[i][j][k].P1[i_p];
//				}
//				delete[] blocks[i][j][k].P1;
				delete[] blocks[i][j][k].Center;
			}
			
		}
		
	}

	
	return 0;
}
//int block_intersection2(solution_struc& solution,Blocks ***&blocks)
//{
//    int n_blocks_no_intersected = 0;
//
//    for (int i = 0; i < solution.nx; i++)
//    {
//        for (int j = 0; j < solution.ny; j++)
//        {
//            for (int k = 0; k < solution.nz; k++)
//            {
//                blocks[i][j][k].intersected = false;
//                if (blocks[i][j][k].inside)
//                {
//                    for (int i_t = 0; i_t < 12; i_t++)
//                    {
//                        for (int i_mesh_t = 0; i_mesh_t < mesh_to_approx.size(); i_mesh_t++)
//                        {
//                            int test_tri_tri = tri_tri_overlap_test_3d(blocks[i][j][k].P1[cube_tri_faces[i_t][0]], blocks[i][j][k].P1[cube_tri_faces[i_t][1]], blocks[i][j][k].P1[cube_tri_faces[i_t][2]], mesh_to_approx[i_mesh_t].p, mesh_to_approx[i_mesh_t].q, mesh_to_approx[i_mesh_t].r);
//                            if (test_tri_tri)
//                            {
//                                blocks[i][j][k].intersected = true;
//                                blocks[i][j][k].rock_type = 2;
//                                break;
//                            }
//                        }
//                    }
//                    if (!blocks[i][j][k].intersected)n_blocks_no_intersected++;
//                }
//
//            }
//        }
//    }
//    return n_blocks_no_intersected;
//
//
//}

int block_intersection3(solution_struc& solution, std::vector<std::vector<std::vector<Blocks>>>&blocks)
{
	int n_blocks_no_intersected = 0;

	for (int i = 0; i < solution.nx; i++)
	{
		for (int j = 0; j < solution.ny; j++)
		{
			for (int k = 0; k < solution.nz; k++)
			{
				blocks[i][j][k].intersected = false;
				if (blocks[i][j][k].inside)
				{
					for (int i_t = 0; i_t < 12; i_t++)
					{
						for (int i_mesh_t = 0; i_mesh_t < mesh_to_approx.size(); i_mesh_t++)
						{
							
							double *p1 = new double[3];
							double *p2 = new double[3];
							double *p3 = new double[3];
							
							for (int ixyz = 0; ixyz < 3; ixyz++)
							{
								p1[ixyz] = blocks[i][j][k].Center[ixyz] + solution.P1[cube_tri_faces[i_t][0]][ixyz];
								p2[ixyz] = blocks[i][j][k].Center[ixyz] + solution.P1[cube_tri_faces[i_t][1]][ixyz];
								p3[ixyz] = blocks[i][j][k].Center[ixyz] + solution.P1[cube_tri_faces[i_t][2]][ixyz];
							}
							
							int test_tri_tri = tri_tri_overlap_test_3d(p1,p2,p3, mesh_to_approx[i_mesh_t].p, mesh_to_approx[i_mesh_t].q, mesh_to_approx[i_mesh_t].r);
							free(p1);
							free(p2);
							free(p3);
							if (test_tri_tri)
							{
								blocks[i][j][k].intersected = true;
								blocks[i][j][k].rock_type = 2;
								break;
							}
						}
					}
					if (!blocks[i][j][k].intersected)n_blocks_no_intersected++;
				}

			}
		}
	}
	return n_blocks_no_intersected;


}


bool isInside(double x, double y, double z, solution_struc &solution)
{
    if (x<solution.xmin - bound_tolerance || x>solution.xmax + bound_tolerance)return false;
    if (y<solution.ymin - bound_tolerance || y>solution.ymax + bound_tolerance)return false;
    if (z<solution.zmin - bound_tolerance || z>solution.zmax + bound_tolerance)return false;
    return true;
}
bool isInside(double *&P, solution_struc &solution)
{
    if (P[0]<solution.xmin - bound_tolerance || P[0]>solution.xmax + bound_tolerance)return false;
    if (P[1]<solution.ymin - bound_tolerance || P[1]>solution.ymax + bound_tolerance)return false;
    if (P[2]<solution.zmin - bound_tolerance || P[2]>solution.zmax + bound_tolerance)return false;
    return true;
}

int rotate_point(double tetha, double phi, double csi, double *&point) 
{
    double x1 = point[0];
    double y1 = point[1];
    double z1 = point[2];
    
    double e1a = cos(phi)*cos(csi) - cos(tetha)*sin(phi)*sin(csi);
    double e1b = cos(phi)*sin(csi) + cos(tetha)*sin(phi)*cos(csi);
    double e1c = sin(tetha)*sin(phi);
    double e2a = -(sin(phi)*cos(csi) + cos(tetha)*cos(phi)*sin(csi));
    double e2b = -(sin(phi)*sin(csi) - cos(tetha)*cos(phi)*cos(csi));
    double e2c = sin(tetha)*cos(phi);
    double e3a = sin(tetha)*sin(csi);
    double e3b = -sin(tetha)*cos(csi);
    double e3c = cos(tetha);
    point[0] = e1a*x1 + e1b*y1 + e1c*z1;
    point[1] = e2a*x1 + e2b*y1 + e2c*z1;
    point[2] = e3a*x1 + e3b*y1 + e3c*z1;
    return 0;
}

int rotate_point2(double tetha, double phi, double csi, double *&point)
{
    double x1 = point[0];
    double y1 = point[1];
    double z1 = point[2];
    
    //rotation around y
    double x2 = x1*cos(tetha) - z1*sin(tetha);
    double y2 = y1;
    double z2 = x1*sin(tetha) + z1*cos(tetha);


    //rotation around z
    double x3 = x2*cos(phi) - y2*sin(phi);
    double y3 = x2*sin(phi) + y2*cos(phi);
    double z3 = z2;
    
    //rotation around x
    double x4 = x3;
    double y4 = y3*cos(csi) - z3*sin(csi);
    double z4 = y3*sin(csi) + z3*cos(csi);



    point[0] = x4;
    point[1] = y4;
    point[2] = z4;
    return 0;
}


int MoveBlocks3(solution_struc &solution, std::vector<std::vector<std::vector<Blocks>>>&blocks)
{
	//rotation:the center of rotation is X0 Y0 Z0


	for (int i_p = 0; i_p < 8; i_p++)
	{
		if (rotation_method == 1)
		{
			rotate_point(solution.tetha, solution.phi, solution.csi, solution.P1[i_p]);
		}
		else
		{
			rotate_point2(solution.tetha, solution.phi, solution.csi, solution.P1[i_p]);
		}

		/*blocks[i][j][k].P1[i_p][0] += solution.X0 + solution.dx;
		blocks[i][j][k].P1[i_p][1] += solution.Y0 + solution.dy;
		blocks[i][j][k].P1[i_p][2] += solution.Z0 + solution.dz;*/
	}

	for (int i = 0; i < solution.nx; i++)
	{
		for (int j = 0; j < solution.ny; j++)
		{
			for (int k = 0; k < solution.nz; k++)
			{

				if (rotation_method == 1)
				{
					rotate_point(solution.tetha, solution.phi, solution.csi, blocks[i][j][k].Center);
				}
				else
				{
					rotate_point2(solution.tetha, solution.phi, solution.csi, blocks[i][j][k].Center);
				}

				blocks[i][j][k].Center[0] += solution.X0 + solution.dx;
				blocks[i][j][k].Center[1] += solution.Y0 + solution.dy;
				blocks[i][j][k].Center[2] += solution.Z0 + solution.dz;
				
			}
		}
	}
	return 0;
}



int pnpoly(int nvert, double *vertx, double *verty, double testx, double testy)
{
    int i, j, c = 0;
    for (i = 0, j = nvert - 1; i < nvert; j = i++) {
        if (((verty[i]>testy) != (verty[j]>testy)) &&
            (testx < (vertx[j] - vertx[i]) * (testy - verty[i]) / (verty[j] - verty[i]) + vertx[i]))
            c = !c;
    }
    return c;
}

//int CreateBlocks4(solution_struc &solution, std::vector<std::vector<std::vector<Blocks>>>&blocks)
//{
//	//da modificare: ogni volta ne dobbiamo usare un numero diverso?
//	int nx=1;
//	int ny=1;
//	int nz=1;
//	double **bounds = new double*[8];
//	for (int i = 0; i < 8; i++) bounds[i] = new double[3];
//	bounds[0][0] = 0;
//	bounds[0][1] = 0;
//	bounds[0][2] = 0;
//
//	bounds[1][0] = solution.dim_block_x;
//	bounds[1][1] = 0;
//	bounds[1][2] = 0;
//
//	bounds[2][0] = solution.dim_block_x;
//	bounds[2][1] = solution.dim_block_y;
//	bounds[2][2] = 0;
//
//	bounds[3][0] = 0;
//	bounds[3][1] = solution.dim_block_y;
//	bounds[3][2] = 0;
//	for (int i = 0; i < 4; i++)
//	{
//		bounds[i + 4][0] = bounds[i][0];
//		bounds[i + 4][1] = bounds[i][1];
//		bounds[i + 4][2] = solution.dim_block_z;
//	}
//	double pigrecomezzi = atan(1.0) * 2.0;
//	for (int i = 0; i < 8; i++)
//	{
//		if (rotation_method == 1)
//		{
//			rotate_point(solution.tetha, solution.phi, solution.csi, bounds[i]);
//		}
//
//		else
//		{
//			rotate_point2(solution.tetha, solution.phi, solution.csi, bounds[i]);
//		}
//
//	}
//
//
//	int index_x = 0;
//	int index_y = 1;
//	int index_z = 2;
//	double dim_x_tmp = abs(bounds[1][0]);
//	double dim_y_tmp = abs(bounds[3][1]);
//	double dim_z_tmp = abs(bounds[4][2]);
//	if (abs(bounds[3][0])>dim_x_tmp)
//	{
//		dim_x_tmp = abs(bounds[3][0]);
//		index_x = 1;
//	}
//	if (abs(bounds[4][0])>dim_x_tmp)
//	{
//		dim_x_tmp = abs(bounds[4][0]);
//		index_x = 2;
//	}
//	/////////////////////////////////////////
//	if (abs(bounds[1][1])>dim_y_tmp)
//	{
//		dim_y_tmp = abs(bounds[1][1]);
//		index_y = 0;
//	}
//	if (abs(bounds[4][1])>dim_y_tmp)
//	{
//		dim_y_tmp = abs(bounds[4][1]);
//		index_y = 2;
//	}
//	/////////////////////////////////////////////7
//	if (abs(bounds[1][2])>dim_z_tmp)
//	{
//		dim_z_tmp = abs(bounds[1][2]);
//		index_z = 0;
//	}
//	if (abs(bounds[3][2])>dim_z_tmp)
//	{
//		dim_z_tmp = abs(bounds[3][2]);
//		index_z = 1;
//	}
//
//	int check = index_x + index_y + index_z;
//	if (check != 3) 
//	{
//		int some_error_here = 1;
//	}
//
//
//	if (index_x == 0)
//	{
//		nx = (solution.xmax - solution.xmin + dx_max) / dim_x_tmp + 2;
//	}
//	else if (index_x == 1)
//	{
//		ny = (solution.xmax - solution.xmin + dx_max) / dim_x_tmp + 2;
//	}
//	else
//	{
//		nz = (solution.xmax - solution.xmin + dx_max) / dim_x_tmp + 2;
//	}
//	//
//	if (index_y == 0)
//	{
//		nx = (solution.ymax - solution.ymin + dy_max) / dim_y_tmp + 2;
//	}
//	else if (index_y == 1)
//	{
//		ny = (solution.ymax - solution.ymin + dy_max) / dim_y_tmp + 2;
//	}
//	else
//	{
//		nz = (solution.ymax - solution.ymin + dy_max) / dim_y_tmp + 2;
//	}
//	//
//	if (index_z == 0)
//	{
//		nx = (solution.xmax - solution.xmin + dx_max) / dim_x_tmp + 2;
//	}
//	else if (index_z == 1)
//	{
//		ny = (solution.ymax - solution.ymin + dy_max) / dim_y_tmp + 2;
//	}
//	else
//	{
//		nz = (solution.zmax - solution.zmin + dz_max) / dim_z_tmp + 2;
//	}
//	/*nx = (solution.xmax - solution.xmin + dx_max) / dim_x_tmp + 2;
//	ny = (solution.ymax - solution.ymin + dy_max) / dim_y_tmp + 2;
//	nz = (solution.zmax - solution.zmin + dz_max) / dim_z_tmp + 2;*/
//
//	if (nx - ((int)(nx / 2) * 2)) nx++;
//	if (ny - ((int)(ny / 2) * 2)) ny++;
//	if (nz - ((int)(nz / 2) * 2)) nz++;
//
//	if (BiDimensional) nz = 1;
//	
//	solution.nx = nx;
//	solution.ny = ny;
//	solution.nz = nz;
//
//	
//	blocks.resize(nx);
//	for (int i = 0; i < nx; i++)
//	{
//		blocks[i].resize(ny);
//		for (int j = 0; j < ny; j++)
//			blocks[i][j].resize(nz);
//	}
//
//
//
//	
//	for (int x = 0; x < nx; ++x) {
//		
//		for (int y = 0; y < ny; ++y) {
//			
//			for (int z = 0; z < nz; ++z) 
//			{ // initialize the values to whatever you want the default to be
//				blocks[x][y][z].Center = new double[3];
//				blocks[x][y][z].Center[0] = +(x - nx / 2)*solution.dim_block_x + solution.dim_block_x / 2.0;
//				blocks[x][y][z].Center[1] = +(y - ny / 2)*solution.dim_block_y + solution.dim_block_y / 2.0;
//				blocks[x][y][z].Center[2] = +(z - nz / 2)*solution.dim_block_z + solution.dim_block_z / 2.0;
//				blocks[x][y][z].P1 = new double *[8];
//				for (int i_p = 0; i_p < 8; ++i_p)
//				{
//					blocks[x][y][z].P1[i_p] = new double[3];
//				}
//				//first point
//				blocks[x][y][z].P1[0][0] = blocks[x][y][z].Center[0] - solution.dim_block_x / 2.0 + alf_cut_saw_thickness;
//				blocks[x][y][z].P1[0][1] = blocks[x][y][z].Center[1] - solution.dim_block_y / 2.0 + alf_cut_saw_thickness;
//				blocks[x][y][z].P1[0][2] = blocks[x][y][z].Center[2] - solution.dim_block_z / 2.0 + alf_cut_saw_thickness;
//				//2nd point
//				blocks[x][y][z].P1[1][0] = blocks[x][y][z].Center[0] + solution.dim_block_x / 2.0 - alf_cut_saw_thickness;
//				blocks[x][y][z].P1[1][1] = blocks[x][y][z].Center[1] - solution.dim_block_y / 2.0 + alf_cut_saw_thickness;
//				blocks[x][y][z].P1[1][2] = blocks[x][y][z].Center[2] - solution.dim_block_z / 2.0 + alf_cut_saw_thickness;
//				//3nd point
//				blocks[x][y][z].P1[2][0] = blocks[x][y][z].Center[0] + solution.dim_block_x / 2.0 - alf_cut_saw_thickness;
//				blocks[x][y][z].P1[2][1] = blocks[x][y][z].Center[1] + solution.dim_block_y / 2.0 - alf_cut_saw_thickness;
//				blocks[x][y][z].P1[2][2] = blocks[x][y][z].Center[2] - solution.dim_block_z / 2.0 + alf_cut_saw_thickness;
//				//4nd point
//				blocks[x][y][z].P1[3][0] = blocks[x][y][z].Center[0] - solution.dim_block_x / 2.0 + alf_cut_saw_thickness;
//				blocks[x][y][z].P1[3][1] = blocks[x][y][z].Center[1] + solution.dim_block_y / 2.0 - alf_cut_saw_thickness;
//				blocks[x][y][z].P1[3][2] = blocks[x][y][z].Center[2] - solution.dim_block_z / 2.0 + alf_cut_saw_thickness;
//				//5nd point
//				blocks[x][y][z].P1[4][0] = blocks[x][y][z].Center[0] - solution.dim_block_x / 2.0 + alf_cut_saw_thickness;
//				blocks[x][y][z].P1[4][1] = blocks[x][y][z].Center[1] - solution.dim_block_y / 2.0 + alf_cut_saw_thickness;
//				blocks[x][y][z].P1[4][2] = blocks[x][y][z].Center[2] + solution.dim_block_z / 2.0 - alf_cut_saw_thickness;
//				//6nd point
//				blocks[x][y][z].P1[5][0] = blocks[x][y][z].Center[0] + solution.dim_block_x / 2.0 - alf_cut_saw_thickness;
//				blocks[x][y][z].P1[5][1] = blocks[x][y][z].Center[1] - solution.dim_block_y / 2.0 + alf_cut_saw_thickness;
//				blocks[x][y][z].P1[5][2] = blocks[x][y][z].Center[2] + solution.dim_block_z / 2.0 - alf_cut_saw_thickness;
//				//7nd point
//				blocks[x][y][z].P1[6][0] = blocks[x][y][z].Center[0] + solution.dim_block_x / 2.0 - alf_cut_saw_thickness;
//				blocks[x][y][z].P1[6][1] = blocks[x][y][z].Center[1] + solution.dim_block_y / 2.0 - alf_cut_saw_thickness;
//				blocks[x][y][z].P1[6][2] = blocks[x][y][z].Center[2] + solution.dim_block_z / 2.0 - alf_cut_saw_thickness;
//				//8nd point
//				blocks[x][y][z].P1[7][0] = blocks[x][y][z].Center[0] - solution.dim_block_x / 2.0 + alf_cut_saw_thickness;
//				blocks[x][y][z].P1[7][1] = blocks[x][y][z].Center[1] + solution.dim_block_y / 2.0 - alf_cut_saw_thickness;
//				blocks[x][y][z].P1[7][2] = blocks[x][y][z].Center[2] + solution.dim_block_z / 2.0 - alf_cut_saw_thickness;
//
//
//			}
//		}
//	}
//	return 0;
//}
int CreateBlocks5(solution_struc &solution, std::vector<std::vector<std::vector<Blocks>>>&blocks)
{
	//da modificare: ogni volta ne dobbiamo usare un numero diverso?
	int nx = 1;
	int ny = 1;
	int nz = 1;
	

	double diagonal = sqrt(pow((solution.xmax - solution.xmin + dx_max), 2) + pow(solution.ymax - solution.ymin + dy_max, 2) + pow(solution.zmax - solution.zmin + dz_max , 2));
	


	nx = (diagonal) / solution.dim_block_x + 2;
	ny = (diagonal) / solution.dim_block_y + 2;
	nz = (diagonal) / solution.dim_block_z + 2;
	



	//we want not even..
	if (nx - ((int)(nx / 2) * 2)) nx++;
	if (ny - ((int)(ny / 2) * 2)) ny++;
	if (nz - ((int)(nz / 2) * 2)) nz++;

	if (BiDimensional) nz = 1;

	solution.nx = nx;
	solution.ny = ny;
	solution.nz = nz;
	blocks.resize(nx);
	for (int i = 0; i < nx; i++)
	{
		blocks[i].resize(ny);
		for (int j = 0; j < ny; j++)
			blocks[i][j].resize(nz);
	}
	solution.P1 = new double *[8];
	for (int i_p = 0; i_p < 8; ++i_p)
	{
		solution.P1[i_p] = new double[3];
	}
	//first point
	solution.P1[0][0] =  - solution.dim_block_x / 2.0 + alf_cut_saw_thickness;
	solution.P1[0][1] =  - solution.dim_block_y / 2.0 + alf_cut_saw_thickness;
	solution.P1[0][2] =  - solution.dim_block_z / 2.0 + alf_cut_saw_thickness;
	//2nd point
	solution.P1[1][0] =  + solution.dim_block_x / 2.0 - alf_cut_saw_thickness;
	solution.P1[1][1] =  - solution.dim_block_y / 2.0 + alf_cut_saw_thickness;
	solution.P1[1][2] =  - solution.dim_block_z / 2.0 + alf_cut_saw_thickness;
	//3nd point
	solution.P1[2][0] =  + solution.dim_block_x / 2.0 - alf_cut_saw_thickness;
	solution.P1[2][1] =  + solution.dim_block_y / 2.0 - alf_cut_saw_thickness;
	solution.P1[2][2] =  - solution.dim_block_z / 2.0 + alf_cut_saw_thickness;
	//4nd point
	solution.P1[3][0] =  - solution.dim_block_x / 2.0 + alf_cut_saw_thickness;
	solution.P1[3][1] =  + solution.dim_block_y / 2.0 - alf_cut_saw_thickness;
	solution.P1[3][2] =  - solution.dim_block_z / 2.0 + alf_cut_saw_thickness;
	//5nd point
	solution.P1[4][0] =  - solution.dim_block_x / 2.0 + alf_cut_saw_thickness;
	solution.P1[4][1] =  - solution.dim_block_y / 2.0 + alf_cut_saw_thickness;
	solution.P1[4][2] =  + solution.dim_block_z / 2.0 - alf_cut_saw_thickness;
	//6nd point
	solution.P1[5][0] =  + solution.dim_block_x / 2.0 - alf_cut_saw_thickness;
	solution.P1[5][1] =  - solution.dim_block_y / 2.0 + alf_cut_saw_thickness;
	solution.P1[5][2] =  + solution.dim_block_z / 2.0 - alf_cut_saw_thickness;
	//7nd point
	solution.P1[6][0] =  + solution.dim_block_x / 2.0 - alf_cut_saw_thickness;
	solution.P1[6][1] =  + solution.dim_block_y / 2.0 - alf_cut_saw_thickness;
	solution.P1[6][2] =  + solution.dim_block_z / 2.0 - alf_cut_saw_thickness;
	//8nd point
	solution.P1[7][0] =  - solution.dim_block_x / 2.0 + alf_cut_saw_thickness;
	solution.P1[7][1] =  + solution.dim_block_y / 2.0 - alf_cut_saw_thickness;
	solution.P1[7][2] =  + solution.dim_block_z / 2.0 - alf_cut_saw_thickness;



	for (int x = 0; x < nx; ++x) {

		for (int y = 0; y < ny; ++y) {

			for (int z = 0; z < nz; ++z)
			{ // initialize the values to whatever you want the default to be
				blocks[x][y][z].Center = new double[3];
				blocks[x][y][z].Center[0] = +(x - nx / 2)*solution.dim_block_x + solution.dim_block_x / 2.0;
				blocks[x][y][z].Center[1] = +(y - ny / 2)*solution.dim_block_y + solution.dim_block_y / 2.0;
				blocks[x][y][z].Center[2] = +(z - nz / 2)*solution.dim_block_z + solution.dim_block_z / 2.0;
				


			}
		}
	}
	return 0;
}
int CheckBlocks3(solution_struc&solution, std::vector<std::vector<std::vector<Blocks>>>&blocks)
{
	int n_blocks_inside = 0;
	for (int i = 0; i < solution.nx; i++)
	{
		for (int j = 0; j < solution.ny; j++)
		{
			for (int k = 0; k < solution.nz; k++)
			{
				
				blocks[i][j][k].inside = true;
				blocks[i][j][k].rock_type = 1;
				for (int i_p = 0; i_p < 8; i_p++)
				{
					
					
					if (!isInside(blocks[i][j][k].Center[0]+ solution.P1[i_p][0], blocks[i][j][k].Center[1] + solution.P1[i_p][1], blocks[i][j][k].Center[2] + solution.P1[i_p][2],solution))
					{
						blocks[i][j][k].inside = false;
						blocks[i][j][k].rock_type = -2;
						break;
					}
				}
				//
				if (read_bound)
				{
					for (int i_p = 0; i_p < 8; i_p++)
					{
						if (blocks[i][j][k].inside)
						{
							double x = blocks[i][j][k].Center[0] + solution.P1[i_p][0];
							double y = blocks[i][j][k].Center[1] + solution.P1[i_p][1];
							if (!pnpoly(n_verts_bound, verts_bound_x, verts_bound_y, x, y))
							{
								blocks[i][j][k].inside = false;
								blocks[i][j][k].rock_type = -1;
								break;
							}
						}

					}
				}
				
				if (blocks[i][j][k].inside)
					n_blocks_inside++;
			}
		}
	}
	return n_blocks_inside;
}

double min(double x, double y)
{
    if (x<y)return x;
    return y;
}
double max(double x, double y)
{
    if (x>y)return x;
    return y;
}

double distance(double x1, double y1, double z1, double x2, double y2, double z2) {
    double dist = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2));
    return dist;
}
#define SMALL_NUM   0.00000001 // anything that avoids division overflow
// dot product (3D) which allows vector operations in arguments
#define dot(u,v)   ((u).x * (v).x + (u).y * (v).y + (u).z * (v).z)

// intersect3D_RayTriangle(): find the 3D intersection of a ray with a triangle
//    Input:  a ray R, and a triangle T
//    Output: *I = intersection point (when it exists)
//    Return: -1 = triangle is degenerate (a segment or point)
//             0 =  disjoint (no intersect)
//             1 =  intersect in unique point I1
//             2 =  are in the same plane

int intersect3D_RayTriangle(Ray R, Triangle T, Point* I, Vector *N)
{
    Vector    u, v, n;              // triangle vectors
    Vector    dir, w0, w, d1;           // ray vectors
    double     r, a, b;              // params to calc ray-plane intersect

                                     // get triangle edge vectors and plane normal
                                     /*u = T.V1 - T.V0;
                                     v = T.V2 - T.V0;*/
    u = Difference(T.V1, T.V0);
    v = Difference(T.V2, T.V0);
    // n = u * v;              // cross product
    n = crossP(u, v);
    if (n.x == 0.0&& n.y == 0.0&&n.z == 0.0)             // triangle is degenerate
        return -1;                  // do not deal with this case
    (*N) = NormalizeV(n);
    //N->x=n.x;
    //N->y=n.y;
    //N->z=n.z;

    //dir = R.P1 - R.P0;              // ray direction vector
    dir = Difference(R.P1, R.P0);
    //w0 = R.P0 - T.V0;
    w0 = Difference(R.P0, T.V0);
    a = -dot(n, w0);
    b = dot(n, dir);
    if (fabs(b) < SMALL_NUM) {     // ray is  parallel to triangle plane
        if (a == 0)                 // ray lies in triangle plane
            return 2;
        else return 0;              // ray disjoint from plane
    }

    // get intersect point of ray with triangle plane
    r = a / b;
    if (r < 0.0)                    // ray goes away from triangle
        return 0;                   // => no intersect
                                    // for a segment, also test if (r > 1.0) => no intersect

                                    //*I = R.P0 + r * dir;            // intersect point of ray and plane
    Point temp = SumPoint(R.P0, ScalarVector(r, dir));
    I->x = temp.x;
    I->y = temp.y;
    I->z = temp.z;
    // is I inside T?
    double    uu, uv, vv, wu, wv, D;
    uu = dot(u, u);
    uv = dot(u, v);
    vv = dot(v, v);
    //w = *I - T.V0;
    w = Difference(*I, T.V0);
    wu = dot(w, u);
    wv = dot(w, v);
    D = uv * uv - uu * vv;

    // get and test parametric coords
    double s, t;
    s = (uv * wv - vv * wu) / D;
    if (s < 0.0 || s > 1.0)         // I is outside T
        return 0;
    t = (uv * wu - uu * wv) / D;
    if (t < 0.0 || (s + t) > 1.0)  // I is outside T
        return 0;
    //warning: the ray is in the triangle, but is in the segme?
    d1 = Difference(*I, R.P0);
    double s1;
    double testx = dir.x;
    double testy = dir.y;
    double testz = dir.z;
    if (testx != 0.0) {
        s1 = d1.x / testx;
    }
    else if (dir.y != 0.0)
    {
        s1 = d1.y / testy;
    }
    else
    {
        s1 = d1.z / testz;
    }
    if (s1 < 0.0 || s1 > 1.0)         // I is outside R
        return 0;
    return 1;                       // I is in T
}
Vector Difference(Point v1, Point v2) {
    Vector v3;
    v3.x = v1.x - v2.x;
    v3.y = v1.y - v2.y;
    v3.z = v1.z - v2.z;
    return v3;
}
Vector Difference(Vector v1, Vector v2) {
    Vector v3;
    v3.x = v1.x - v2.x;
    v3.y = v1.y - v2.y;
    v3.z = v1.z - v2.z;
    return v3;
}

//Vector Difference(Point* v1, Point v2){
//     Vector v3;
//     v3.x=v1.x-v2.x;
//     v3.y=v1.y-v2.y;
//     v3.z=v1.z-v2.z;
//     return v3;
//}
Vector Sum(Point v1, Point v2) {
    Vector v3;
    v3.x = v1.x + v2.x;
    v3.y = v1.y + v2.y;
    v3.z = v1.z + v2.z;
    return v3;
}
Point SumPoint(Point v1, Point v2) {
    Point v3;
    v3.x = v1.x + v2.x;
    v3.y = v1.y + v2.y;
    v3.z = v1.z + v2.z;
    return v3;
}
Point SumPoint(Vector v1, Vector v2) {
    Point v3;
    v3.x = v1.x + v2.x;
    v3.y = v1.y + v2.y;
    v3.z = v1.z + v2.z;
    return v3;
}
Vector SumVector(Vector v1, Vector v2) {
    Vector v3;
    v3.x = v1.x + v2.x;
    v3.y = v1.y + v2.y;
    v3.z = v1.z + v2.z;
    return v3;
}
Point SumPoint(Point v1, Vector v2) {
    Point v3;
    v3.x = v1.x + v2.x;
    v3.y = v1.y + v2.y;
    v3.z = v1.z + v2.z;
    return v3;
}
Vector NormalizeV(Vector V1) {
    Vector v3;

    double normp = norm(V1);
    v3.x = V1.x / normp;
    v3.y = V1.y / normp;
    v3.z = V1.z / normp;
    return v3;
}
Vector ScalarVector(double r, Point v1) {
    Vector v3;
    v3.x = v1.x*r;
    v3.y = v1.y*r;
    v3.z = v1.z*r;
    return v3;
}
double norm(Vector v)

{
    double d;
    d = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);


    return d;
}
Vector crossP(Vector v1, Vector v2) {


    Vector v3;
    v3.x = v1.y*v2.z - v2.y*v1.z;
    v3.y = v2.x*v1.z - v1.x*v2.z;
    v3.z = v1.x*v2.y - v1.y*v2.x;
    return v3;
}
double scalarP(Vector v1, Vector v2) {
    double r = v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
    return r;
}
Vector ScalarVector(double r, Vector v1) {
    Vector v3;
    v3.x = v1.x*r;
    v3.y = v1.y*r;
    v3.z = v1.z*r;
    return v3;
}


int readParameter(const char* Nomefilein) {
    //18/08/2013
    //Reads data from an ASCII file.
    //character "!" used as comment symbol, all character following "!" are ignored.
    //data have a name and a value, e.g.: "Max=100".Symbol "=" is used to separate data name and data value.
    FILE *f_in;
    int EndFile;
    std::string head, value;
    //char str[1000];
    errno_t err;

    err=fopen_s(&f_in,Nomefilein, "rt");
    if (err) {
        printf("\n\n\n\t\t\t\t OPEN INPUT %s FILE ERROR", Nomefilein);
        return 1;
    }
    while (true) {
        EndFile = fgetstringline(lines, f_in);

        lines = EliminateCommentAndExtract(lines, head, value);
        if (lines.empty()) { continue; }
        if (head == "x_max")
        {
            xmax = atof(value.c_str());
        }
        else if (head == "x_min")
        {
            xmin = atof(value.c_str());
        }
        else if (head == "y_max")
        {
            ymax = atof(value.c_str());
        }
        else if (head == "y_min")
        {
            ymin = atof(value.c_str());
        }
        else if (head == "z_max")
        {
            zmax = atof(value.c_str());
        }
        else if (head == "z_min")
        {
            zmin = atof(value.c_str());
        }
        
        else if (head == "tetha_step")
        {
            tetha_step = atof(value.c_str());
        }
        else if (head == "tetha_max")
        {
            tetha_max = atof(value.c_str());
        }
        else if (head == "phi_max")
        {
            phi_max = atof(value.c_str());
        }
        else if (head == "phi_step")
        {
            phi_step = atof(value.c_str());
        }
        else if (head == "csi_max")
        {
            csi_max = atof(value.c_str());
        }
        else if (head == "csi_step")
        {
            csi_step = atof(value.c_str());
        }
		
		else if (head == "n_x_division")
		{
			n_x_division = atoi(value.c_str());
		}
		else if (head == "n_y_division")
		{
			n_y_division = atoi(value.c_str());
		}
		else if (head == "n_z_division")
		{
			n_z_division = atoi(value.c_str());
		}

        else if (head == "dx_step")
        {
            dx_step = atof(value.c_str());
        }
        else if (head == "dy_step")
        {
            dy_step = atof(value.c_str());
        }
        else if (head == "dz_step")
        {
            dz_step = atof(value.c_str());
        }
        else if (head == "dx_max")
        {
            dx_max = atof(value.c_str());
        }
        else if (head == "dy_max")
        {
            dy_max = atof(value.c_str());
        }
        else if (head == "dz_max")
        {
            dz_max = atof(value.c_str());
        }
        else if (head == "dim_block_x")
        {
            dim_block_x = atof(value.c_str());
        }
        else if (head == "dim_block_y")
        {
            dim_block_y = atof(value.c_str());
        }
        else if (head == "dim_block_z")
        {
            dim_block_z = atof(value.c_str());
        }
        
        else if (head == "read_block_dimension")
        {
            read_block_dimension = atoi(value.c_str());
        }

        else if (head == "read_bound")
        {
            read_bound = atoi(value.c_str());
        }
        
        else if (head == "BiDimensional")
        {
            BiDimensional = atoi(value.c_str());
        }
        else if (head == "write_vtu")
        {
            write_vtu = atoi(value.c_str());
        }
        else if (head == "read_PLY_FileList")
        {
            read_PLY_FileList = atoi(value.c_str());
        }
		else if (head == "n_x_division")
		{
			n_x_division = atoi(value.c_str());
		}
		else if (head == "n_y_division")
		{
			n_y_division = atoi(value.c_str());
		}
		else if (head == "n_z_division")
		{
			n_z_division = atoi(value.c_str());
		}
		else if (head == "rotation_method")
		{
            rotation_method = atoi(value.c_str());
        }
        else if (head == "cut_saw_thickness")
        {
            cut_saw_thickness = atof(value.c_str());
            alf_cut_saw_thickness = cut_saw_thickness / 2.0;
        }
		else if (head == "angle_type")
		{
			angle_type = atoi(value.c_str());
			
		}
		else if (head == "bound_tolerance")
		{
			bound_tolerance = atof(value.c_str());

		}
		
        else if (head == "end")
        {
            break;
        }
        //printf("\n%s", lines.c_str());
        //printf("\n");
        if (EndFile)break;
    }
    fclose(f_in);
    return 0;
}
int readPLY_FileList(void)
{
    int verbose = 1;
    FILE *fp_PLY_FileList;
    errno_t err;
    std::string lines(10000, ' ');
    std::string head, value;
    int EndFile;
    err = fopen_s(&fp_PLY_FileList, "PLY_FileList.dat", "r");
    if (err)
    {
        printf("PLY_FileList.dat file open error, program will exit with code 5\n");
        return 5;
    }
    int n_line = 0;
    //READING 
    //Example of file (no header, 3 double for each line).
    //0.0 0.0 0.0
    //1.0 0.0 0.0
    //...
    while (true)
    {


        EndFile = fgetstringline(lines, fp_PLY_FileList);
        if (lines.length() < 2)break;
        n_line++;
        if (verbose)printf("Reading line number: %d %s\n", n_line, lines.c_str());


        //string s = "What is the right way to split a string into a vector of strings";
        
        ply_file_name.push_back(lines);
    }

    fclose(fp_PLY_FileList);
    
    return 0;
}
int read_blocks_dimensions()
{
    if (verbose)printf("read_ending_points()\n");
    FILE *fp_ending_points;
    errno_t err;
    std::string lines(10000, ' ');
    std::string head, value;
    int EndFile;
    err = fopen_s(&fp_ending_points, "slab_dimensions.dat", "r");
    if (err)
    {
        printf("seed_points.dat file not found\n");
        return 5;
    }
    int n_line = 0;
    //READING 
    //Example of file (no header, 3 double for each line).
    //0.0 0.0 0.0
    //1.0 0.0 0.0
    //...
    while (true)
    {


        EndFile = fgetstringline(lines, fp_ending_points);
        if (lines.length() < 2)break;
        n_line++;
        if (verbose>2)printf("Reading line number: %d %s\n", n_line, lines.c_str());


        //string s = "What is the right way to split a string into a vector of strings";
        std::stringstream ss(lines);
        std::istream_iterator<std::string> begin(ss);
        std::istream_iterator<std::string> end;
        std::vector<std::string> vstrings(begin, end);
        Block_dimensions tmp_block_dim;
        tmp_block_dim.dim_x = atof(vstrings[0].c_str());
        tmp_block_dim.dim_y = atof(vstrings[1].c_str());
        tmp_block_dim.dim_z = atof(vstrings[2].c_str());
        blocks_dimensions_vector.push_back(tmp_block_dim);
        
    }

    fclose(fp_ending_points);

    return 0;
}
int fgetstringline(std::string &line, FILE* f_in) {
    //A file created under windows but read under linux has problem
    //with escape command. Linux ends text line with '0A' (\n, NL or LF ASCII), 
    //whereas Windows ends text line with '0D0A' (\r\n, CRNL ASCII). 
    int verbose = 0;
    int c;
    line.clear();
    do
    {
        c = fgetc(f_in);
        if (c == '\n')
        {
            if (verbose == 1)printf("fgetstringline_01\n");
            line.append(1, '\x0');
            if (line.size()>0)
            {                                           //under linux an empty text line has only NL char
                if (line[line.size() - 1] == '\r')
                {
                    if (verbose == 1)printf("fgetstringline_02\n");
                    line[line.size() - 1] = '\x0';
                } //for a text file created under Windows
            }
            break;
        }
        else
        {
            if (verbose == 1)printf("fgetstringline_03\n");
            if (c == EOF)
            {
                if (verbose == 1)printf("fgetstringline_04\n");
                line.append(1, '\x0');
                if (verbose == 1)printf("fgetstringline_04_A\n");
                return 1;
            }
            else
            {
                //line.append<int>(1,c);
                if (verbose == 1)printf("fgetstringline_05\n");
                line.append(1, (char)c);
            }
        }
    } while (true);
    return 0;
}

std::string EliminateCommentAndExtract(std::string str, std::string &head, std::string &value) {
    size_t pos1, pos2;
    std::string lines;
    //erases possible initial blanks
    pos1 = 0;
    for (unsigned int j = 0; j<str.size(); j++) {
        if (str[j] == ' ') { pos1++; }
        else { break; }
    }
    pos2 = str.find(" ", pos1);
    str = str.substr(pos1, (pos2 - pos1));
    //
    pos1 = str.find("#");
    if (pos1 == 0) return "";
    if (pos1 != str.npos) {
        str = str.substr(pos1, str.size() - pos1);
        str.append(1, '\x0');
    }
    //
    pos1 = str.find("=");
    if (pos1 == str.npos) return "";
    head = str.substr(0, pos1);
    value = str.substr(pos1 + 1, str.size() - pos1);
    //     return lines.substr(pos1,(pos2-pos1)); //exports the time value
    return str;
}

std::vector<std::string> ParseLineSpace(std::string str) {
    std::vector<std::string> tmp_string;
    tmp_string.reserve(4);
    std::size_t pos1, pos2;
    std::string lines;
    //erases possible initial blanks
    pos1 = 0;
    for (unsigned long int j = 0; j<str.size(); j++) {
        if (str[j] == ' ') { pos1++; }
        else { break; }
    }

    while (pos1<str.length()) {
        ///
        pos2 = str.find(" ", pos1);
        if (pos2 == std::string::npos)
        {
            tmp_string.push_back(str.substr(pos1, str.length()));
            break;
        }
        tmp_string.push_back(str.substr(pos1, (pos2 - pos1)));
        //     i++;
        pos1 = pos2 + 1;
    }
    //     return lines.substr(pos1,(pos2-pos1)); //exports the time value
    return tmp_string;
}


void read_BOUNDS_file(FILE *f_in)
{
    fgetstringline(lines, f_in);
    double x,y;
    n_verts_bound = atof(lines.c_str());
    verts_bound_x = new double[n_verts_bound];
    verts_bound_y = new double[n_verts_bound];
    for (int i = 0; i < n_verts_bound; i++) 
    {
        fgetstringline(lines, f_in);
        if (2 == sscanf(lines.c_str(), "%lf%lf", &x, &y))
        {
            verts_bound_x[i] = x;
            verts_bound_y[i] = y;
            
        }
    }

}

void read_PLY_file(FILE *f_in)
{

    std::string lines(10000, ' ');
    bool end_header = false;
    long n_vertex = 0;
    long n_faces = 0;
    std::string str_n_vertex("element vertex");
    std::string str_n_faces("element face");
    std::string str_end_header("end_header");
    std::string str_comment("comment");
    while (!end_header)
    {
        fgetstringline(lines, f_in);
        std::size_t found = lines.find(str_comment);
        if (found != std::string::npos)
        {
            continue;
        }
        
        found = lines.find(str_n_vertex);
        if(found != std::string::npos)
        {
            lines = SplitLine(lines, str_n_vertex);
            n_vertex=atoi(lines.c_str());
        }

        found = lines.find(str_n_faces);
        if (found != std::string::npos)
        {
            lines = SplitLine(lines, str_n_faces);
            n_faces = atoi(lines.c_str());
            n_triangles = n_faces;
        }
        found = lines.find(str_end_header);
        if (found != std::string::npos)
        {
            end_header = true;
        }
    }
    
    
    
    double x, y, z;
    Point *vertex_aray = new Point[n_vertex];
    for (int i = 0; i < n_vertex; i++)
    {
        fgetstringline(lines, f_in);
        if (3 == sscanf(lines.c_str(), "%lf%lf%lf", &x, &y, &z)) 
        {
            vertex_aray[i].x = x;
            vertex_aray[i].y = y;
            vertex_aray[i].z = z;
        }
    }
    //mesh_to_approx = new Triangle2[n_faces];
    
    int n_f,p_i,r_i,q_i;
    for (int i = 0; i < n_faces; i++)
    {
        fgetstringline(lines, f_in);
        if (4 == sscanf(lines.c_str(), "%d%d%d%d",&n_f, &p_i, &q_i,&r_i))
        {
            Triangle2 tmp_mesh_to_approx;
            tmp_mesh_to_approx.p[0] = vertex_aray[p_i].x;
            tmp_mesh_to_approx.p[1] = vertex_aray[p_i].y;
            tmp_mesh_to_approx.p[2] = vertex_aray[p_i].z;
            tmp_mesh_to_approx.q[0] = vertex_aray[q_i].x;
            tmp_mesh_to_approx.q[1] = vertex_aray[q_i].y;
            tmp_mesh_to_approx.q[2] = vertex_aray[q_i].z;
			tmp_mesh_to_approx.r[0] = vertex_aray[r_i].x;
            tmp_mesh_to_approx.r[1] = vertex_aray[r_i].y;
            tmp_mesh_to_approx.r[2] = vertex_aray[r_i].z;
            mesh_to_approx.push_back(tmp_mesh_to_approx);
        }
        
    }




}
std::string SplitLine(std::string input_line, std::string split_str)
{
    long int pos2 = input_line.find(split_str, 0);
    input_line = input_line.substr(pos2+ split_str.size(), input_line.size() - (pos2 + split_str.size()));
    return input_line;

}