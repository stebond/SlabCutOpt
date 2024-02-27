#pragma once
#include <vector>
#include <math.h>
#include <string>



//#include "ply.h"

double xmin, xmax, ymin, ymax, zmin, zmax;
//Parameters to be optimized
double pigreco= atan(1.0) * 4.0;
double csi_step   = 10.0; 
double tetha_step = 10.0;
double phi_step   = 10.0;//z axis
double tetha_max  = pigreco/4.0;//improve this definition...


double phi_max = pigreco;
double csi_max = pigreco;

double dx_step = 10.0;
double dy_step = 10.0;
double dz_step = 10.0;
double dx_max = 0.2;
double dy_max = 0.2;
double dz_max = 0.2;

double dim_block_x;
double dim_block_y;
double dim_block_z;
int n_x_division = 1;
int n_y_division = 1;
int n_z_division = 1;

struct Block_dimensions
{
	double dim_x;
	double dim_y;
	double dim_z;

};

std::vector<Block_dimensions> blocks_dimensions_vector;

double bound_tolerance = 0.001;
double cut_saw_thickness = 0.0f;
double alf_cut_saw_thickness = 0.0f;
int n_triangles;
int read_bound = 0;
/////////////////////////////////////////
//options
int BiDimensional = 0;
int write_vtu = 1;
int read_PLY_FileList = 0;
//Rotation method: =1 is the classical  Eulero rotation. If =2, then are 3 succesive rotation around Y, around Z and around X axes.
int rotation_method = 1;
int read_block_dimension = 0;
//angle_type=0: radiant, default: angle_type=1, deg
int angle_type = 0;
std::vector<std::string> ply_file_name;
int verbose = 0;
struct Blocks {
	//double x, y, z;
	//double **P1;
	double *Center;
	bool intersected;
	bool inside;
	int rock_type;//=1 is inside;=0: outside;=2:intersected
};

//FUNCTIONS


//STRUCTURES
struct Point {
	double x, y, z;
};
struct Vector {
	double x, y, z;
};
struct Ray {
	Point P0, P1;
};
struct Triangle {
	Point V0, V1, V2;
};

struct Triangle2 {
	double p[3];
	double q[3];
	double r[3];
};
struct solution_struc {
	double dim_block_x, dim_block_y, dim_block_z, tetha, phi, csi, dx, dy, dz;
	int nx, ny, nz;
	double **P1;
	int n_block_inside = 0;
	int n_block_no_intersect = 0;
	int n_dimension;
	//to divide the domain in several  zones, the dimension of the space to be investigated can be
	//smaller than the whole domain
	double xmin, xmax, ymin, ymax, zmin, zmax;
	double X0;
	double Y0;
	double Z0;

	//to store in which domain the solution belongs, we have to store it:
	int i_domain;
	int j_domain;
	int k_domain;

};
struct best_results_struc {
	int solution_number;
	int n_block_inside;
	int n_block_no_intersect;
};
int intersect3D_RayTriangle(Ray R, Triangle T, Point *I, Vector *N);
Triangle *mesh_of_fractures;
int intersect3D_RayTriangle(Ray R, Triangle T, Point *I, Vector *N);
Vector Difference(Point v1, Point v2);
Vector Difference(Vector v1, Vector v2);
Vector crossP(Vector v1, Vector v2);
Vector ScalarVector(double r, Vector v1);
Point SumPoint(Point v1, Point v2);
Point SumPoint(Vector v1, Vector v2);
Point SumPoint(Point v1, Vector v2);
Vector SumVector(Vector v1, Vector v2);
double norm(Vector V);
Vector NormalizeV(Vector V);
bool isInside(double, double, double, solution_struc &);
bool isInside(double *&, solution_struc &);
double distance(double, double, double, double, double, double);
double min(double, double);
double max(double, double);
struct Blocks ***myBlocks;
//To eliminate:


//new functions:
int CreateBlocks3(solution_struc&, Blocks ***&);
int CreateBlocks5(solution_struc&, std::vector<std::vector<std::vector<Blocks>>>&);

int CheckBlocks3(solution_struc &, std::vector<std::vector<std::vector<Blocks>>>&);

int MoveBlocks3(solution_struc&, std::vector<std::vector<std::vector<Blocks>>>&);

int block_intersection2(solution_struc&,Blocks***&);
int block_intersection3(solution_struc&, std::vector<std::vector<std::vector<Blocks>>>&);

int freemem2(solution_struc&,Blocks***&);
int freemem3(solution_struc&, std::vector<std::vector<std::vector<Blocks>>>&);
int printVTU_Blocks2(char *, solution_struc&, Blocks ***&);
int printVTU_Blocks3(char *, solution_struc&, std::vector<std::vector<std::vector<Blocks>>>&);
int initialize(void);
int rotate_point(double, double, double,double *&/*in-out*/);
int rotate_point2(double, double, double, double *&/*in-out*/);

int readParameter(const char*);
int readPLY_FileList(void);
int read_blocks_dimensions(void);
std::string EliminateCommentAndExtract(std::string, std::string&, std::string&);
std::vector<std::string> ParseLineSpace(std::string);
std::string lines(10000, ' ');
int fgetstringline(std::string &, FILE*);
int CreateTriangles3(void);
int write_parameters(void);
int **cube_rect_faces;
int **cube_tri_faces;
std::vector<Triangle2> mesh_to_approx;
//Triangle2 *mesh_to_approx;
void read_PLY_file(FILE *);
std::string SplitLine(std::string, std::string);
int pnpoly(int , double *, double *, double, double);
long n_verts_bound;
double *verts_bound_x;
double *verts_bound_y;
void read_BOUNDS_file(FILE *);