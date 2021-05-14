#include "Definitions.h"
#include <cmath>

using namespace std;

double Distance(Point p1, Point p2)  // Функция, считает расстояние между двумя точками
{
	double deltaX = p1.x - p2.x;
	double deltaY = p1.y - p2.y;
	double deltaZ = p1.z - p2.z;
	return sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);
	//return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z));
}

double Distance(Coordinates p1, Coordinates p2)  // Функция, считает расстояние между двумя точками
{
	double deltaX = p1.x - p2.x;
	double deltaY = p1.y - p2.y;
	double deltaZ = p1.z - p2.z;
	return sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);
	//return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z));
}

double RandGauss(double M, double sigma)
{
	double a, b, r, Sq;		// В Frenkel Smit алгоритм ровно такой же
	do {
		a = (double)rand() / RAND_MAX * 2 - 1;
		b = (double)rand() / RAND_MAX * 2 - 1;
		r = a * a + b * b;
	} while (r > 1);
	Sq = sqrt(-2 * log(r) / r);
	return M + sigma * a * Sq;
}