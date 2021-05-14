#include "Definitions.h"
#include "GlobalVariables.h"
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

double Distance(Point p1, Point p2);  // Функция, считает расстояние между двумя точками
double Distance(Coordinates p1, Coordinates p2);  // Функция, считает расстояние между двумя точками
double RandGauss(double M, double sigma);

double Potential(double dist)  // Функция, возвращающая модуль потенциала Леннарда-Джонса по расстоянию между точками
{
	double one_dist = 1 / dist;
	double pp = one_dist * one_dist;
	double ppp = pp * pp * pp;
	return 4 * ppp * (ppp - 1);
	//return 4 * (pow((1 / dist), 12) - pow((1 / dist), 6));
}

double Force(double dist)  // Функция, возвращающая модуль силы взаимодействия частиц по расстоянию между ними
{
	double one_dist = 1 / dist;
	double pp = one_dist * one_dist;
	double ppp = pp * pp * pp;
	return 24 * ppp * one_dist * (1 - 2 * ppp);
	//return 24 * (pow(1 / dist, 7) - 2 * pow(1 / dist, 13));
}

void AddWallsForce(Point* p, int CutRadius) // Добавляем силу от стен
{
	double one_dist;
	double pp;
	double ppp;
	if (p->x <= CutRadius) {
		one_dist = 1 / p->x;
		pp = one_dist * one_dist;
		ppp = pp * pp * pp * pp * pp;
		p->ax += 9 * ppp;
	}
	if (Dim - p->x <= CutRadius) {
		one_dist = 1 / (Dim - p->x);
		pp = one_dist * one_dist;
		ppp = pp * pp * pp * pp * pp;
		p->ax -= 9 * ppp;
	}
	if (p->y <= CutRadius) {
		one_dist = 1 / p->y;
		pp = one_dist * one_dist;
		ppp = pp * pp * pp * pp * pp;
		p->ay += 9 * ppp;
	}
	if (Dim - p->y <= CutRadius) {
		one_dist = 1 / (Dim - p->y);
		pp = one_dist * one_dist;
		ppp = pp * pp * pp * pp * pp;
		p->ay -= 9 * ppp;
	}
	if (p->z <= CutRadius) {
		one_dist = 1 / p->z;
		pp = one_dist * one_dist;
		ppp = pp * pp * pp * pp * pp;
		p->az += 9 * ppp;
	}
	if (Dim - p->z <= CutRadius) {
		one_dist = 1 / (Dim - p->z);
		pp = one_dist * one_dist;
		ppp = pp * pp * pp * pp * pp;
		p->az -= 9 * ppp;
	}
}

double AddWallsPotential(Point* p, int CutRadius) // Считаем потенциальную энергию взаимодействия частицы со стенками
{
	double potential = 0;
	double one_dist, pp, ppp;
	if (p->x <= CutRadius) {
		one_dist = 1 / p->x;
		pp = one_dist * one_dist * one_dist;
		ppp = pp * pp * pp;
		potential += ppp;
	}
	if (Dim - p->x <= CutRadius) {
		one_dist = 1 / (Dim - p->x);
		pp = one_dist * one_dist * one_dist;
		ppp = pp * pp * pp;
		potential += ppp;
	}
	if (p->y <= CutRadius) {
		one_dist = 1 / p->y;
		pp = one_dist * one_dist * one_dist;
		ppp = pp * pp * pp;
		potential += ppp;
	}
	if (Dim - p->y <= CutRadius) {
		one_dist = 1 / (Dim - p->y);
		pp = one_dist * one_dist * one_dist;
		ppp = pp * pp * pp;
		potential += ppp;
	}
	if (p->z <= CutRadius) {
		one_dist = 1 / p->z;
		pp = one_dist * one_dist * one_dist;
		ppp = pp * pp * pp;
		potential += ppp;
	}
	if (Dim - p->z <= CutRadius) {
		one_dist = 1 / (Dim - p->z);
		pp = one_dist * one_dist * one_dist;
		ppp = pp * pp * pp;
		potential += ppp;
	}
	return potential;
}

Force_Calculation_Output Calculate_Force(Point *arr, double CutRadius, int EnergyMode, long long int Step)	// Функция подсчёта силы на текущем шаге
{
	double DimDiv2 = Dim / 2;
	double CurrentPotentialEnergy = 0;
	for (int i = 0; i < PointNumber; i++)         // Обнуляем ускорения перед каждым новым шагом
		arr[i].ax = arr[i].ay = arr[i].az = 0;
	for (int i = 0; i < PointNumber; i++) {              // Один шаг
		for (int j = i + 1; j < PointNumber; j++)
		{
			double dist = Distance(arr[i], arr[j]);         // Расстояние между частицами
			if (dist <= CutRadius)			// Считаем силы, если расстояние между частицами меньше радиуса обрезки потенциала
			{
				double force = Force(dist);
				arr[i].ax -= force * (arr[i].x - arr[j].x) / dist;	// Вычисляем изменение ускорения под действием частицы j
				arr[i].ay -= force * (arr[i].y - arr[j].y) / dist;
				arr[i].az -= force * (arr[i].z - arr[j].z) / dist;

				arr[j].ax += force * (arr[i].x - arr[j].x) / dist;	// Добавляем сразу же такую же силу к второй частице, но в обратном направлении
				arr[j].ay += force * (arr[i].y - arr[j].y) / dist;
				arr[j].az += force * (arr[i].z - arr[j].z) / dist;
				
				if (EnergyMode)
					CurrentPotentialEnergy += Potential(dist);  // Считаем потенциальную энергию взаимодействия
			}
		}
		AddWallsForce(&arr[i], CutRadius);
	}
	if (EnergyMode)
		for (int i = 0; i < PointNumber; i++) {
			CurrentPotentialEnergy += AddWallsPotential(&arr[i], CutRadius);
		}

	Force_Calculation_Output ret = { CurrentPotentialEnergy };
	return ret;
}