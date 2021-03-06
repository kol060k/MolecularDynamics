#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <string>

using namespace std;

double Density = 0.5;
double Temperature = 1.5;
int PointNumber = 1000;
int LatticeLength = 10;
double Model = 2;	// 1 - Периодически-граничные условия // 2 - Мягкие стенки //
double MinDistance = 0.8;
string LatticeType = "rand";
// Доступные: rand, sc, fcc
// Необходимые условия на количество частиц: 
	// sc: PointNumber = k^3, k - натуральное число
	// fcc: PointNumber = 4 * k^3, k - натуральное число

double Wall_MinDistance_Factor = 0.9; // Множитель минимального расстояния до стен для генерируемых частиц относительно переменной MinDistance

struct Point {
	double x, y, z;
	double vx, vy, vz;
};

double RandGauss(double M, double sigma)
{
	double a, b, r, Sq;
	do {
		a = (double)rand() / RAND_MAX * 2 - 1;
		b = (double)rand() / RAND_MAX * 2 - 1;
		r = a * a + b * b;
	} while (r > 1);
	Sq = sqrt(-2 * log(r) / r);
	return M + sigma * a * Sq;
	/*	double a = (double)rand() / RAND_MAX;
	while (a == 0.0) // Если a оказалось равно нулю
	a = (double)rand() / RAND_MAX;
	double b = (double)rand() / RAND_MAX;
	return sigma * sqrt(-log(a)) * cos(2 * M_PI * b) + M; */
	/*	double a, b; // M = 0, sigma = 1; иначе нужно пересчитывать
	do
	{
	a = ((double)rand() / RAND_MAX) * 0.4; // Значение по y - от 0 до 0.4
	b = ((double)rand() / RAND_MAX - 0.5) * 6; // Значение по x - от -3 до +3
	} while (1 / sqrt(2 * M_PI) * exp(-b * b / 2) < a); // Если сгенерированное в прямоугольнике число лежит под кривой, то выводим b. Т.о. случайная величина b имеет Гауссово распределение
	return b; */
}

double Distance(Point p1, Point p2)
{
	return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z));
}

void RandomCoordinates(Point *a, double Dim)
{
	a->x = (double)rand() / RAND_MAX * Dim;
	a->y = (double)rand() / RAND_MAX * Dim;
	a->z = (double)rand() / RAND_MAX * Dim;
}

void scCoordinates(Point *a)
{

}

void fccCoordinates(Point *a)
{

}

void RandomSpeeds(Point *a)
{
/*	a->vx = ((double)rand() / RAND_MAX - 0.5) * 4;
	a->vy = ((double)rand() / RAND_MAX - 0.5) * 4;
	a->vz = ((double)rand() / RAND_MAX - 0.5) * 4; */
	a->vx = RandGauss(0, 1);
	a->vy = RandGauss(0, 1);
	a->vz = RandGauss(0, 1);
}


int main()
{
	if (LatticeType == "sc")
		PointNumber = pow(LatticeLength, 3);
	if (LatticeType == "fcc")
		PointNumber = 4 * pow(LatticeLength, 3);

	srand(time(NULL));
	Point *arr;
	arr = (Point*)malloc(PointNumber * sizeof(Point));
	double Dim = 0;

//// Задание координат ////
	if (LatticeType == "rand") {
		Dim = pow(PointNumber / Density, (double)1 / 3);
		double DimDiv2 = Dim / 2;
		for (int i = 0; i < PointNumber; i++)  // Генерируем частицы
		{
			int err;
			do
			{
				RandomCoordinates(&arr[i], Dim);
				err = 0;

				if (Model == 1) {	// Для ПГУ
					for (int j = 0; j < i; j++)
					{
						Point image;              // Для каждой точки i находим ближайший из образов каждой точки j, чтобы проверить расстояние
						if (abs(arr[i].x - arr[j].x) <= (DimDiv2)) image.x = arr[j].x;  // По X
						else if (arr[i].x < arr[j].x) image.x = arr[j].x - Dim;
						else image.x = arr[j].x + Dim;
						if (abs(arr[i].y - arr[j].y) <= (DimDiv2)) image.y = arr[j].y;  // По Y
						else if (arr[i].y < arr[j].y) image.y = arr[j].y - Dim;
						else image.y = arr[j].y + Dim;
						if (abs(arr[i].z - arr[j].z) <= (DimDiv2)) image.z = arr[j].z;  // По Z
						else if (arr[i].z < arr[j].z) image.z = arr[j].z - Dim;
						else image.z = arr[j].z + Dim;

						double dist = Distance(arr[i], image);		// Вычисляем расстояние между одной частицей и ближайшим изображением второй
						if (dist < MinDistance)
							err = 1;
					}
				}

				if (Model == 2) {	// Для Мягких стенок
					if ((arr[i].x < MinDistance * Wall_MinDistance_Factor) || (Dim - arr[i].x < MinDistance * Wall_MinDistance_Factor)
						|| (arr[i].y < MinDistance * Wall_MinDistance_Factor) || (Dim - arr[i].y < MinDistance * Wall_MinDistance_Factor)
						|| (arr[i].z < MinDistance * Wall_MinDistance_Factor) || (Dim - arr[i].z < MinDistance * Wall_MinDistance_Factor))
						err = 1;
					else for (int j = 0; j < i; j++)
						if (Distance(arr[i], arr[j]) < MinDistance)
							err = 1;
				}
			} while (err);
		}
	}

	if (LatticeType == "sc") {
		Dim = LatticeLength / pow(Density, 1.0 / 3);
		double LatticePeriod = Dim / LatticeLength;
		for (int i = 0; i < LatticeLength; i++)
			for (int j = 0; j < LatticeLength; j++)
				for (int k = 0; k < LatticeLength; k++) {
					int num = i * LatticeLength * LatticeLength + j * LatticeLength + k;
					arr[num].x = LatticePeriod * i;
					arr[num].y = LatticePeriod * j;
					arr[num].z = LatticePeriod * k;
				}
	}

	if (LatticeType == "fcc") {
		Dim = pow(4 * pow(LatticeLength, 3) / Density, 1.0 / 3);
		double LatticePeriod = Dim / LatticeLength;
		for (int i = 0; i < LatticeLength; i++)
			for (int j = 0; j < LatticeLength; j++)
				for (int k = 0; k < LatticeLength; k++) {
					int num = 4 * (i * LatticeLength * LatticeLength + j * LatticeLength + k);
					arr[num].x = LatticePeriod * i;
					arr[num].y = LatticePeriod * j;
					arr[num].z = LatticePeriod * k;
					arr[num + 1].x = LatticePeriod * (i + 0.5);
					arr[num + 1].y = LatticePeriod * (j + 0.5);
					arr[num + 1].z = LatticePeriod * k;
					arr[num + 2].x = LatticePeriod * (i + 0.5);
					arr[num + 2].y = LatticePeriod * j;
					arr[num + 2].z = LatticePeriod * (k + 0.5);
					arr[num + 3].x = LatticePeriod * i;
					arr[num + 3].y = LatticePeriod * (j + 0.5);
					arr[num + 3].z = LatticePeriod * (k + 0.5);
				}
	}

//// Задание скоростей ////
	for (int i = 0; i < PointNumber; i++)
		RandomSpeeds(&arr[i]);

	long double CenterOfMass_Speed_X = 0, CenterOfMass_Speed_Y = 0, CenterOfMass_Speed_Z = 0;
	for (int i = 0; i < PointNumber; i++)  // Скорость центра масс
	{
		CenterOfMass_Speed_X += arr[i].vx;
		CenterOfMass_Speed_Y += arr[i].vy;
		CenterOfMass_Speed_Z += arr[i].vz;
	}
	CenterOfMass_Speed_X = CenterOfMass_Speed_X / PointNumber;
	CenterOfMass_Speed_Y = CenterOfMass_Speed_Y / PointNumber;
	CenterOfMass_Speed_Z = CenterOfMass_Speed_Z / PointNumber;

	for (int i = 0; i < PointNumber; i++)  // Вычитаем скорость центра масс
	{
		arr[i].vx -= CenterOfMass_Speed_X;
		arr[i].vy -= CenterOfMass_Speed_Y;
		arr[i].vz -= CenterOfMass_Speed_Z;
	}

	// Масштабируем скорости чтобы получить необходимую температуру
	double CurrentKineticEnergy = 0;
	for (int i = 0; i < PointNumber; i++) {
		CurrentKineticEnergy += (arr[i].vx * arr[i].vx + arr[i].vy * arr[i].vy + arr[i].vz * arr[i].vz) / 2;
	}
	double CurrentTemperature = 2 * CurrentKineticEnergy / 3 / PointNumber;
	double scale = sqrt(Temperature / CurrentTemperature);
	for (int i = 0; i < PointNumber; i++) {
		arr[i].vx *= scale;
		arr[i].vy *= scale;
		arr[i].vz *= scale;
	}

	ofstream fout("points.txt", ios_base::out | ios_base::trunc);
	fout << Dim << endl;
	fout << PointNumber << endl;
	fout << MinDistance << endl;
	for (int i = 0; i < PointNumber; i++)
		fout << fixed << std::setprecision(15) << arr[i].x << " " << arr[i].y << " " << arr[i].z << " " << arr[i].vx << " " << arr[i].vy << " " << arr[i].vz << " " << endl;
	fout.close();
	return 0;
}