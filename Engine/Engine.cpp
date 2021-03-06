#include <iostream>
#include <cmath>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <string>
#include "Definitions.h"

using namespace std;

double Dim = 0;
int PointNumber = 0;
double MinDistance = 0;
double PointRadius = 0.20;

double StepTime = 0.001;				// Для малого числа использовалось 0.0001
long long int StepNumber = 1000;
long long int BalanceStep = 1000;

long long int PointsOutPeriod = 100;
long long int RealPointsOutPeriod = 100;
long long int ScreenOutPeriod = 100;
long long int DiffusionOutPeriod = 100;
long long int EnergyAveragePeriod = 100;	
long long int TempAveragePeriod = 100;

double Distance(Point p1, Point p2);  // Функция, считает расстояние между двумя точками
double Distance(Coordinates p1, Coordinates p2);  // Функция, считает расстояние между двумя точками
double Potential(double dist);  // Функция, возвращающая модуль потенциала Леннарда-Джонса по расстоянию между точками
double Force(double dist);  // Функция, возвращающая модуль силы взаимодействия частиц по расстоянию между ними
void AddWallsForce(Point* p, int CutRadius); // Добавляем силу от стен
double AddWallsPotential(Point* p, int CutRadius); // Считаем потенциальную энергию взаимодействия частицы со стенками
double RandGauss(double M, double sigma);
Force_Calculation_Output Calculate_Force(Point *arr, double CutRadius, int EnergyMode, long long int Step);	// Функция подсчёта силы на текущем шаге


int main()
{
	// Модификаторы
	int EnergyMode = 0;           // Подсчёт полной энергии системы
	int TemperatureMode = 0;      // Подсчёт температуры системы
	int OutputMode = 0;           // Вывод результатов работы программы (пошагово)
	int MaxwellMode = 0;          // Обработка распределения Максвелла
	int DiffusionMode = 0;        // Обработка диффузии
	double CutRadius = 0;         // Радиус обрезания потенциала (0 - нет обрезания)
	int RealOutputMode = 0;       // Выводить реальные координаты частиц
	cout << "Mode settings (1 - ON, 0 - OFF):" << endl;
	cout << "Create output file: ";
	cin >> OutputMode;
	cout << "Energy calculation mode: ";
	cin >> EnergyMode;
	if (EnergyMode) {
		cout << "Temperature calculation mode: ";
		cin >> TemperatureMode;
	}
	else {
		cout << "Temperature calculation mode: 0 (Energy mode is necessary)" << endl;
	}
	cout << "Maxwell distribution calculation mode: ";
	cin >> MaxwellMode;
	cout << "Diffusion calculation mode: ";
	cin >> DiffusionMode;
	cout << "Potential cutoff radius (0 - no cutoff): ";
	cin >> CutRadius;
	cout << "Create output file with real coordinates: ";
	cin >> RealOutputMode;
	cout << endl;


// Считываем начальные данные точек газа. Другие начальные условия
	string points_str = "points.txt";
	ifstream fin(points_str.c_str(), ios_base::in); // Начальные положения             
	fin >> Dim;
	fin >> PointNumber;
	fin >> MinDistance;
	Point *arr;
	arr = new Point[PointNumber];
	for (int i = 0; i < PointNumber; i++)
		fin >> arr[i].x >> arr[i].y >> arr[i].z >> arr[i].vx >> arr[i].vy >> arr[i].vz;
	fin.close();

	if (CutRadius == 0)		// Если радиус обрезки не указан, берём его больше, чем область моделирования
		CutRadius = 5 * Dim;

	double CurrentEnergySumm = 0;	// Переменные, использующиеся при усреднении энергии
	double KineticEnergySumm = 0;
	double PotentialEnergySumm = 0;
	ofstream Energy_file;
	if (EnergyMode) {
		Energy_file.open("Output-Energy.txt", ios_base::out | ios_base::trunc);
	}

	double TemperatureSumm = 0;	// Переменная, использующаяся при усреднении температуры
	ofstream Temperature_file;
	if (TemperatureMode) {
		string Temperature_str = "Output-Temperature.txt";
		Temperature_file.open(Temperature_str.c_str(), ios_base::out | ios_base::trunc);
	}

	Coordinates *arr0 = NULL; // Для диффузии. Сохраняем начальные положения точек
	Coordinates *arr1 = NULL; // Для диффузии. Настоящие координаты точек
	Coordinates *arr01 = NULL; // То же самое. Используем для усреднения диффузии
	Coordinates *arr11 = NULL;
	Coordinates *arr02 = NULL;
	Coordinates *arr12 = NULL;
	ofstream Diffusion_file;


	long int *Speed_Dist = NULL; // Для Максвелла. Записываем количество частиц, имеющих скорость [Vx; Vx + dVx);
	double Speed_Dist_Step = 0.02; // Шаг dVx
	int Speed_Dist_Size = 10 * 50; // Количество шагов dVx от нуля до предельной скорости
	if (MaxwellMode) {
		Speed_Dist = new long int[Speed_Dist_Size];
		for (int i = 0; i < Speed_Dist_Size; i++) {
			Speed_Dist[i] = 0;
		}
	}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Создаём .xyz файл с нужными опциями, вносим в него начальные положения точек (если необходимо)
	ofstream fout;
	if (OutputMode) {
		fout.open("process.xyz", ios_base::out | ios_base::trunc);
		fout << PointNumber << endl;
		fout << "Lattice = \"1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0\" ";
		fout << "Properties=pos:R:3:velo:R:3:radius:R:1 ";
		fout << "Time = 0" << endl;
		for (int i = 0; i < PointNumber; i++)
			fout << arr[i].x << " " << arr[i].y << " " << arr[i].z << " " << arr[i].vx << " " << arr[i].vy << " " << arr[i].vz << " " << PointRadius << " " << endl;
	}
	
// .xyz файл с реальными координатами частиц. Вносим начальные данные, создаём массив для реальных координат
	Coordinates *arrReal = NULL;
	ofstream Realfout;
	if (RealOutputMode) {
		string Realfout_str = "RealPoints.xyz";
		Realfout.open(Realfout_str.c_str(), ios_base::out | ios_base::trunc);
		Realfout << PointNumber << endl;
		Realfout << "Lattice = \"1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0\" ";
		Realfout << "Properties=pos:R:3:velo:R:3:radius:R:1 ";
		Realfout << "Time = 0" << endl;
		for (int i = 0; i < PointNumber; i++)
			Realfout << arr[i].x << " " << arr[i].y << " " << arr[i].z << " " << arr[i].vx << " " << arr[i].vy << " " << arr[i].vz << " " << PointRadius << " " << endl;
		arrReal = new Coordinates[PointNumber];
		for (int i = 0; i < PointNumber; i++) {
			arrReal[i].x = arr[i].x;
			arrReal[i].y = arr[i].y;
			arrReal[i].z = arr[i].z;
		}
	}
//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Алгоритм
	srand(time(NULL));

	long double MinEnergy = 1000000000;
	long double MaxEnergy = -1000000000;

	double StepTimeDiv2 = StepTime / 2; // Для ускорения расчётов (уменьшение количества операций)
	double StepTimePower2 = StepTime * StepTime;
	double DimDiv2 = Dim / 2;

	unsigned int BeginTime = clock(); // Посчитаем время выполнения непосредственно большого цикла

	int PointOutOfCube_flag = 0; // Используется в модели с мягкими стенками. Флаг ошибочного выхода точки за границу

	Calculate_Force(arr, CutRadius, EnergyMode, 0);	// Первоначальный подсчёт силы (на нулевом шаге)

////////////////// БОЛЬШОЙ ЦИКЛ ОБРАБОТКИ КАЖДОГО ШАГА //////////////////////////////////////////
	for (long long int Step = 1; Step <= StepNumber; Step++)
	{
		long double CurrentKineticEnergy = 0;
		long double CurrentPotentialEnergy = 0;
		long double CurrentEnergy = 0;
		double CurrentTemperature = 0;

//// Для диффузии. Первоначальное задание необходимых массивов
		if (DiffusionMode) {
			if (Step == BalanceStep) {
				string Diffusion_str = "Output-Diffusion.txt";
				Diffusion_file.open(Diffusion_str.c_str(), ios_base::out | ios_base::trunc);
				Diffusion_file << "0\t0" << endl;
				arr0 = new Coordinates[PointNumber];
				arr1 = new Coordinates[PointNumber];
				for (int i = 0; i < PointNumber; i++) {
					arr0[i].x = arr1[i].x = arr[i].x;
					arr0[i].y = arr1[i].y = arr[i].y;
					arr0[i].z = arr1[i].z = arr[i].z;
				}
			}
		}

//// Схема Верле. Часть 1. r(t + dt)
		for (int i = 0; i < PointNumber; i++)
		{
			arr[i].x += arr[i].vx * StepTime + arr[i].ax * StepTimePower2 / 2;
			arr[i].y += arr[i].vy * StepTime + arr[i].ay * StepTimePower2 / 2;
			arr[i].z += arr[i].vz * StepTime + arr[i].az * StepTimePower2 / 2;
		}

	//// Проверяем, чтобы молекулы не выходили за границы куба
		for (int i = 0; i < PointNumber; i++)
			if ((arr[i].x < 0) || (arr[i].x > Dim) || (arr[i].y < 0) || (arr[i].y > Dim) || (arr[i].z < 0) || (arr[i].z > Dim))
			{
				PointOutOfCube_flag = 1;
				cout << arr[i].x << endl;
				cout << arr[i].y << endl;
				cout << arr[i].z << endl;
				goto ErrorExit;
			}

//// Для вывода реальных координат частиц (без сдвига через ПГУ)
		if (RealOutputMode) {
			for (int i = 0; i < PointNumber; i++) {
				arrReal[i].x += arr[i].vx * StepTime + arr[i].ax * StepTimePower2 / 2;
				arrReal[i].y += arr[i].vy * StepTime + arr[i].ay * StepTimePower2 / 2;
				arrReal[i].z += arr[i].vz * StepTime + arr[i].az * StepTimePower2 / 2;
			}
		}
//// Для диффузии. Расчёт настоящих положений точек без пересчёта выхода за границу
		//// Предложение по оптимизации - не считать два раза "arr[i].vx * StepTime + arr[i].ax * StepTimePower2 / 2" ////
		if (DiffusionMode) {
			if (Step >= BalanceStep) {
				long double Summ_SquareDist = 0;
				for (int i = 0; i < PointNumber; i++) {
					arr1[i].x += arr[i].vx * StepTime + arr[i].ax * StepTimePower2 / 2;
					arr1[i].y += arr[i].vy * StepTime + arr[i].ay * StepTimePower2 / 2;
					arr1[i].z += arr[i].vz * StepTime + arr[i].az * StepTimePower2 / 2;
					double deltaX = arr1[i].x - arr0[i].x; // Уменьшим количество операций
					double deltaY = arr1[i].y - arr0[i].y;
					double deltaZ = arr1[i].z - arr0[i].z;
					Summ_SquareDist += deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ;
				}
				if (Step % DiffusionOutPeriod == 0)
					Diffusion_file << (Step - BalanceStep + 1) * StepTime << "\t" << Summ_SquareDist / PointNumber << endl;
			}
		}

//// Схема Верле. Часть 2. v(t + dt)
		for (int i = 0; i < PointNumber; i++)	// Первая часть обновления скоростей
		{
			arr[i].vx += arr[i].ax * StepTimeDiv2;
			arr[i].vy += arr[i].ay * StepTimeDiv2;
			arr[i].vz += arr[i].az * StepTimeDiv2;
		}

	//// Подсчёт a(t + dt)
		Force_Calculation_Output CurrentOutput;
		CurrentOutput = Calculate_Force(arr, CutRadius, EnergyMode, Step);
		CurrentPotentialEnergy = CurrentOutput.CurrentPotentialEnergy;

		for (int i = 0; i < PointNumber; i++)	// Вторая часть обновления скоростей
		{
			arr[i].vx += arr[i].ax * StepTimeDiv2;
			arr[i].vy += arr[i].ay * StepTimeDiv2;
			arr[i].vz += arr[i].az * StepTimeDiv2;
		}

		if (EnergyMode)
			for (int i = 0; i < PointNumber; i++) {
				CurrentKineticEnergy += (arr[i].vx * arr[i].vx + arr[i].vy * arr[i].vy + arr[i].vz * arr[i].vz) / 2;  // Добавляем кинетическую энергию
			}

//// Для Максвелла. Находим распределение скоростей на этом шаге
		if ((MaxwellMode) && (Step >= BalanceStep)) {
			for (int i = 0; i < PointNumber; i++)
			{
				int j;
				double speed = sqrt(arr[i].vx * arr[i].vx + arr[i].vy * arr[i].vy + arr[i].vz * arr[i].vz); // Распределение по модулю скорости
				j = trunc(speed / Speed_Dist_Step);
				Speed_Dist[j]++;
			}
		}

		if (OutputMode) {
			if (Step % PointsOutPeriod == 0)  // Сохраняем не все кадры, а только каждый n-й
			{
				fout << PointNumber << endl;
				fout << "Lattice = \"1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0\" ";
				fout << "Properties=pos:R:3:velo:R:3:radius:R:1 ";
				fout << "Time = " << Step * StepTime << endl;
				for (int i = 0; i < PointNumber; i++)
					fout << arr[i].x << " " << arr[i].y << " " << arr[i].z << " " << arr[i].vx << " " << arr[i].vy << " " << arr[i].vz << " " << PointRadius << " " << endl;
			}
		}

		if (RealOutputMode) {
			if (Step % RealPointsOutPeriod == 0)  // Сохраняем не все кадры, а только каждый n-й
			{
				Realfout << PointNumber << endl;
				Realfout << "Lattice = \"1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0\" ";
				Realfout << "Properties=pos:R:3:velo:R:3:radius:R:1 ";
				Realfout << "Time = " << Step * StepTime << endl;
				for (int i = 0; i < PointNumber; i++)
					Realfout << arrReal[i].x << " " << arrReal[i].y << " " << arrReal[i].z << " " << arr[i].vx << " " << arr[i].vy << " " << arr[i].vz << " " << PointRadius << " " << endl;
			}
		}

		if (EnergyMode)
		{
			CurrentEnergy = CurrentKineticEnergy + CurrentPotentialEnergy;	// Вывод средней энергии на экран
			if (Step % ScreenOutPeriod == 0) {
				cout << "Current Energy = " << CurrentEnergy << endl;
				cout << "Current Kinetic Energy = " << CurrentKineticEnergy << endl;
				cout << "Current Potential Energy = " << CurrentPotentialEnergy << endl;
			}
			if (CurrentEnergy < MinEnergy) MinEnergy = CurrentEnergy; // Нахождение min и max энергии
			if (CurrentEnergy > MaxEnergy) MaxEnergy = CurrentEnergy;
			CurrentEnergySumm += CurrentEnergy;		// Усреднение энергии для ввода в файл
			KineticEnergySumm += CurrentKineticEnergy;
			PotentialEnergySumm += CurrentPotentialEnergy;
			if (Step % EnergyAveragePeriod == 0)
			{
				double AverageCurrentEnergy = CurrentEnergySumm / EnergyAveragePeriod;
				double AverageKineticEnergy = KineticEnergySumm / EnergyAveragePeriod;
				double AveragePotentialEnergy = PotentialEnergySumm / EnergyAveragePeriod;
				Energy_file << Step * StepTime << "\t" << AverageCurrentEnergy << "\t" << AverageKineticEnergy << "\t" << AveragePotentialEnergy << endl;
				CurrentEnergySumm = KineticEnergySumm = PotentialEnergySumm = 0;
			}
		}

//// Расчёт температуры
		if (TemperatureMode)
		{
			CurrentTemperature = CurrentKineticEnergy / PointNumber * 2 / 3;
			//cout << "Current Temperature = " << CurrentTemperature << endl;
			TemperatureSumm += CurrentTemperature;	// Для усреднения температуры
			if (Step % TempAveragePeriod == 0)
			{
				double AverageTemperature = TemperatureSumm / TempAveragePeriod;
				Temperature_file << Step * StepTime << "\t" << AverageTemperature << endl;
				TemperatureSumm = 0;
			} 
		}

		if (Step % ScreenOutPeriod == 0)
		{
			cout << "Step " << Step << " calculated" << endl;
			cout << (double(Step) / StepNumber) * 100 << "% ready" << endl;
			unsigned int CurrentTime = (clock() - BeginTime) / 1000;
			cout << CurrentTime / 3600 << ":" << (CurrentTime % 3600) / 60 << ":" << CurrentTime % 60 << " passed" << endl;
			cout << endl;
		}
	}

ErrorExit: // Сюда же выходим, если появилась какая-то ошибка
	
	if (PointOutOfCube_flag)
		cout << "Error!" << endl << "Some points are out of cube!" << endl;

	unsigned int EndTime = clock();
	unsigned int RunTime = (EndTime - BeginTime) / 1000; // Algorithm time in seconds

	if (EnergyMode)
	{
		cout << endl << "Minimal Energy in process: " << MinEnergy << endl;
		cout << "Maximal Energy in process: " << MaxEnergy << endl;
		cout << "Energy fluctuations during the process: " << MaxEnergy - MinEnergy << endl;
		cout << "Fluctuations of energy in comparison with average energy (average = (min + max) / 2): " << (MaxEnergy - MinEnergy) / (MaxEnergy + MinEnergy) << endl;
	}

	cout << "Runtime is " << RunTime / 3600 << ":" << (RunTime % 3600) / 60 << ":" << RunTime % 60 << endl;

	ofstream Maxwell_file;
	if (MaxwellMode) {
		Maxwell_file.open("Output-Maxwell.txt", ios_base::out | ios_base::trunc);
		for (int i = 0; i < Speed_Dist_Size; i++)
			Maxwell_file << (i * Speed_Dist_Step) * (i * Speed_Dist_Step) << "\t" << /*log((long double)Speed_Dist[i] / 3 / StepNumber)*/ Speed_Dist[i] << endl;
	}

	ofstream Points_out_file;
	Points_out_file.open("Output-points.txt", ios_base::out | ios_base::trunc);
	Points_out_file << Dim << endl;
	Points_out_file << PointNumber << endl;
	Points_out_file << MinDistance << endl;
	for (int i = 0; i < PointNumber; i++)
		Points_out_file << fixed << std::setprecision(15) << arr[i].x << " " << arr[i].y << " " << arr[i].z << " " << arr[i].vx << " " << arr[i].vy << " " << arr[i].vz << " " << endl;
	
	delete[] arr;
	delete[] arr0;
	delete[] arr1;
	delete[] arr01;
	delete[] arr11;
	delete[] arr02;
	delete[] arr12;
	delete[] Speed_Dist;
	delete[] arrReal;
	Energy_file.close();
	Diffusion_file.close();
	Maxwell_file.close();
	Temperature_file.close();
	fout.close();
	Realfout.close();
	Points_out_file.close();

	return 0;
}