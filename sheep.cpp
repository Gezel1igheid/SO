#include <iostream>
#include <set>
#include "sheep.h"


double GPu(double x, double a, double k, double m) {
	if (x > a)
		return k*pow((x - a), m);
	else if (x<-a) {
		return k*pow((-x - a), m);
	}
	else
		return 0;
}
double GPy(double x) {
	return 1 + 0.25 * (x + 1);
}
void mul_matrix(double z[DIM], double M[DIM][DIM]) {
	double z2[DIM];
	int i = 0, j = 0;
	for (int i = 0; i < DIM; i++) {
		z2[i] = z[i];
		z[i] = 0.0;
	}
	for (j = 0;j < DIM;j++) {
		for (i = 0;i < DIM;i++) {
			z[j] += z2[i] * M[i][j];
		}
	}
}

double Computafitness(double X[]) {
#if six_hump_camel_back_function
	//����ֵ�ձ��������ú�����6���ֲ���С��,����ȫ����С�� F��-0.0898��0.7126��=F��0.0898��-0.7126��=-1.031628
	return (4.0 - 2.1*X[0] * X[0] + (1.0 / 3.0)*X[0] * X[0] * X[0] * X[0])*X[0] * X[0] + X[0] * X[1] + (4 * X[1] * X[1] - 4)*X[1] * X[1];
#endif

#if Martin_and_Gaddy
	//Martin and Gaddy��������СֵΪ0.00000����ռ�Xi����[-20 20]
	return (X[0] - X[1])*(X[0] - X[1]) + ((X[0] + X[1] - 10) / 3.0)* ((X[0] + X[1] - 10) / 3.0);
#endif

#if Goldstein_and_Price
	return (1 + (X[0] + X[1] + 1)*(X[0] + X[1] + 1)*(19 - 14 * X[0] + 3 * X[0] * X[0] - 14 * X[1] + 6 * X[0] * X[1] + 3 * X[1] * X[1]))*(30 + (2 * X[0] - 3 * X[1])*(2 * X[0] - 3 * X[1])*(18 - 32 * X[0] + 12 * X[0] * X[0] + 48 * X[1] - 36 * X[0] * X[1] + 27 * X[1] * X[1]));
#endif
#if Demo
	double y = 0;
	y = X[0] * sin(10 * pai*X[0]) + 2.0;
	return y;
#endif

	//���Ե�Nά����
#if Sphere
	double f = 0;
	for (int i = 0; i < DIM; i++) {
		f = f + X[i] * X[i];
	}
	return f;
#endif

#if Sphere_Shifted
	double f = 0;
	extern double Sh_O[100];
	double z[DIM];
	for (int i = 0;i < DIM;i++) {
		z[i] = X[i] - Sh_O[i];
	}
	for (int i = 0; i < DIM; i++) {
		f = f + z[i] * z[i];
	}
	return f + -450.0;
#endif

#if Rastrigin
	double f = 0;
	for (int i = 0; i < DIM; i++) {
		f = f + (X[i] * X[i] - 10 * cos(2 * pai*X[i]) + 10);
	}
	return f;
#endif

#if Rastrigin_Shifted
	double f = 0;
	extern double Ra_O[100];
	double z[DIM];
	for (int i = 0;i < DIM;i++) {
		z[i] = X[i] - Ra_O[i];
	}
	for (int i = 0; i < DIM; i++){
		f = f + (z[i] * z[i] - 10 * cos(2 * pai*z[i]) + 10);
	}
	return f + -330;
#endif


#if Rastrigin_Shifted_Rotated
	double f = 0;
	extern double Ra_O[100];
	extern double Ra_M10[10][10], Ra_M30[30][30], Ra_M50[50][50];
	double z[DIM];
	for (int i = 0;i < DIM;i++) {
		z[i] = X[i] - Ra_O[i];
	}
#if DIM == 10
		mul_matrix(z, Ra_M10);
#endif
#if DIM == 30
		mul_matrix(z, Ra_M30);
#endif
#if DIM == 50
		mul_matrix(z, Ra_M50);
#endif
	for (int i = 0; i < DIM; i++){
		f = f + (z[i] * z[i] - 10 * cos(2 * pai*z[i]) + 10);
	}
	return f + -330.0;
#endif

#if Schwefel
	double  f = 0;

	for (int i = 0; i < DIM; i++) {
		f = f + (X[i] * sin(sqrt(fabs(X[i]))));
	}
	f = f + DIM*418.9829;
	return f;
#endif
#if Ackley
	double result;
	double f1 = 0, f2 = 0;
	for (int i = 0; i < DIM; i++) {
		f1 += X[i] * X[i];
		f2 += cos(2 * pai * X[i]);
	}
	result = -20 * exp(-0.2*sqrt((1.0 / DIM) * f1)) - exp((1.0 / DIM) * f2) + exp(1) + 20;
	return result;
#endif
#if Ackley_Shifted_Rotated
	double result;
	double f1 = 0, f2 = 0;
	extern double Ac_O[100];
	extern double Ac_M10[10][10], Ac_M30[30][30], Ac_M50[50][50];
	double z[DIM];
	for (int i = 0;i < DIM;i++) {
		z[i] = X[i] - Ac_O[i];
	}
#if DIM == 10
	mul_matrix(z, Ac_M10);
#endif
#if DIM == 30
	mul_matrix(z, Ac_M30);
#endif
#if DIM == 50
	mul_matrix(z, Ac_M50);
#endif
	for (int i = 0; i < DIM; i++) {
		f1 += z[i] * z[i];
		f2 += cos(2 * pai * z[i]);
	}
	result = -20 * exp(-0.2*sqrt((1.0 / DIM) * f1)) - exp((1.0 / DIM) * f2) + exp(1) + 20;
	return result + -140.0;
#endif

#if Griewank
	double result;
	double f1 = 0, f2 = 1;
	for (int i = 0; i < DIM; i++) {
		f1 += X[i] * X[i] / 4000;
		f2 *= cos(X[i] / sqrt(i + 1));
	}
	result = f1 - f2 + 1;
	return result;
#endif
#if Rosenbrock
	double result;
	double f1 = 0;
	for (int i = 0; i < DIM - 1; i++) {
		f1 += 100 * (X[i] * X[i] - X[i + 1]) * (X[i] * X[i] - X[i + 1]) + (X[i] - 1) * (X[i] - 1);
	}
	result = f1;
	return result;
#endif
#if Generalized_Penalized_1
	double result;
	double f1 = 0, f2 = 0;
	for (int i = 0; i < DIM - 1; i++) {
		f1 += ((GPy(X[i]) - 1)*(GPy(X[i]) - 1)*(1 + 10 * (sin(pai*GPy(X[i + 1]))) * (sin(pai*GPy(X[i + 1])))));
		f2 += GPu(X[i], 10, 100, 4);
	}
	f2 += GPu(X[DIM - 1], 10, 100, 4);
	result = pai / DIM*(10 * sin(pai*GPy(X[0]))*sin(pai*GPy(X[0])) + f1 + (GPy(X[DIM - 1]) - 1) * (GPy(X[DIM - 1]) - 1)) + f2;
	return result;
#endif
#if Generalized_Penalized_2
	double result;
	double f1 = 0, f2 = 0;
	for (int i = 0; i < DIM - 1; i++) {
		f1 += ((X[i] - 1)*(X[i] - 1)*(1 + 10 * (sin(pai*X[i + 1] * 3)) * (sin(pai*X[i + 1] * 3))));
		f2 += GPu(X[i], 5, 100, 4);
	}
	f2 += GPu(X[DIM - 1], 5, 100, 4);
	result = 0.1*(sin(3 * pai*X[0])*sin(3 * pai*X[0]) + f1 + (X[DIM - 1] - 1) * (X[DIM - 1] - 1) * (1 + (sin(2 * pai*X[DIM - 1]))*(sin(2 * pai*X[DIM - 1])))) + f2;
	return result;
#endif

#if Axis_parallel_hyper_ellipsoid
	double f = 0;
	for (int i = 0; i < DIM; i++) {
		f = f + X[i] * X[i] * (i + 1);
	}
	return f;
#endif
}

//�޲����Ĺ��캯����ȷ����ռ�ķ�Χ
GROUPSHEEP::GROUPSHEEP() {
	//ÿά��������޺�����
	left  = left_range;
	right = right_range;
}

void GROUPSHEEP::initofgroup() {					//��ʼ����Ⱥ
	generationTimes = 1;
	for (int i = 0; i < SNUM; i++) {				//��ʼ����Ⱥ
		for (int j = 0; j < DIM; j++) {				//��ʼ��ÿά�����ֵ
			sheep[i].coordinate[j] = rand() / (double)RAND_MAX*(right - left) + left;	//�����ÿһֻ�������
		}
		sheep[i].fitness = Computafitness(sheep[i].coordinate);							//�����һ����Ⱥÿ������ʶ�ֵ
		sheep[i].number = i + 1;					//����Ⱥ��ÿ�������ţ���ű������±��1
	}
}
void GROUPSHEEP::leader() {
	double loc[DIM], res;
	for (int i = 0; i<SNUM; i++) {
		sheep[i].scatter = 0;
		for (int j = 0; j<DIM; j++) {
			loc[j] = sheep[i].coordinate[j] + rand() / (double)RAND_MAX*(sheep[bellwethernumber - 1].coordinate[j] - sheep[i].coordinate[j]);
			if (loc[j]<left)
				loc[j] = right - left + loc[j];
			else if (loc[j]>right)
				loc[j] = loc[j] - right + left;
		}
		res = Computafitness(loc);
		if (res<sheep[i].fitness) {					//��λ�����λ�ý��бȽϣ�ȡ����Ϊ��
			for (int k = 0; k<DIM; k++) {
				sheep[i].coordinate[k] = loc[k];
			}
			sheep[i].fitness = res;
		}
	}
	generationTimes++;
}
void GROUPSHEEP::wander() {
	int i, j, random;
	double loc[DIM], res;
	for (i = 0; i<SNUM; i++) {						//ÿֻ����������ѧϰһ��
		do {
			random = (int)(rand() / (double)RAND_MAX*SNUM);
		} while (random == i || random == SNUM);	//��ѧϰ����������Ǽ���

		if (sheep[i].fitness < sheep[random].fitness) {
			for (j = 0; j < DIM; j++) {
				loc[j] = sheep[random].coordinate[j] + rand() / (double)RAND_MAX*(sheep[i].coordinate[j] - sheep[random].coordinate[j]);
				if (loc[j]<left)					//��������ռ�
					loc[j] = right - left + loc[j];
				else if (loc[j]>right)
					loc[j] = loc[j] - right + left;
			}
			res = Computafitness(loc);
			if (res<sheep[random].fitness) {		//��λ�����λ�ý��бȽ�
				for (j = 0; j<DIM; j++) {
					sheep[random].coordinate[j] = loc[j];
				}
				sheep[random].fitness = res;
			}
			for (j = 0; j < DIM; j++) {
				loc[j] = sheep[i].coordinate[j] + rand() / (double)RAND_MAX*(sheep[i].coordinate[j] - sheep[random].coordinate[j]);
				if (loc[j]<left)					//��������ռ�
					loc[j] = right - left + loc[j];
				else if (loc[j]>right)
					loc[j] = loc[j] - right + left;
			}
			res = Computafitness(loc);
			if (res<sheep[i].fitness) {				//��λ�����λ�ý��бȽ�
				for (j = 0; j<DIM; j++) {
					sheep[i].coordinate[j] = loc[j];
				}
				sheep[i].fitness = res;
			}
		}
		else {
			for (j = 0; j < DIM; j++) {
				loc[j] = sheep[i].coordinate[j] + rand() / (double)RAND_MAX*(sheep[random].coordinate[j] - sheep[i].coordinate[j]);
				if (loc[j]<left)					//��������ռ�
					loc[j] = right - left + loc[j];
				else if (loc[j]>right)
					loc[j] = loc[j] - right + left;
			}
			res = Computafitness(loc);
			if (res<sheep[i].fitness) {				//��λ�����λ�ý��бȽ�
				for (j = 0; j<DIM; j++) {
					sheep[i].coordinate[j] = loc[j];
				}
				sheep[i].fitness = res;
			}

			for (j = 0; j < DIM; j++) {
				loc[j] = sheep[random].coordinate[j] + rand() / (double)RAND_MAX*(sheep[random].coordinate[j] - sheep[i].coordinate[j]);
				if (loc[j]<left)					//��������ռ�
					loc[j] = right - left + loc[j];
				else if (loc[j]>right)
					loc[j] = loc[j] - right + left;
			}
			res = Computafitness(loc);
			if (res<sheep[random].fitness) {		//��λ�����λ�ý��бȽ�
				for (j = 0; j<DIM; j++){
					sheep[random].coordinate[j] = loc[j];
				}
				sheep[random].fitness = res;
			}
		}
	}
}

void GROUPSHEEP::bellwether() {						//�ҳ���һ�ε�������ͷ��
	double bellwetherfitness;
	double sum = 0;
	bellwetherfitness = sheep[0].fitness;
	bellwethernumber = 1;
	for (int k = 1; k < SNUM; k++) {
		if (sheep[k].fitness < bellwetherfitness) {
			bellwetherfitness = sheep[k].fitness;
			bellwethernumber = k + 1;				//������ͷ��ı��
		}
	}
	//��ƽ���ʶ�ֵ
	for (int i = 0; i < SNUM; i++)
		sum += sheep[i].fitness;
	meanfitness = sum / SNUM;
	//������ʶ�ֵ
	worstfitness = sheep[0].fitness;
	for (int j = 1; j<SNUM; j++) {
		if (sheep[j].fitness>worstfitness) {
			worstfitness = sheep[j].fitness;
		}
	}
}

int GROUPSHEEP::huntaway()
{
	static int lastHunt = 0;
	if (fabs(sheep[bellwethernumber - 1].fitness - oldbellwetherfitness) >= U) return 0;
	if (generationTimes < lastHunt) lastHunt = 0;
	if (generationTimes - lastHunt <= TIMES_HUNTAWAY) return 0; //������ϴ�����Ȯ�����������С�ڵ�������Ȯ�����������ִ������Ȯ����

	int i, j, random;
	double loc[DIM], res;

	int isscatter = 0;

	lastHunt = generationTimes;		//�������һ������Ȯ�������

	for (i = 0; i < SNUM; i++)
	{
		if (i + 1 != bellwethernumber)			//������ͷ�򲻱�
		{
			if (rand() / (double)RAND_MAX < DEGREE_HUNTAWAY)
			{
				sheep[i].scatter = 1;
				isscatter++;
				//printf("%d�����´�ɢ\n", i);
				for (j = 0; j < DIM; j++)		//��ɢ��Ⱥ
				{
					sheep[i].coordinate[j] = rand() / (double)RAND_MAX*(right - left) + left;
				}
				sheep[i].fitness = Computafitness(sheep[i].coordinate);
			}
		}
	}
	//fprintf(fp,"%dֻ�򱻴�ɢ\n", isscatter);

	for (i = 0; i < SNUM; i++)
	{
		if (isscatter == 0)
			return 0;
		if (sheep[i].scatter == 1)
			continue;
		while (1)
		{
			random = (int)(rand() / (double)RAND_MAX*(SNUM));
			if (sheep[random].scatter == 1)
				break;
		}
		if (sheep[i].fitness <= sheep[random].fitness)
		{
			//���Ҫѧϰ������Լ����ʶ�ֵ�����һ������ô����������
			for (j = 0; j < DIM; j++)
			{
				//	printf("%.6lf  %d\n", sheep[i].fitness,sheep[i].scatter);
				//	printf("%.6lf  %d\n", sheep[random].fitness, sheep[random].scatter);
				loc[j] = sheep[random].coordinate[j] + rand() / (double)RAND_MAX*(sheep[i].coordinate[j] - sheep[random].coordinate[j]);
				//printf("%.6lf", rand() / (double)RAND_MAX);
				if (loc[j]<left)			//����λ�ó�����ռ�����
					loc[j] = right - left + loc[j];
				else if (loc[j]>right)
					loc[j] = loc[j] - right + left;
			}
			res = Computafitness(loc);		//�����ʶ�ֵ
			if (res<sheep[random].fitness)	//��λ�����λ�ý��бȽ�
			{
				for (j = 0; j<DIM; j++)
				{
					sheep[random].coordinate[j] = loc[j];
				}
				//sheep[i].fitness = res;
				sheep[random].fitness = res;
			}


			for (j = 0; j < DIM; j++)
			{
				loc[j] = sheep[i].coordinate[j] + rand() / (double)RAND_MAX*(sheep[i].coordinate[j] - sheep[random].coordinate[j]);
				if (loc[j]<left)					//��������ռ�
					loc[j] = right - left + loc[j];
				else if (loc[j]>right)
					loc[j] = loc[j] - right + left;
			}
			res = Computafitness(loc);
			if (res<sheep[i].fitness)			//��λ�����λ�ý��бȽ�
			{
				for (j = 0; j<DIM; j++)
				{
					sheep[i].coordinate[j] = loc[j];
				}
				sheep[i].fitness = res;
			}

		}
		else	//���Ҫѧϰ������Լ����ʶ�ֵҪ�ã���ô������
		{
			for (j = 0; j < DIM; j++)
			{
				loc[j] = sheep[i].coordinate[j] + rand() / (double)RAND_MAX*(sheep[random].coordinate[j] - sheep[i].coordinate[j]);
				//printf("%.6lf", rand() / (double)RAND_MAX);
				if (loc[j]<left)			//����λ�ó�����ռ�����
					loc[j] = right - left + loc[j];
				else if (loc[j]>right)
					loc[j] = loc[j] - right + left;
			}
			res = Computafitness(loc);	//�����ʶ�ֵ
			if (res < sheep[i].fitness)	//��λ�����λ�ý��бȽ�
			{
				for (j = 0; j < DIM; j++)
				{
					sheep[i].coordinate[j] = loc[j];
				}
				sheep[i].fitness = res;
			}

			for (j = 0; j < DIM; j++)
			{
				loc[j] = sheep[random].coordinate[j] + rand() / (double)RAND_MAX*(sheep[random].coordinate[j] - sheep[i].coordinate[j]);
				if (loc[j]<left)					//��������ռ�
					loc[j] = right - left + loc[j];
				else
					if (loc[j]>right)
						loc[j] = loc[j] - right + left;
			}
			res = Computafitness(loc);
			if (res<sheep[random].fitness)			//��λ�����λ�ý��бȽ�
			{
				for (j = 0; j<DIM; j++)
				{
					sheep[random].coordinate[j] = loc[j];
				}
				sheep[random].fitness = res;
			}

		}
	}
	return 1;
}
