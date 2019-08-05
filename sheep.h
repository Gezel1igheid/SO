#ifndef SHEEP_h
#define SHEEP_h

#include <fstream>
#include <cstring>
#include <cstdio>
#include <ctime>
#include <cmath>

/********************************************************************************
**��Ⱥ�㷨V15.1
**����Ȯ�׶�Ϊ����Ⱥ��һ����������ֻ�������ɢ
**���汾����Ȯ�׶���δ��ɢ��Ⱥ���ɢ��Ⱥ���εĹ���
*********************************************************************************/
#define SHEEP_VERSION "15.1"				//��ǰ�㷨�İ汾��Ϣ������ʱע��ʹ��˫����

using  namespace std;
#define pai acos(-1.0)

#define DIM					30					//����ά��
#define SNUM				40					//��Ⱥ��ģ
#define ITE					100000				//��������
#define N					50					//ÿ��Ĳ�����Ŀ
#define SET					1					//��������
#define DEGREE_HUNTAWAY		0.2				//����Ȯ����̶�

#define U					0.000001			//��С��ֵ
#define PRINT_BELLWETHER	1					//�Ƿ����ļ�1�����ÿһ����ͷ��
#define PRINT_HUNTAWAY		1					//�Ƿ��������Ȯ����
#define TIMES_HUNTAWAY		0					//����Ȯ���������С���
#define PRINTF_RESULT_FILE	1					//�Ƿ�����ļ�2 �����ÿ�����ս��

/*������Ӧ�Ⱥ���ѡ��Ҫ���ĸ���������Ϊ1,��ֻ����һ��Ϊ1*/
#define TEST_FUN_CHOICE		"Griewank"	
//�ı���Ժ�������ͬʱ���˺궨���Ϊ��Ӧ���Ժ�����
#define Best_fitness		0.0					//�ı���Ժ�������ͬʱ���Ĵ˺�����ȫ�����Ž�
#define left_range			-600.00			//ÿһά�����귶Χ(����������������дΪ**.0����������600дΪ600.0)
#define right_range			600.00			//ÿһά�����귶Χ(����������������дΪ**.0����������600дΪ600.0)

//ע�⣺
//��ά����ֻ�ܲ��Զ�ά
//������ͨ�����ɲ�������ά��
//��תƽ�ƺ���ֻ�ܲ���10ά��30ά��50ά

//��ά����
#define Goldstein_and_Price				0		//range = 2.0,		Best_fitness = 3
#define Martin_and_Gaddy				0		//left_range = 0.0, right_range = 10.0,		Best_fitness = 0
#define six_hump_camel_back_function	0		//range 3.0,		Best_fitness = -1.031628
//���庯��
#define Sphere							0		//range = 100.0,	Best_fitness = 0
#define Sphere_Shifted					0		//range = 100.0,	Best_fitness = -450.0,	fbias = -450.0	ƽ�ƺ���
#define Axis_parallel_hyper_ellipsoid	0		//range = 5.12,		Best_fitness = 0
//��庯��
#define Rastrigin						0		//range = 5.12,		Best_fitness = 0
#define Rastrigin_Shifted				1		//range = 5.0,		Best_fitness = -330.0,	fbias = -330.0	ƽ�ƺ���
#define Rastrigin_Shifted_Rotated		0		//range = 5.0,		Best_fitness = -330.0,	fbias = -330.0	��תƽ�ƺ���
#define Ackley							0		//range = 32.0,		Best_fitness = 0
#define Ackley_Shifted_Rotated			0		//range = 32.0,		Best_fitness = -140.0,	fbias = -140.0	��תƽ�ƺ���
#define Griewank						1		//range = 600.0,	Best_fitness = 0
#define Schwefel						0		//range = 500.0,	Best_fitness = 0
#define Generalized_Penalized_1			0		//range = 50.0,		Best_fitness = 0
#define Generalized_Penalized_2			0		//range = 50.0,		Best_fitness = 0

//����ʽ���Ժ���
#define Demo							0
#define Rosenbrock						0		//range = 100.0,	Best_fitness = 0
class SHEEP
{
public:
	double	coordinate[DIM];					//�洢ÿֻ�������
	double	fitness;							//�洢�ʶ�ֵ
	int		number;								//�洢��ı��
	int		scatter;							//�ж���ֻ���Ƿ񱻴�ɢ 1����ɢ��0û����ɢ
public:
	friend class GROUPSHEEP;
};
class GROUPSHEEP
{
private:										//������Χ��Ҳ���ǽ�ռ䷶Χ
	double	left;								//�洢ÿά����ķ�Χ
	double	right;
public:
	SHEEP	sheep[SNUM];						//��Ⱥ
	int		bellwethernumber;					//��ͷ����
	double	worstfitness;						//��Ⱥ����ʶ�ֵ
	double	meanfitness;						//������Ⱥ���ʶ�ֵ��ƽ��ֵ
	double	oldbellwetherfitness;				//�洢��һ����ͷ����ʶ�ֵ
	int		generationTimes;
public:
	GROUPSHEEP();
	void	initofgroup();						//��ʼ����Ⱥ
	void	leader();							//��ͷ��׶�
	void	wander();							//��Ⱥ���ν׶�
	void	bellwether();						//����һ����ͷ��,���������Ⱥ��ƽ���ʶ�ֵ������ʶ�ֵ
	int		huntaway();							//����Ȯ�׶�
};

double Computafitness(double a[]);				//���Ժ���

#endif
