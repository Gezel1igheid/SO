#ifndef SHEEP_h
#define SHEEP_h

#include <fstream>
#include <cstring>
#include <cstdio>
#include <ctime>
#include <cmath>

/********************************************************************************
**羊群算法V15.1
**牧羊犬阶段为将羊群中一定数量的整只羊随机打散
**本版本牧羊犬阶段有未打散羊群向打散羊群漫游的过程
*********************************************************************************/
#define SHEEP_VERSION "15.1"				//当前算法的版本信息，更改时注意使用双引号

using  namespace std;
#define pai acos(-1.0)

#define DIM					30					//粒子维度
#define SNUM				40					//种群规模
#define ITE					100000				//迭代次数
#define N					50					//每组的测试数目
#define SET					1					//测试组数
#define DEGREE_HUNTAWAY		0.2				//牧羊犬介入程度

#define U					0.000001			//最小阀值
#define PRINT_BELLWETHER	1					//是否在文件1中输出每一代领头羊
#define PRINT_HUNTAWAY		1					//是否输出牧羊犬介入
#define TIMES_HUNTAWAY		0					//牧羊犬介入代数最小间隔
#define PRINTF_RESULT_FILE	1					//是否输出文件2 仅输出每次最终结果

/*各种适应度函数选择，要用哪个，就设置为1,但只能有一个为1*/
#define TEST_FUN_CHOICE		"Griewank"	
//改变测试函数，请同时将此宏定义改为对应测试函数名
#define Best_fitness		0.0					//改变测试函数，请同时更改此函数的全局最优解
#define left_range			-600.00			//每一维度坐标范围(如果是整数，请务必写为**.0，比如整数600写为600.0)
#define right_range			600.00			//每一维度坐标范围(如果是整数，请务必写为**.0，比如整数600写为600.0)

//注意：
//二维函数只能测试二维
//其余普通函数可测试任意维度
//旋转平移函数只能测试10维，30维和50维

//二维函数
#define Goldstein_and_Price				0		//range = 2.0,		Best_fitness = 3
#define Martin_and_Gaddy				0		//left_range = 0.0, right_range = 10.0,		Best_fitness = 0
#define six_hump_camel_back_function	0		//range 3.0,		Best_fitness = -1.031628
//单峰函数
#define Sphere							0		//range = 100.0,	Best_fitness = 0
#define Sphere_Shifted					0		//range = 100.0,	Best_fitness = -450.0,	fbias = -450.0	平移函数
#define Axis_parallel_hyper_ellipsoid	0		//range = 5.12,		Best_fitness = 0
//多峰函数
#define Rastrigin						0		//range = 5.12,		Best_fitness = 0
#define Rastrigin_Shifted				1		//range = 5.0,		Best_fitness = -330.0,	fbias = -330.0	平移函数
#define Rastrigin_Shifted_Rotated		0		//range = 5.0,		Best_fitness = -330.0,	fbias = -330.0	旋转平移函数
#define Ackley							0		//range = 32.0,		Best_fitness = 0
#define Ackley_Shifted_Rotated			0		//range = 32.0,		Best_fitness = -140.0,	fbias = -140.0	旋转平移函数
#define Griewank						1		//range = 600.0,	Best_fitness = 0
#define Schwefel						0		//range = 500.0,	Best_fitness = 0
#define Generalized_Penalized_1			0		//range = 50.0,		Best_fitness = 0
#define Generalized_Penalized_2			0		//range = 50.0,		Best_fitness = 0

//非正式测试函数
#define Demo							0
#define Rosenbrock						0		//range = 100.0,	Best_fitness = 0
class SHEEP
{
public:
	double	coordinate[DIM];					//存储每只羊的坐标
	double	fitness;							//存储适度值
	int		number;								//存储羊的编号
	int		scatter;							//判断这只羊是否被打散 1被打散，0没被打散
public:
	friend class GROUPSHEEP;
};
class GROUPSHEEP
{
private:										//牧场范围，也就是解空间范围
	double	left;								//存储每维坐标的范围
	double	right;
public:
	SHEEP	sheep[SNUM];						//羊群
	int		bellwethernumber;					//领头羊编号
	double	worstfitness;						//种群最差适度值
	double	meanfitness;						//整个羊群的适度值的平均值
	double	oldbellwetherfitness;				//存储上一代领头羊的适度值
	int		generationTimes;
public:
	GROUPSHEEP();
	void	initofgroup();						//初始化种群
	void	leader();							//领头羊阶段
	void	wander();							//羊群漫游阶段
	void	bellwether();						//更新一次领头羊,并且求出羊群的平均适度值和最差适度值
	int		huntaway();							//牧羊犬阶段
};

double Computafitness(double a[]);				//测试函数

#endif
