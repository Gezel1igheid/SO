#include <iostream>
#include "sheep.h"

FILE *fp;
FILE *fp2;

int main() {
	struct tm *newtime;
	char tmpbuf[128];
	char tmpbuf2[128];
	time_t lt1;
	lt1 = time(NULL);
	newtime = localtime(&lt1);
	srand((unsigned)time(NULL));

	sprintf(tmpbuf, "SHEEPv%s_U%-.2EDIM%d_SNUM%d_ITE%d_DH%.2f_TH%d_%d%d_%s_%d%02d%02d_%02d%02d%02d_File1.txt", SHEEP_VERSION,U,
		DIM, SNUM, ITE, DEGREE_HUNTAWAY, TIMES_HUNTAWAY, PRINT_BELLWETHER, PRINT_HUNTAWAY, TEST_FUN_CHOICE,
		1900 + newtime->tm_year, 1 + newtime->tm_mon, newtime->tm_mday, newtime->tm_hour, newtime->tm_min, newtime->tm_sec);
	sprintf(tmpbuf2, "SHEEPv%s_U%-.2EDIM%d_SNUM%d_ITE%d_DH%.2f_TH%d_%d%d_%s_%d%02d%02d_%d%02d%02d_File2.txt", SHEEP_VERSION,U,
		DIM, SNUM, ITE, DEGREE_HUNTAWAY, TIMES_HUNTAWAY, PRINT_BELLWETHER, PRINT_HUNTAWAY, TEST_FUN_CHOICE,
		1900 + newtime->tm_year, 1 + newtime->tm_mon, newtime->tm_mday, newtime->tm_hour, newtime->tm_min, newtime->tm_sec);

#if PRINT_BELLWETHER
	fp = fopen(tmpbuf, "wt");
	fprintf(fp,"%s\n\n", tmpbuf);
	fprintf(fp, "粒子维度:\t%d\n种群规模:\t%d\n迭代次数:\t%d\n左坐标范围:\t%f\n右坐标范围:\t%f\n每组测试数目:\t%d\n测试组数:\t%d\n牧羊犬介入程度:\t%.2f\n牧羊犬介入最小间隔:\t%d\n牧羊犬介入阀值:\t%lf\n测试函数名称:\t%s\n\n\n",
		DIM, SNUM, ITE, left_range, right_range, N, SET, DEGREE_HUNTAWAY, TIMES_HUNTAWAY, U, TEST_FUN_CHOICE);
#endif
#if PRINTF_RESULT_FILE
	fp2 = fopen(tmpbuf2, "wt");
	fprintf(fp2, "%s\n\n", tmpbuf2);
	fprintf(fp2, "粒子维度:\t%d\n种群规模:\t%d\n迭代次数:\t%d\n左坐标范围:\t%f\n右坐标范围:\t%f\n每组测试数目:\t%d\n测试组数:\t%d\n牧羊犬介入程度:\t%.2f\n牧羊犬介入最小间隔:\t%d\n牧羊犬介入阀值:\t%lf\n测试函数名称:\t%s\n\n\n",
		DIM, SNUM, ITE, left_range, right_range, N, SET, DEGREE_HUNTAWAY, TIMES_HUNTAWAY, U, TEST_FUN_CHOICE);//向文件中输出本次测试的备注信息
#endif
	for (int i = 0; i < SET; i++) {
#if PRINT_BELLWETHER
		fprintf(fp, "=======================================================================================================\n");
		fprintf(fp, "第%d组\n", i + 1);
#endif
#if PRINTF_RESULT_FILE
		fprintf(fp2, "======================================================================================================\n");
		fprintf(fp2, "第%d组\n", i + 1);
		fprintf(fp2, "测试次数\t达到误差允许阀值时的代数\t函数值\t函数误差\t最终编号\t最终函数值\t最终函数误差\n");
		fprintf(fp2, "注：达到误差允许阀值时的代数 列中数据为%5d时为未在最大迭代次数限制内找到误差允许的值。\n",ITE);
#endif
		long iteSUM = 0, bestITE = LONG_MAX, worstITE = LONG_MIN;
		clock_t start, finish;			//计算程序运行时间
		double sum = 0, best = LONG_MAX, worst = LONG_MIN, duration, satisfactFitness;
		double satisfactSum = 0, satisfactBest = LONG_MAX, satisfactWorst = LONG_MIN;
		start = clock();
		int k = 0;
		int isHuntaway, satisfactGnrtTm, satisfactTms = 0;
		for (k = 0; k < N; k++) {		//等于运行N 次程序 得到N个数据
			GROUPSHEEP group;
			isHuntaway = 0;
			satisfactGnrtTm = 0;
			group.initofgroup();		//初始化粒子群
			group.bellwether();			//找出领头羊
#if PRINT_BELLWETHER
			fprintf(fp, "代数\t编号\t函数值\t函数误差\t是否牧羊犬介入\t坐标\t\n");
			fprintf(fp, "%6d\t%3d\t%-E\t%-E\t", group.generationTimes, group.bellwethernumber, group.sheep[group.bellwethernumber - 1].fitness, group.sheep[group.bellwethernumber - 1].fitness - Best_fitness);
			fprintf(fp, "        \t");
			//for (int dim = 0;dim < DIM;dim++)
			//	fprintf(fp, "%-E\t", group.sheep[group.bellwethernumber - 1].coordinate[dim]);
			fprintf(fp, "\n");
#endif
			while (group.generationTimes < ITE) {
				group.leader();					//领头羊阶段
				group.wander();					//羊群漫游
				group.oldbellwetherfitness = group.sheep[group.bellwethernumber - 1].fitness;
				group.bellwether();				//找出领头羊
				isHuntaway = group.huntaway();	//牧羊犬阶段开始
				group.bellwether();				//重新更新领头羊
				//printf("牧羊犬阶段领头羊：编号%d    %.6lf\n", group.bellwethernumber, group.sheep[group.bellwethernumber - 1].fitness);
				//printf("worstfitness:%.6lf\nmeanfitness:%.6lf\n\n\n", group.worstfitness, group.meanfitness);
#if PRINT_BELLWETHER
				fprintf(fp, "%6d\t%3d\t%-E\t%-E\t", group.generationTimes, group.bellwethernumber, group.sheep[group.bellwethernumber - 1].fitness, group.sheep[group.bellwethernumber - 1].fitness - Best_fitness);
#if PRINT_HUNTAWAY
				if (isHuntaway == 1)
					fprintf(fp, "Huntaway");
#endif
				if (isHuntaway == 1)
					fprintf(fp, "\t");
				else 
					fprintf(fp, "        \t");
				//for (int dim = 0;dim < DIM;dim++)
					//fprintf(fp, "%-E\t", group.sheep[group.bellwethernumber - 1].coordinate[dim]);
				fprintf(fp, "\n");
#endif
				if (group.sheep[group.bellwethernumber - 1].fitness - Best_fitness < 1E-06 && satisfactGnrtTm == 0 || group.generationTimes == (ITE) && satisfactGnrtTm == 0) {
					iteSUM = iteSUM + (group.generationTimes);											//记录迭代次数的和
					if ((group.generationTimes) < bestITE) bestITE = group.generationTimes;				//记录最优迭代次数
					if ((group.generationTimes) > worstITE) worstITE = group.generationTimes;			//记录最优迭代次数
					satisfactGnrtTm = group.generationTimes;
					satisfactFitness = group.sheep[group.bellwethernumber - 1].fitness;
					if (satisfactGnrtTm < ITE) satisfactTms++;
				}
			}
			satisfactSum = satisfactFitness;			//计算20组数据的和
			if (satisfactFitness<satisfactBest)			//找出N组数据里的最好的值
				satisfactBest = satisfactFitness;
			if (satisfactFitness>satisfactWorst)		//找出N组数据里的最差值
				satisfactWorst = satisfactFitness;

			sum = sum + group.sheep[group.bellwethernumber - 1].fitness;	//计算20组数据的和
			if (group.sheep[group.bellwethernumber - 1].fitness<best)		//找出N组数据里的最好的值
				best = group.sheep[group.bellwethernumber - 1].fitness;
			if (group.sheep[group.bellwethernumber - 1].fitness>worst)		//找出N组数据里的最差值
				worst = group.sheep[group.bellwethernumber - 1].fitness;
#if PRINT_BELLWETHER
			if (satisfactGnrtTm == ITE) {
				fprintf(fp, "第%d次测试：\n在\t%5d\t代停止未能在最大迭代次数内找到误差允许的函数值：", k + 1,satisfactGnrtTm);
			}
			else
				fprintf(fp, "第%d次测试：\n在\t%5d\t代函数值达到误差允许：", k + 1, satisfactGnrtTm);
			fprintf(fp, "函数值：\t%E\t函数误差：\t%E\n", satisfactFitness, satisfactFitness - Best_fitness);
			fprintf(fp, "在\t%5d\t代达到终止条件：", group.generationTimes);
			fprintf(fp, "领头羊的编号\t%3d\t函数值：\t%E\t函数误差：\t%E\n",group.bellwethernumber, group.sheep[group.bellwethernumber - 1].fitness, group.sheep[group.bellwethernumber - 1].fitness - Best_fitness);
#endif
#if PRINTF_RESULT_FILE
			fprintf(fp2, "%d\t%5d\t%E\t%E\t", k + 1, satisfactGnrtTm, satisfactFitness, satisfactFitness - Best_fitness);
			fprintf(fp2, "%3d\t%E\t%E\n", group.bellwethernumber, group.sheep[group.bellwethernumber - 1].fitness, group.sheep[group.bellwethernumber - 1].fitness - Best_fitness);
#endif
			//printf("%d 编号%d    %.6lf\n", k + 1, group.bellwethernumber, group.sheep[group.bellwethernumber - 1].fitness);
			printf("Set %d/%d %d'/%d' is OK.\n", i + 1, SET, k + 1, N);
		}
#if PRINT_BELLWETHER
		if (satisfactTms == 0) {
			fprintf(fp, "\n未能在有限迭代次数内找到误差允许范围内的值。\n");
		}
		else{
			fprintf(fp, "\n在\t%d\t次测试中，有\t%d\t次测试找到了误差允许的值，本组测试数据综合统计如下：\n", N, satisfactTms);
			fprintf(fp, "最优函数值、函数误差：\t%E\t%E\n", satisfactBest, satisfactBest - Best_fitness);							//N组组组数据里的最好值存入文档
			fprintf(fp, "最差函数值、函数误差：\t%E\t%E\n", satisfactWorst, satisfactWorst - Best_fitness);						//N组组组数据里的最差值存入文档
			fprintf(fp, "平均函数值、函数误差：\t%E\t%E\n", satisfactSum / N, satisfactSum / N - Best_fitness);					//N组组组数据里的平均值存入文档
			fprintf(fp, "最小迭代次数：\t%6d\n", bestITE);				//N组组组数据里的最好迭代次数存入文档
			fprintf(fp, "最大迭代次数：\t%6d\n", worstITE);				//N组组组数据里的最差迭代次数存入文档
			fprintf(fp, "平均迭代次数：\t%.1lf\n", (double)iteSUM / k);	//N组组组数据里的平均迭代次数存入文档
		}
#endif
#if PRINTF_RESULT_FILE
		if (satisfactTms == 0) {
			fprintf(fp2, "\n未能在有限迭代次数内找到误差允许范围内的值。\n");
		}
		else {
			fprintf(fp2, "\n在\t%d\t次测试中，有\t%d\t次测试找到了误差允许的值，本组测试数据综合统计如下：\n", N, satisfactTms);
			fprintf(fp2, "最优函数值、函数误差：\t%E\t%E\n", satisfactBest, satisfactBest - Best_fitness);						//N组组组数据里的最好值存入文档
			fprintf(fp2, "最差函数值、函数误差：\t%E\t%E\n", satisfactWorst, satisfactWorst - Best_fitness);						//N组组组数据里的最差值存入文档
			fprintf(fp2, "平均函数值、函数误差：\t%E\t%E\n", satisfactSum / N, satisfactSum / N - Best_fitness);					//N组组组数据里的平均值存入文档
			fprintf(fp2, "最小迭代次数：\t%6d\n", bestITE);				//N组组组数据里的最好迭代次数存入文档
			fprintf(fp2, "最大迭代次数：\t%6d\n", worstITE);				//N组组组数据里的最差迭代次数存入文档
			fprintf(fp2, "平均迭代次数：\t%.1lf\n", (double)iteSUM / k);	//N组组组数据里的平均迭代次数存入文档
		}
#endif
		finish = clock();
		duration = (double)(finish - start) / CLOCKS_PER_SEC;
#if PRINT_BELLWETHER
		fprintf(fp, "\n每次测试都迭代到最大迭代次数的统计数据如下：\n");
		fprintf(fp, "最优函数值、函数误差：\t%E\t%E\n", best, best - Best_fitness);						//N组组组数据里的最好值存入文档
		fprintf(fp, "最差函数值、函数误差：\t%E\t%E\n", worst, worst - Best_fitness);						//N组组组数据里的最差值存入文档
		fprintf(fp, "平均函数值、函数误差：\t%E\t%E\n", sum / N, sum / N - Best_fitness);					//N组组组数据里的平均值存入文档
		fprintf(fp, "运行总时间：\t%lf秒\n", duration);
		fprintf(fp, "======================================================================================================\n\n\n\n");
#endif
#if PRINTF_RESULT_FILE
		fprintf(fp2, "\n每次测试都迭代到最大迭代次数的统计数据如下：\n");
		fprintf(fp2, "最优函数值、函数误差：\t%E\t%E\n", best, best - Best_fitness);						//N组组组数据里的最好值存入文档
		fprintf(fp2, "最差函数值、函数误差：\t%E\t%E\n", worst, worst - Best_fitness);					//N组组组数据里的最差值存入文档
		fprintf(fp2, "平均函数值、函数误差：\t%E\t%E\n", sum / N, sum / N - Best_fitness);				//N组组组数据里的平均值存入文档
		fprintf(fp2, "运行总时间：\t%lf秒\n", duration);
		fprintf(fp2, "======================================================================================================\n\n\n\n");
#endif
		//printf("best :%.6lf\n", best);							//N组组组数据里的最好值
		//printf("worst :%.6lf\n", worst);							//N组组组数据里的最差值
		printf("平均适应度: %E\n", sum / N - Best_fitness);			//N组组组数据里的平均值
		printf("Set %d is OK.\n", i + 1);
	}
	fclose(fp);
	fclose(fp2);
	cout << "Finish" << endl;
	return 0;
}
