#include<iostream>
#include<random>
#include<cmath>
#include<vector>
#include<list>
#include<iomanip>
#include<iterator>
#include<algorithm>
#include<fstream>
#include<string>
#include"AnnealingAlgorithm.h"
#include"EasyBMP.h"
#include"EasyBMP_Geometry.h"
#include "EasyBMP_Font.h"

using namespace std;
uniform_real_distribution<> dis(0,1);

/*random function related*/
void Random_init();
float Random_realVal();//0~1, float
int Random_intVal(int min, int max);//min~max, int
unsigned int seed_val;
myclock::time_point beginning = myclock::now();
myclock::duration d;
mt19937 Rand_gen;//keep one global instance


AnnealMethod::AnnealMethod()
{
	//Random_init();
}
void Random_init()
{
	d = myclock::now() - beginning;
	seed_val = d.count();
	Rand_gen.seed(seed_val);
}
/*return 0~1, float*/
float Random_realVal()
{
	return dis(Rand_gen);
}
/*return min~max, int*/
int Random_intVal(int min, int max)
{
	float rd;
	rd = dis(Rand_gen)*(max - min + 1) + min;// +1 to make it possible for the function to get max val.
	return (int)rd;
}
//clear and rebuild data after begin , if begin = 0, all data will be cleared and rebuilded 
float AnnealMethod::Anneal_cost_Cal(int begin)
{
	//cost function: COST = AREA * (punish coefficient)// + coe * WIRE_LENGTH;
	int seq_id2, seq_id1;//seq id start from 1~2N-1
	int i = begin;
	struct node tempnode;
	tempnode.node_style = 0;
	for (i; i < (Expression.N - 1); i++)
	{
		Expression.OperatorSeq.at(i).ComCorner.clear();//clear all member in the list and rebuild all data after that
		tempnode.SelCornerPts1 = 0;//CLEAR TEMPNODE AFTER EVERY NODE CALCULATION
		tempnode.SelCornerPts2 = 0;
		tempnode.SelNode1 = 0;
		tempnode.SelNode2 = 0;
		tempnode.Operand1 = NULL;
		tempnode.Operand2 = NULL;
		tempnode.Operator1 = NULL;
		tempnode.Operator2 = NULL;
		tempnode.next = NULL;
		tempnode.check_flag = 0;
		tempnode.head_flag = 0;
		tempnode.oper_style = 0;
		seq_id2 = Expression.OperatorSeq.at(i).Oper_Seq_ID - 1;//seq_id2 is the pos of operator need to be calculated in polish expression
		int operatorflag = Expression.OperatorSeq.at(i).OperatorType;
//		int module_index2;
		tempnode.oper_style = operatorflag;
		if (Expression.Polish_Express.at(seq_id2 - 1)>0)//means #, not '*,+'
		{
			seq_id1 = seq_id2 - 1;
			if (Expression.Polish_Express.at(seq_id1 - 1)>0)
			{
				tempnode.node_style = 1;
				tempnode.Operand2 = &Expression.OperandSeq.at(seq_id2 - i - 1);
				tempnode.Operand1 = &Expression.OperandSeq.at(seq_id1 - i - 1);
				for (int j = 0; j < tempnode.Operand1->CornerPts.size(); j++)
				{
					tempnode.SelCornerPts1 = j;
					for (int k = 0; k < tempnode.Operand2->CornerPts.size(); k++)
					{
						tempnode.SelCornerPts2 = k;
						if (operatorflag == -1)
						{
							tempnode.w = tempnode.Operand1->CornerPts.at(j).w + tempnode.Operand2->CornerPts.at(k).w;
							tempnode.h = max(tempnode.Operand1->CornerPts.at(j).h, tempnode.Operand2->CornerPts.at(k).h);
						}
						else
						{
							tempnode.w = max(tempnode.Operand1->CornerPts.at(j).w, tempnode.Operand2->CornerPts.at(k).w);
							tempnode.h = tempnode.Operand1->CornerPts.at(j).h + tempnode.Operand2->CornerPts.at(k).h;
						}
						Expression.OperatorSeq.at(i).ComCorner.push_back(tempnode);
					}
				}
			}
			else
			{
				tempnode.node_style = 4;
				tempnode.Operand1 = &Expression.OperandSeq.at(seq_id2 - i - 1);
				tempnode.Operator1 = &Expression.OperatorSeq.at(i - 1);
				list<struct node>::iterator it;
				for (int j = 0; j < tempnode.Operand1->CornerPts.size(); j++)
				{
					it = tempnode.Operator1->ComCorner.begin();
					tempnode.SelCornerPts1 = j;
					for (int k = 0; k < tempnode.Operator1->ComCorner.size(); k++)
					{
						tempnode.SelNode1 = k;
						if (operatorflag == -1)
						{
							tempnode.w = tempnode.Operand1->CornerPts.at(j).w + it->w;
							tempnode.h = max(tempnode.Operand1->CornerPts.at(j).h, it->h);
						}
						else
						{
							tempnode.w = max(tempnode.Operand1->CornerPts.at(j).w, it->w);
							tempnode.h = tempnode.Operand1->CornerPts.at(j).h + it->h;
						}
						Expression.OperatorSeq.at(i).ComCorner.push_back(tempnode);
						it++;
					}
				}
			}
		}
		else
		{
			list<struct node>::iterator it;
			tempnode.node_style = 2;
			tempnode.Operator1 = &Expression.OperatorSeq.at(i - 1);
			it = tempnode.Operator1->ComCorner.begin();
			int seq_delt=0;
			int flag=1;
			while (flag)
			switch (it->node_style)
			{
			case 1:seq_id1 = it->Operand1->Oper_Seq_ID - 2;//seq_delt += 4; 
				flag = 0;//get out of loop
				break;
			case 2:seq_id1 = it->Operand1->Oper_Seq_ID - 2;//pos in polish expression
				flag = 0;//flag for getting out of loop
				break;
			case 3: int temp; 
				it = it->Operator1->ComCorner.begin();
				break;
			case 4:it = it->Operator1->ComCorner.begin();
				break;//break;
			}
			if (Expression.Polish_Express.at(seq_id1)>0)
			{
				tempnode.node_style = 2;
				int cn;
				for (cn = 0; cn < (Expression.N - 1); cn++)
				{
					if (Expression.OperandSeq.at(cn).OperandID == Expression.Polish_Express.at(seq_id1))
						break;
				}
				tempnode.Operand1 = &Expression.OperandSeq.at(cn);
				tempnode.Operator1 = &Expression.OperatorSeq.at(i - 1);
				list<struct node>::iterator itt;
				for (int j = 0; j < tempnode.Operand1->CornerPts.size(); j++)
				{
					itt = tempnode.Operator1->ComCorner.begin();
					tempnode.SelCornerPts1 = j;
					for (int k = 0; k < tempnode.Operator1->ComCorner.size(); k++)
					{
						tempnode.SelNode1 = k;
						if (operatorflag == -1)
						{
							tempnode.w = tempnode.Operand1->CornerPts.at(j).w + itt->w;
							tempnode.h = max(tempnode.Operand1->CornerPts.at(j).h, itt->h);
						}
						else
						{
							tempnode.w = max(tempnode.Operand1->CornerPts.at(j).w, itt->w);
							tempnode.h = tempnode.Operand1->CornerPts.at(j).h + itt->h;
						}
						Expression.OperatorSeq.at(i).ComCorner.push_back(tempnode);
						itt++;
					}
				}
			}
			else
			{
				tempnode.node_style = 3;
				int cn;
				for (cn = 0; cn < (Expression.N - 2); cn++)
				{
					if (Expression.OperatorSeq.at(cn).Oper_Seq_ID == (seq_id1+1))
						break;
				}
				tempnode.Operator1 = &Expression.OperatorSeq.at(cn);
				tempnode.Operator2 = &Expression.OperatorSeq.at(i - 1);
				list<struct node>::iterator it1 = tempnode.Operator1->ComCorner.begin();
				list<struct node>::iterator it2;
				for (int j = 0; j < tempnode.Operator1->ComCorner.size(); j++)
				{
					it2 = tempnode.Operator2->ComCorner.begin();
					tempnode.SelNode1 = j;
					for (int k = 0; k < tempnode.Operator2->ComCorner.size(); k++)
					{
						tempnode.SelNode2 = k;
						if (operatorflag == -1)
						{
							tempnode.w = it1->w + it2->w;
							tempnode.h = max(it1->h, it2->h);
						}
						else
						{
							tempnode.w = max(it1->w, it2->w);
							tempnode.h = it1->h + it2->h;
						}
						Expression.OperatorSeq.at(i).ComCorner.push_back(tempnode);
						it2++;
					}
					it1++;
				}
			}
		}
		CornerPtsMerge(i);
	}
	/*area cal and decide which val should be return*/
	list<struct node>::iterator checker, node_smallest_area, node_smallest_area_fit_pq;
	checker = node_smallest_area_fit_pq = node_smallest_area = Expression.OperatorSeq.back().ComCorner.begin();
	float aspect_ratio = checker->h / checker->w;
	float area = checker->h*checker->w;
	float area_fit_pq = area;
	for (int x = 1; x < Expression.OperatorSeq.back().ComCorner.size(); x++)
	{
		checker++;
		float area_comp = checker->h*checker->w;
		aspect_ratio = checker->h / checker->w;
		if (area_comp < area)
		{
			node_smallest_area = checker;
			area = area_comp;
			if ((aspect_ratio >= Expression.p) && (aspect_ratio <= Expression.q))
			{
				node_smallest_area_fit_pq = checker;
				area_fit_pq = area_comp;
			}
		}
	}
	aspect_ratio = node_smallest_area_fit_pq->h / node_smallest_area_fit_pq->w;
	Expression.Sel_node_final = node_smallest_area;
	/*return logic: if ratio fits the module constrain exists, then return area of that node, otherwise return area with lowest value with punished coefficient*/
	if ((aspect_ratio >= Expression.p) && (aspect_ratio <= Expression.q))
	{
		return area_fit_pq;
	}
	else
	{
		return (punish_coefficient*area);//1.5 is used to give the program a trend to area_fit_pq
	}
}
void AnnealMethod::CornerPtsMerge(int i)
{
	list<struct node>::iterator it1, it2, itend, itout;
	float h1, w1, h2, w2;
	it1 = Expression.OperatorSeq.at(i).ComCorner.begin();
	itout = itend = Expression.OperatorSeq.at(i).ComCorner.end();
	itend--;
	while (it1 != itend)
	{ 
		int breakflag = 0;
		h1 = it1->h; w1 = it1->w;
		it2 = it1;
		it2++;//start at next data
		while(it2!=itout)
		{
			h2 = it2->h; w2 = it2->w;
			if ((h1 >= h2)&&(w1 >= w2))
			{
				it1 = Expression.OperatorSeq.at(i).ComCorner.erase(it1);
				breakflag = 1;
				break;
			}
			else
			{
				if ((h1 <= h2) && (w1 <= w2))
				{
					it2 = Expression.OperatorSeq.at(i).ComCorner.erase(it2);
					continue;//it2++ is not performed here 'cause the it2 will auto point to next data .
				}
			}
			it2++;
		}
		if (breakflag == 0)
		{
			it1++;
			if (it1 == itout)
				break;
			itend = Expression.OperatorSeq.at(i).ComCorner.end();
			itend--;
		}
	}
}

bool AnnealMethod::Anneal_move1()
{
	int rand = Random_intVal(0, Expression.N-2);//0~43, swap 1~44;
	int id1 = Expression.OperandSeq.at(rand).Oper_Seq_ID;
	int id2 = Expression.OperandSeq.at(rand + 1).Oper_Seq_ID;
	Expression.OperandSeq.at(rand).Oper_Seq_ID = id2;
	Expression.OperandSeq.at(rand + 1).Oper_Seq_ID = id1;
	vector<Def_Operand>::iterator first = Expression.OperandSeq.begin();
	vector<Def_Operand>::iterator second = Expression.OperandSeq.begin();
	iter_swap(first + rand, second + rand + 1);//swap
	Anneal_expression_updata(true);//updata complete polish expression
	return true;
}
bool AnnealMethod::Anneal_move2()
{
	int rand = Random_intVal(0, Expression.Operator_chain.size() - 1);//0~size-1
	list<Def_Operator_chain_begin>::iterator it = Expression.Operator_chain.begin();
	advance(it, rand);
	int p = it->Operator_Serial;
	int pos = Expression.OperatorSeq.at(p).Oper_Seq_ID-1;
	while (pos < (2*Expression.N-1))
	{
		if (Expression.Polish_Express.at(pos) == -1)
		{
			Expression.OperatorSeq.at(p).OperatorType = -2;
			p++;
		}
		else
		{
			if (Expression.Polish_Express.at(pos) == -2)
			{
				Expression.OperatorSeq.at(p).OperatorType = -1;
				p++;
			}
			else
				break;//come with 1,2,3...operands
		}
		pos++;
	}
	Anneal_expression_updata(true);
	return true;
}
bool AnnealMethod::Anneal_move3()
{
	bool flag = false;
	int falsecnt = 0;
	while ((flag != true) && falsecnt < (10 * Expression.N))
	{	
		int randchoice = Random_intVal(1, 2);//1,2
		int randpos;
		list<Def_Operator_chain_begin>::iterator it = Expression.Operator_chain.begin();//todo : check if the programm here is suitable
		switch (randchoice)
		{
		case 1://operand & operator
			randpos = Random_intVal(0, Expression.Operator_chain.size() - 1);
			advance(it, randpos);
			if (Expression.Polish_Express.at(it->Oper_Seq_ID - 2) != Expression.Polish_Express.at(it->Oper_Seq_ID))
			{
				int dk = 2 * (it->Operator_Serial + 1);
				if (dk < it->Oper_Seq_ID)
				{
					int id1, id2;
					id1 = Expression.OperandSeq.at(it->Operand_Serial).Oper_Seq_ID;
					id2 = Expression.OperatorSeq.at(it->Operator_Serial).Oper_Seq_ID;
					Expression.OperandSeq.at(it->Operand_Serial).Oper_Seq_ID = id2;
					Expression.OperatorSeq.at(it->Operator_Serial).Oper_Seq_ID = id1;
					flag = true;
				}
			}
			break;
		case 2://operator & operand
			if (Expression.Operator_chain.size() < 2)//todo
				break;
			randpos = Random_intVal(0, Expression.Operator_chain.size() - 2);
			advance(it, randpos);
			if (Expression.Polish_Express.at(it->Oper_Seq_ID_R) != Expression.Polish_Express.at(it->Oper_Seq_ID_R + 2))
			{
				int id1, id2;
				id1 = Expression.OperandSeq.at(it->Operand_Serial_R).Oper_Seq_ID;
				id2 = Expression.OperatorSeq.at(it->Operator_Serial_R).Oper_Seq_ID;
				Expression.OperandSeq.at(it->Operand_Serial_R).Oper_Seq_ID = id2;
				Expression.OperatorSeq.at(it->Operator_Serial_R).Oper_Seq_ID = id1;
				flag = true;
			}
			break;
		}
		falsecnt++;
	}
	if (flag == true)
	{
		Anneal_expression_updata(true);
	}
	return flag;
}

bool AnnealMethod::Anneal_move()
{
	//random produce move #1#2#3
	bool flag = false;
	while (flag==false)
	{
		int rand = Random_intVal(1, 3);
		switch (rand)
		{
		case 1:flag = Anneal_move1(); mov1++;
			break;
		case 2:flag = Anneal_move2(); mov2++;
			break;
		case 3:flag = Anneal_move3(); mov3++;
			break;
		}
	}
	return true;
}
void AnnealMethod::AnnealProcedure()
{
	float rej_cnt;
	float rej_ratio = 0;
	float cost, delt_cost;
	Def_ChartData ChartData;
	Def_Expression_ Expression_backup = Expression;//copy Expression for recovery after one move was rejected.
	while ((T > T_frozen)&&(rej_ratio < 0.99))
	{
		int i;
		rej_cnt = 0;
		ChartData.area = Expression.Sel_node_final->h*Expression.Sel_node_final->w;
		ChartData.t = T;
		Chart.push_back(ChartData);
		for (i = 0; i < K_move; i++)
		{
			Anneal_move();
			cost = Anneal_cost_Cal(0);//todo
			delt_cost = cost - Cost_best;
			float randreal;
			if ((delt_cost > 0) && ((randreal = Random_realVal()) > exp(-delt_cost / T)))
			{
				Expression = Expression_backup;//recovery to previous data after rejected
				rej_cnt++;
			}
			else
			{
				Expression_backup = Expression;
				Cost_best = cost;
			}
			if (delt_cost < 0)//todo, what if =0
			{
				Expression_best = Expression;//updata best expression to the current expression
			}
		}
		Anneal_cost_Cal(0);//updata all the relationship between nodes, must be pressented here
		rej_ratio = rej_cnt / K_move;
		T = Dec_ratio*T;
	}
	Expression = Expression_best;
	Anneal_cost_Cal(0);//updata all the relationship between nodes, must be pressented here
	printf("rej_ratio = %.2f\n", rej_ratio);
}

bool AnnealMethod::InputfileRead()
{
	tmp_ = 0;
	//temp_ = 1;
	printf("input the name of module info file.(include file name extension)\n");
	string filename;
	cin >> filename;
	ifstream netlist_file(filename);
	if (!netlist_file.is_open())
	{
		printf("file open failed!\n");
		return false;//error
	}
	string line;
	getline(netlist_file, line);
	if (Getheader(line) < 1)
		return false;
	StructCreat();
	while (getline(netlist_file, line))
	{
		if (line.size()<5)//less than 5 means the line contain 'enter'
			continue;
		Getbody(line);
		tmp_++;
	}
	netlist_file.close();
	return true;
}
float AnnealMethod::getdatanbr(string line, string::size_type pos)
{
	string::size_type tmp;
	string str;
	str.clear();
	for (tmp = pos; tmp < line.size(); tmp++)
	{
		if ((line[tmp] >= '0') && (line[tmp] <= '9') || (line[tmp] == '.'))
			str += line[tmp];
		else
			break;
	}
	return stof(str);
}
char AnnealMethod::Getheader(string line)
{
	string::size_type pos=0;
	Expression.N = (int)getdatanbr(line, pos);
	pos = line.find(" ", pos + 1);
	if (pos == string::npos)
	{
		printf("header error!\n");
		return -1;//error
	}
	Expression.p = getdatanbr(line, pos + 1);
	pos = line.find(" ", pos + 1);
	if (pos == string::npos)
	{
		printf("header error!\n");
		return -1;//error
	}
	Expression.q = getdatanbr(line, pos + 1);
	return 1;//succeed
}
char AnnealMethod::Getbody(string line)
{
	string::size_type pos = 0;
	int i = 0;
	Expression.OperandSeq.at(tmp_).OperandID = (int)getdatanbr(line, pos);
	pos = line.find(" ", pos + 1);
	if (pos == string::npos)
	{
		printf("body error!\n");
		return -1;//error
	}
	Expression.OperandSeq.at(tmp_).block_A = getdatanbr(line, pos + 1);
	pos = line.find(" ", pos + 1);
	if (pos == string::npos)
	{
		printf("body error!\n");
		return -1;//error
	}
	Expression.OperandSeq.at(tmp_).block_r = getdatanbr(line, pos + 1);
	pos = line.find(" ", pos + 1);
	if (pos == string::npos)
	{
		printf("body error!\n");
		return -1;//error
	}
	Expression.OperandSeq.at(tmp_).block_s = getdatanbr(line, pos + 1);
	pos = line.find(" ", pos + 1);
	if (pos == string::npos)
	{
		printf("body error!\n");
		return -1;//error
	}
	Expression.OperandSeq.at(tmp_).S = getdatanbr(line, pos + 1);
	return 1;//succeed
}
void AnnealMethod::StructCreat()
{
	Expression.OperandSeq.resize(Expression.N);
	Expression.OperatorSeq.resize(Expression.N - 1);
	Expression.Polish_Express.resize(2 * Expression.N - 1, 0);
}
void AnnealMethod::RefineStruct()
{
	int seq, i, operatype;
	float net_area_val = 0;
	seq = 2;
	operatype = -1;//todo
	Expression.OperandSeq.at(0).Oper_Seq_ID = 1;
	Expression.OperandSeq.at(0).Oper_Seq_ID = 1;
	Expression.OperandSeq.at(0).FinalCornerPts = 0;//init to 0
	Expression.OperandSeq.at(0).SelCornerPts = 0;//init to 0
	net_area_val += Expression.OperandSeq.at(0).block_A;
	CornerPtsCal(&Expression.OperandSeq.at(0));//get corner point for the module
	for (i = 1; i < Expression.N; i++)
	{
		Expression.OperandSeq.at(i).Oper_Seq_ID = seq++;
		Expression.OperandSeq.at(i).FinalCornerPts = 0;//init to 0
		Expression.OperandSeq.at(i).SelCornerPts = 0;//init to 0
		net_area_val += Expression.OperandSeq.at(i).block_A;
		CornerPtsCal(&Expression.OperandSeq.at(i));//get corner point for the module
		Expression.OperatorSeq.at(i - 1).Oper_Seq_ID = seq++;
		Expression.OperatorSeq.at(i - 1).OperatorType = operatype;//-1 means '*', init to be 12*3*4*5*6*...
		if (operatype == -1)//todo
			operatype = -2;
		else
			operatype = -1;
	}
//	uselesstest();//useless todo
	printf("Area of All module: %.2f\n", net_area_val);
	Anneal_expression_updata(true);
}
void AnnealMethod::Anneal_expression_updata(bool ctrl)//todo
{	//updata polish expression according to the swap of operands and operators during move
	//updata operatorchain data
	for (int i = 0; i < Expression.N; i++)
	{
		Expression.Polish_Express.at(Expression.OperandSeq.at(i).Oper_Seq_ID - 1) = Expression.OperandSeq.at(i).OperandID;
		if (i==(Expression.N - 1))
			continue;
		Expression.Polish_Express.at(Expression.OperatorSeq.at(i).Oper_Seq_ID - 1) = Expression.OperatorSeq.at(i).OperatorType;
	}
	if (ctrl == false)
		return;
	Expression.Operator_chain.clear();//clear the original data
	Def_Operator_chain_begin Operatorchain;
	vector<int>::iterator it = Expression.Polish_Express.begin();
	int cnt = 0, operatorcnt = 0, operandcnt = 0;
	while (it != Expression.Polish_Express.end())
	{
		while (*it > 0)
		{
			it++;
			cnt++;
			operandcnt++;
		}
		Operatorchain.Oper_Seq_ID = cnt;
		Operatorchain.Operand_Serial = operandcnt-1;
		Operatorchain.Operator_Serial = operatorcnt;
		Operatorchain.Oper_Seq_ID_R = -1;
		Operatorchain.Operand_Serial_R = -1;
		Operatorchain.Operator_Serial_R = -1;
		Expression.Operator_chain.push_back(Operatorchain);
		while (*it < 0)
		{
			it++;
			cnt++;
			operatorcnt++;
			if (it == Expression.Polish_Express.end())
				break;
		}
		if (it != Expression.Polish_Express.end())
		{
			Expression.Operator_chain.back().Oper_Seq_ID_R = cnt - 1;
			Expression.Operator_chain.back().Operand_Serial_R = operandcnt;
			Expression.Operator_chain.back().Operator_Serial_R = operatorcnt - 1;
		}
	}
}
void AnnealMethod::CornerPtsCal(Def_Operand *operand)
{
	int i=0;
	float delt, cnt = 0;
	Def_CornerPts cornerpts;
	delt = sqrt(operand->block_A / operand->block_r) - sqrt(operand->block_A / operand->block_s);
	delt /= (Pts_lim-1);
	while (i < Pts_lim)
	{
		cornerpts.w = cnt*delt + sqrt(operand->block_A / operand->block_s);
		cornerpts.h = operand->block_A / cornerpts.w;
		cornerpts.cornerID = i;
		operand->CornerPts.push_back(cornerpts);
		cnt++;
		i++;
		if (operand->block_r == operand->block_s)
		{
			break;
		}
	}
	if (operand->S == 1)
		return;
	delt = sqrt(operand->block_A * operand->block_s) - sqrt(operand->block_A * operand->block_r);
	delt /= (Pts_lim - 1);
	cnt = 0;
	while (i < 2 * Pts_lim)
	{
		cornerpts.w = cnt*delt + sqrt(operand->block_A * operand->block_r);
		cornerpts.h = operand->block_A / cornerpts.w;
		cornerpts.cornerID = i;
		operand->CornerPts.push_back(cornerpts);
		cnt++;
		i++;
		if (operand->block_r == operand->block_s)
		{
			break;
		}
	}
}
void AnnealMethod::Anneal_leaves_search(BMP *image)
{	
	float coe = 100;//coefficient for amplifing the plot
	float x1, x2, y2, y1;
	RGBApixel color;//black
	color.Blue = 0;
	color.Green = 0;
	color.Red = 0;
	list<struct node>::iterator it, it_tmp;
	Expression.Sel_node_final->head_flag = 1;//set the head node before starting the searching
	Expression.Sel_node_final->topright[0] = Expression.Sel_node_final->w;
	Expression.Sel_node_final->topright[1] = Expression.Sel_node_final->h;
	int opertype = Expression.Sel_node_final->oper_style;//get the operator type of the head node
	it_tmp = it = Expression.Sel_node_final;
	int nodesel, leafsel;
	while (1)
	{
		switch (it->node_style)
		{
		case 1://Operand1 & Operand2
			opertype = it->oper_style;
			it->Operand1->FinalCornerPts = it->SelCornerPts1;
			it->Operand2->FinalCornerPts = it->SelCornerPts2;

			it->Operand2->topright[0] = it->topright[0];
			it->Operand2->topright[1] = it->topright[1];
			if (opertype == -1)
			{
				it->Operand1->topright[0] = it->topright[0] - it->Operand2->CornerPts.at(it->SelCornerPts2).w;
				it->Operand1->topright[1] = it->topright[1];
				DrawFastLine(*image, it->Operand1->topright[0] * coe, it->topright[1] * coe, it->Operand1->topright[0] * coe, (it->topright[1] - it->h) * coe, color);
				it->Operand2->center[0] = it->topright[0] - it->Operand2->CornerPts.at(it->SelCornerPts2).w / 2;
				it->Operand2->center[1] = it->topright[1] - it->h / 2;
				it->Operand1->center[0] = it->Operand1->topright[0] - it->Operand1->CornerPts.at(it->SelCornerPts1).w / 2;
				it->Operand1->center[1] = it->topright[1] - it->h / 2;
			}
			else
			{
				it->Operand1->topright[1] = it->topright[1] - it->Operand2->CornerPts.at(it->SelCornerPts2).h;
				it->Operand1->topright[0] = it->topright[0];
				DrawFastLine(*image, it->topright[0] * coe, it->Operand1->topright[1] * coe, (it->topright[0] - it->w) * coe, it->Operand1->topright[1] * coe, color);
				it->Operand2->center[1] = it->topright[1] - it->Operand2->CornerPts.at(it->SelCornerPts2).h / 2;
				it->Operand2->center[0] = it->topright[0] - it->w / 2;
				it->Operand1->center[1] = it->Operand1->topright[1] - it->Operand1->CornerPts.at(it->SelCornerPts1).h / 2;
				it->Operand1->center[0] = it->topright[0] - it->w / 2;
			}
			it->check_flag = 0x11;
			if (it->head_flag == 1)//only happen when "12+" is the whole expression
				return;
			it = it->parentnode;//goback to its parent
			break;
		case 2://Operand1 & Operator1
			if (it->check_flag == 0x11)
			{
				if (it->head_flag != 1)
					it = it->parentnode;//go back to its parent node
				else
					return;//its already the head node then return.
			}
			else
			{
				opertype = it->oper_style;
				it->Operand1->FinalCornerPts = it->SelCornerPts1;
				it->check_flag = 0x11;
				it_tmp = it;
				it = it->Operator1->ComCorner.begin();
				advance(it, it_tmp->SelNode1);
				it->topright[0] = it_tmp->topright[0];
				it->topright[1] = it_tmp->topright[1];
				if (opertype == -1)
				{
					it_tmp->Operand1->topright[0] = it_tmp->topright[0] - it->w;
					it_tmp->Operand1->topright[1] = it_tmp->topright[1];
					DrawFastLine(*image, it_tmp->Operand1->topright[0] * coe, it_tmp->topright[1] * coe, it_tmp->Operand1->topright[0] * coe, (it_tmp->topright[1] - it_tmp->h) * coe, color);
					it_tmp->Operand1->center[0] = it_tmp->Operand1->topright[0] - it_tmp->Operand1->CornerPts.at(it_tmp->SelCornerPts1).w / 2;
					it_tmp->Operand1->center[1] = it_tmp->topright[1] - it_tmp->h / 2;
				}
				else
				{
					it_tmp->Operand1->topright[1] = it_tmp->topright[1] - it->h;
					it_tmp->Operand1->topright[0] = it_tmp->topright[0];
					DrawFastLine(*image, it_tmp->topright[0] * coe, it_tmp->Operand1->topright[1] * coe, (it_tmp->topright[0] - it_tmp->w) * coe, it_tmp->Operand1->topright[1] * coe, color);
					it_tmp->Operand1->center[1] = it_tmp->Operand1->topright[1] - it_tmp->Operand1->CornerPts.at(it_tmp->SelCornerPts1).h / 2;
					it_tmp->Operand1->center[0] = it_tmp->topright[0] - it_tmp->w / 2;
				}
				it->parentnode = it_tmp;
			}
			break;
		case 3://Operator1 & Operator2
			if (it->check_flag == 0x11)
			{
				if (it->head_flag != 1)
					it = it->parentnode;//go back to its parent node
				else
					return;//its already the head node then return.
			}
			else
			{		
				if (it->check_flag==0x00)
				{
					list<struct node>::iterator it_temp;
					it->check_flag = 0x01;
					opertype = it->oper_style;
					it_tmp = it;
					it_temp = it_tmp->Operator1->ComCorner.begin();
					advance(it_temp, it_tmp->SelNode1);
					it = it->Operator2->ComCorner.begin();
					advance(it, it_tmp->SelNode2);
					it->topright[0] = it_tmp->topright[0];
					it->topright[1] = it_tmp->topright[1];
					if (opertype == -1)
					{
						DrawFastLine(*image, (it_tmp->topright[0] - it->w) * coe, it_tmp->topright[1] * coe, (it_tmp->topright[0] - it->w) * coe, (it_tmp->topright[1] - it_tmp->h) * coe, color);
					}
					else
					{
						DrawFastLine(*image, it_tmp->topright[0] * coe, (it_tmp->topright[1] - it->h) * coe, (it_tmp->topright[0] - it_tmp->w) * coe, (it_tmp->topright[1] - it->h) * coe, color);
					}
					it->parentnode = it_tmp;
				}
				else
				{
					opertype = it->oper_style;
					list<struct node>::iterator it_temp;
					it->check_flag = 0x11;
					it_tmp = it;
					it = it->Operator1->ComCorner.begin();
					advance(it, it_tmp->SelNode1);
					it_temp = it_tmp->Operator2->ComCorner.begin();
					advance(it_temp, it_tmp->SelNode2);
					if (opertype == -1)
					{
						it->topright[0] = it_tmp->topright[0] - it_temp->w;
						it->topright[1] = it_tmp->topright[1];
				//		DrawFastLine(*image, it->topright[0] * coe, it_tmp->topright[1] * coe, it->topright[0] * coe, (it_tmp->topright[1] - it_tmp->h) * coe, color);
					}
					else
					{
						it->topright[1] = it_tmp->topright[1] - it_temp->h;
						it->topright[0] = it_tmp->topright[0];
					//	DrawFastLine(*image, it_tmp->topright[0] * coe, it->topright[1] * coe, (it_tmp->topright[0] - it_tmp->w) * coe, it->topright[1] * coe, color);
					}
					it->parentnode = it_tmp;
				}
			}
			break;
		case 4://Operator1 & Operand1
			if (it->check_flag == 0x11)
			{
				if (it->head_flag != 1)
					it = it->parentnode;//go back to its parent node
				else
					return;//its already the head node then return.
			}
			else
			{
				opertype = it->oper_style;
				it->check_flag = 0x11;
				nodesel = it->SelNode1;
				leafsel = it->SelCornerPts1;
				it->Operand1->FinalCornerPts = leafsel;
				it_tmp = it;
				it = it->Operator1->ComCorner.begin();
				advance(it, nodesel);
				it_tmp->Operand1->topright[0] = it_tmp->topright[0];
				it_tmp->Operand1->topright[1] = it_tmp->topright[1];
				if (opertype == -1)
				{
					it->topright[0] = it_tmp->topright[0] - it_tmp->Operand1->CornerPts.at(it_tmp->SelCornerPts1).w;
					it->topright[1] = it_tmp->topright[1];
					DrawFastLine(*image, it->topright[0] * coe, it_tmp->topright[1] * coe, it->topright[0] * coe, (it_tmp->topright[1] - it_tmp->h) * coe, color);
					it_tmp->Operand1->center[0] = it_tmp->topright[0] - it_tmp->Operand1->CornerPts.at(it_tmp->SelCornerPts1).w / 2;
					it_tmp->Operand1->center[1] = it_tmp->topright[1] - it_tmp->h / 2;
				}
				else
				{
					it->topright[0] = it_tmp->topright[0];
					it->topright[1] = it_tmp->topright[1] - it_tmp->Operand1->CornerPts.at(it_tmp->SelCornerPts1).h;
					DrawFastLine(*image, it_tmp->topright[0] * coe, it->topright[1] * coe, (it_tmp->topright[0] - it_tmp->w) * coe, it->topright[1] * coe, color);
					it_tmp->Operand1->center[1] = it_tmp->topright[1] - it_tmp->Operand1->CornerPts.at(it_tmp->SelCornerPts1).h / 2;
					it_tmp->Operand1->center[0] = it_tmp->topright[0] - it_tmp->w / 2;
				}
				it->parentnode = it_tmp;
			}
			break;//break;
		}
	}
}

void AnnealMethod::Anneal_init()
{
	float  cost;
	punish_coefficient = 1.5;
//	int Plot_N = 0;
	T_frozen = 0.5f;//todo: 
	Cost_best = Anneal_cost_Cal(0);//init best cost
//	FloorPLot(&Plot_N);//init plot
//	Rawdataprint();
	Expression_best = Expression;//record original expression
	Anneal_move2();//init a move to cal the init temp
	cost = Anneal_cost_Cal(0);//T = 2000;
	T = (-abs(Cost_best - cost)) / log(0.97);//end temperature init/start temperature init, actually, we could calcuate the start temperature via p=1; delt E=ln();
	printf("Initialize Temperature: %.2f\n", T);
	printf("\n");
	printf("Number of corner pts between r and s: %d\n", Pts_lim);
	printf("\n");
	printf("Initial polish expression: \n");
	printExpression();
	printf("\n");
	if (cost < Cost_best)
	{
		Expression_best = Expression;
		Cost_best = cost;
	}
	else
	{
		Expression = Expression_best;
		Anneal_cost_Cal(0);
	}
	K_move = Expression.N * 20;//2;//
	printf("Initialize cost:%.2f\n", Expression.Sel_node_final->h*Expression.Sel_node_final->w);
	printf("\n");
	Dec_ratio = 0.88f;
}
void AnnealMethod::Doing()
{
	int Plot_N = 1;
	Random_init();
	InputfileRead();
	RefineStruct();
	Anneal_init();
	FloorPLot(&Plot_N);//init plot
	Plot_N++;
	AnnealProcedure();
	FloorPLot(&Plot_N);//final plot
	plot_chart();//temperature vs cost
	printf("Final Polish expression: \n");
	printExpression();
	Rawdataprint();
	printf("\n");
	printf("Final cost is: %.2f\n", Cost_best = Expression.Sel_node_final->h*Expression.Sel_node_final->w);
	printf("\n");
	printf("Final Aspect ratio: %.2f\n", Expression.Sel_node_final->h / Expression.Sel_node_final->w);
	printf("\n");
	printf("Temperature decrease rate: %.2f\n", Dec_ratio);
	printf("\n");
	printf("Number of moves for every Temperature: %d\n", K_move);
	printf("\n");
	printf("mov1: %d, mov2: %d, mov3: %d\n", mov1, mov2, mov3);
	printf("\n");
	printf("Final Temperature: %.2f\n", T);
	printf("\n");
	printf("COMPLETE!\n");
}
void AnnealMethod::printExpression()
{
	for (int i = 0; i < Expression.Polish_Express.size(); i++)
	{
		if (Expression.Polish_Express.at(i)>0)
			printf("%d ", Expression.Polish_Express.at(i));
		else
		{
			if (Expression.Polish_Express.at(i) == -1)
				cout << "* ";
			else
				cout << "+ ";

		}
	}
	cout << "\n";
}
void AnnealMethod::FloorPLot(int *NUM)
{
	BMP image;
	RGBApixel color;//black
	float coe = 100;//coefficient for amplifing the plot
	color.Blue = 0;
	color.Green = 0;
	color.Red = 0;
	image.SetSize(Expression.Sel_node_final->w * coe+100, Expression.Sel_node_final->h * coe+100);
	DrawFastLine(image, 0, 0, 0, Expression.Sel_node_final->h * coe, color);
	DrawFastLine(image, 0, 0, Expression.Sel_node_final->w * coe, 0, color);
	DrawFastLine(image, Expression.Sel_node_final->w * coe, Expression.Sel_node_final->h * coe, 0, Expression.Sel_node_final->h * coe, color);
	DrawFastLine(image, Expression.Sel_node_final->w * coe, Expression.Sel_node_final->h * coe, Expression.Sel_node_final->w * coe, 0, color);

	Anneal_leaves_search(&image);
	plot_node(&image);
	plot_module(&image);
	char filename[32];
	sprintf_s(filename, "FloorPlot_%d.bmp", *NUM);
	image.WriteToFile(filename);
}

void AnnealMethod::plot_module(BMP *image)
{
	char buf[32];
	RGBApixel color;//black
	float coe = 100.f;//coefficient for amplifing the plot
	color.Blue = 0;
	color.Green = 0;
	color.Red = 0;
	RGBApixel color2,color3;//black
	color2.Blue = 0;
	color2.Green = 0;
	color2.Red = 255;
	color3.Blue = 255;
	color3.Green = 0;
	color3.Red = 50;
	float x0, x1, y0, y1;
	int finalsel;
	for (int i = 0; i < Expression.N; i++)
	{
		finalsel = Expression.OperandSeq.at(i).FinalCornerPts;
		y0 = (coe * (Expression.OperandSeq.at(i).center[1] - Expression.OperandSeq.at(i).CornerPts[finalsel].h / 2));
		y1 = (coe * (Expression.OperandSeq.at(i).center[1] + Expression.OperandSeq.at(i).CornerPts[finalsel].h / 2));
		x0 = (coe * (Expression.OperandSeq.at(i).center[0] - Expression.OperandSeq.at(i).CornerPts[finalsel].w / 2));
		x1 = (coe * (Expression.OperandSeq.at(i).center[0] + Expression.OperandSeq.at(i).CornerPts[finalsel].w / 2));
		DrawFastLine(*image, x0, y0, x1, y0, color);
		DrawFastLine(*image, x0, y0, x0, y1, color);
		DrawFastLine(*image, x1, y1, x0, y1, color);
		DrawFastLine(*image, x1, y1, x1, y0, color);
		DrawFastLine(*image, x1, y1, x0, y0, color2);
		DrawFastLine(*image, x1, y0, x0, y1, color2);
		sprintf_s(buf, "M%d\n (%.2f, %.2f)", Expression.OperandSeq.at(i).OperandID,(x0+x1)/(2*coe),(y1+y0)/(2*coe));
		PrintString(*image, buf, (x1+x0)/2-30, (y1+y0)/2-20, 12, color3);
	}
}
void AnnealMethod::plot_node(BMP *image)
{
	float coe = 100;//coefficient for amplifing the plot
	float x1, x2, y2, y1;
	RGBApixel color;//black
	color.Blue = 0;
	color.Green = 0;
	color.Red = 0;
	for (int i = 0; i < (Expression.N-1); i++)
	{
		list<struct node>::iterator it = Expression.OperatorSeq.at(i).ComCorner.begin();
		for (int j = 0; j < Expression.OperatorSeq.at(i).ComCorner.size(); j++)
		{
			if (it->next != NULL)
			{
				x2 = coe*it->topright[0];
				y2 = coe*it->topright[1];
				x1 = coe*(x2-it->w);
				y1 = coe*(y2-it->h);
				DrawFastLine(*image, x1, y2, x2, y2, color);
				DrawFastLine(*image, x2, y2, x2, y1, color);
				DrawFastLine(*image, x1, y1, x1, y2, color);
				DrawFastLine(*image, x1, y1, x2, y1, color);
				break;
			}
			it++;
		}
	}
}
void AnnealMethod::plot_chart()
{
	char buf[32];
	BMP image;
	RGBApixel color1,color2;//black
	float coe = 10, coa = 1;//coefficient for amplifing the plot
	color1.Blue = 0;
	color1.Green = 0;
	color1.Red = 0;
	color2.Blue = 0;
	color2.Green = 0;
	color2.Red = 255;
	image.SetSize(3000+100, 3000+100);
	float x0, x1, y0, y1;
	DrawFastLine(image, 20, 20, 20, 3000, color1);
	DrawFastLine(image, 20, 20, 3000, 20, color1);
	for (int i = 0; i < Chart.size()-1; i++)
	{
		y0 = coa*(Chart[i].area) + 20;
		x0 = Chart[i].t*coe+20;
		y1 = coa*(Chart[i + 1].area) + 20;
		x1 = Chart[i + 1].t*coe + 20;
		DrawFastLine(image, x0, y0, x1, y1, color2);
	}
	for (int i = 0; i < 30; i++)
	{
		sprintf_s(buf, "%d", i*100);
		PrintString(image, buf, 20, i*100+20, 10, color1);
	}
	sprintf_s(buf, "cost");
	PrintString(image, buf, 80, 200, 30, color2);
	for (int i = 0; i < 30; i++)
	{
		sprintf_s(buf, "%d", i*10 );
		PrintString(image, buf, i * 100 + 20, 20, 10, color1);
	}
	sprintf_s(buf, "%s", "T: C");
	PrintString(image, buf, 80, 80, 30, color2);
	image.WriteToFile("Area_VS_T_Plot.bmp");
}
void AnnealMethod::Rawdataprint()
{
	ofstream out;
		out.open("rawdata.txt");
		if (out.is_open())
	    {
			for (int i = 0; i < Expression.N; i++)
			{
				int sel = Expression.OperandSeq[i].FinalCornerPts;
				out << "(" << Expression.OperandSeq[i].OperandID << " , " << Expression.OperandSeq[i].block_A << " , " << Expression.OperandSeq[i].block_r << " , " << Expression.OperandSeq[i].block_s << " , " << Expression.OperandSeq[i].S <<")" <<"		"<< "(" << Expression.OperandSeq[i].CornerPts[sel].h << " , " << Expression.OperandSeq[i].CornerPts[sel].w << ")" << "\n";
			}
	        out.close();
	    }
		else
		{
			cout << "Error opening file";
		}
}

void AnnealMethod::uselesstest()
{}
