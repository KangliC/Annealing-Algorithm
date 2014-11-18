#ifndef _ANNEALINGALGORITHM_H_
#define _ANNEALINGALGORITHM_H_
#include<iostream>
#include<random>
#include<chrono>
#include<vector>
#include<list>
#include<algorithm>
#include<string>
#include"EasyBMP.h"
#include"EasyBMP_Geometry.h"

using namespace std;
/*define parameters*/
#define Pts_lim	11

/*define data struct*/
typedef struct{
	float w;
	float h;
	int cornerID;//the # of corner ID among the bunch of corner points of one particular operand;
}Def_CornerPts;

typedef struct{
	int OperandID;//original ID from file
	int Oper_Seq_ID;//the ID in the expression sequence, 14*5+7*8*3*, seq_id of '4' is 2;
	float block_A;//area of the block
	float block_r;//r<=h/w in S1;
	float block_s;//h/w<=s in S1;
	int S;//S1,S2, S2 means the block can be rotated
	int SelCornerPts;//the selected corner point during the calculation
	int FinalCornerPts;//the final selected corner point after the one calculation of whole expression
	float topright[2];//top right coordinate(x, y)
	float center[2];
	vector<Def_CornerPts> CornerPts;//info of corner point
}Def_Operand;

typedef struct{
	int OperatorType;//'*':-1; '+':-2
	int Oper_Seq_ID;//The ID in the expression sequence, 14*5+7*8*3*, seq_ID of '+' is 5;
	list<struct node> ComCorner;
}Def_Operator;

struct node {
	float w;//combined corner pts
	float h;//combined corner pts
//	struct node *prev;
	struct node *next;
	list<struct node>::iterator parentnode;
	int oper_style;//record the operation of the node
	int node_style;//1,2,3,4: combination of different choices below
	/************************
	from left to right
	1: Operand1 & Operand2
	2: Operand1 & Operator1
	3: Operator1 & Operator2
	4: Operator1 & Operand1
	************************/
	char check_flag;//0x00,unchecked; 0x10, first checked; 0x01, second checked; 0x11, both checked;(first in here means first point, *operand or *operator(accordint to the node_style), second point means the left point)
	char head_flag;

	Def_Operand *Operand1;
	Def_Operand *Operand2;
	int SelCornerPts1;
	int SelCornerPts2;

	Def_Operator *Operator1;
	Def_Operator *Operator2;
	int SelNode1;//the # of selected node in the list of comcorner
	int SelNode2;
	float topright[2];//c1(x,y) top-right corner coordinate
};
typedef struct{
	int Oper_Seq_ID;//seq_id of operator which has operand on its left.
	int Operator_Serial;//pos in operatorseq
	int Operand_Serial;//pos in operandseq(operand on the left of the operator)

	int Oper_Seq_ID_R;//seq_id of operator which has operand on its RIGHT.
	int Operator_Serial_R;//pos in operatorseq
	int Operand_Serial_R;//pos in operandseq(operand on the right of the operator)
}Def_Operator_chain_begin;

typedef struct{
	vector<Def_Operand> OperandSeq;//sequence of Operand in the expression, should be N in size
	vector<Def_Operator> OperatorSeq;//sequence of operator, should be n-1 in size
//	list<Def_Operand_next_Operator> Operand_next_Operator;//record the seq_id of operator which has operand on its left.
	list<Def_Operator_chain_begin> Operator_chain;//store the operatorseq(in polish expression val) of begin operator of chain for M2 move and M3 move.
	vector<int> Polish_Express;//complete expression of problem
	list<struct node>::iterator Sel_node_final;//the final node # of the last operator

	int N;//total # of operands, e.g. # of module
	float p;//aspect ratio constrain of final chip(min)
	float q;//aspect ratio constrain of final chip(max)
}Def_Expression_;//size of a normalized polish expression is 2*n-1;

typedef struct{
	float t;//temperature
	float area;
	float p;//probability
}Def_ChartData;

typedef std::chrono::high_resolution_clock myclock;
class AnnealMethod{
public:
	AnnealMethod();
	void Doing();
private:
	/*file reading*/
	bool InputfileRead();//0:fail,1:succeed
	float getdatanbr(string line, string::size_type pos);
	char Getheader(string line);
	char Getbody(string line);
	void StructCreat();
	void RefineStruct();

	void CornerPtsCal(Def_Operand *operand);
	void CornerPtsMerge(int i);//i is operatorseq
	/*annealing process*/
	void AnnealProcedure();
	void Anneal_init();
	float Anneal_cost_Cal(int begin);//calculate start from begin
	void Anneal_leaves_search(BMP *image);
	bool Anneal_move();
	bool Anneal_move1();
	bool Anneal_move2();
	bool Anneal_move3();
	void Anneal_expression_updata(bool ctrl);//updata the complete expression, ctrl: false dont perform the Operatorchain updata, true perform Operatorchain updata.
	/*test useless*/
	void uselesstest();
	/*test useless*/

	void FloorPLot(int *NUM);
	void plot_module(BMP *image);
	void plot_node(BMP *image);
	void plot_chart();
	void printExpression();
	void Rawdataprint();

	float T_frozen;//end temperature
	float T;//temperature
	float Tp;//for calcuating initial T
	int K_move;//# of move in a particular temperature, Kp N
	int Kp;//
	float Cost_best;//the best cost during the procedure
	float Dec_ratio;//the rate of decrease the temperature
	float aspect_ratio_final;
	vector<Def_ChartData> Chart;//chart : area vs T
	int mov1 = 0, mov2 = 0, mov3 = 0;
//	list<struct node>::iterator Sel_node_final;//the final node # of the last operator
	
	int tmp_;
	int temp_;
	
	/*random function related*/	
	//void Random_init();
	//float Random_realVal();//0~1, float
	//int Random_intVal(int min, int max);//min~max, int
	//unsigned int seed_val;
	//myclock::time_point beginning = myclock::now();
	//myclock::duration d;
	//mt19937 Rand_gen;//keep one global instance
	float punish_coefficient;
	Def_Expression_ Expression;
	Def_Expression_ Expression_best;
};
void Random_init();
#endif