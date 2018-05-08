#pragma once 
#include "network.h"
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include"lib.h"

using namespace std;


//#define mu  1.0
//#define SF  1.0       //到時不用改為 1.0
//#define Pn  2.056*pow(10.0,(-11))     //mW 原先
//#define Pn  (double)4.11*pow(10.0,(-14))
#define Pn  (double) 4.11*pow(10.0,(-11.0))//mW
//#define Pn  1.0*pow(10.0,(-10.0))     //mW 根據後來論文使用的
//#define Keq 2.75*pow(10.0,(15))
#define par (int) 3//越大越準
#define nsum  3.0    //方法1(paper 切割法)***
//#define Alpha (float) 1.0
#define W  (int) 10.0 //M
//#define Gainthreshold 0.0 //***
//#define bsmaxpower 40.0      //mW ***
#define ga 3.0 // 方法1 ga可帶3~10 
//#define datarate 0.2 //0.25 //孝浩的Rmin
//#define PsCO 200.0 
//#define fa 2.3*1000
//#define DGA 100.0
//#define dm 1000.0 後來改掉了 by 11/22
//#define dp 500.0
//#define dcomp 250.0
//#define Pmax 30.0  //tier 1 Pmax
//#define Pmax2 15.0  //tier 2 Pmax
#define Pmax 5000.0  //tier 1 Pmax mW
#define Pmax2 100.0 //tier 2 Pmax mW
#define Pmax3 100.0 //tier 2 Pmax mW
#define Pmax4 100.0 //tier 2 Pmax mW
#define Pmax5 100.0 //tier 2 Pmax mW
#define Pmax6 100.0 //tier 2 Pmax mW
#define threshold 2 //comp門檻
#define Rmin 0.49 //Mbps
#define afa1 (double)-3.5 //path loss a1
#define afa2 (double)-3.5 //path loss a2
//----------------------------------------------------------------------------//



void NETWORK::AddNode_T(NODE *node) {
	if (Node_T == NULL) {
		Node_T = new LIST<NODE*>(node);
		Node_T->GetValue()->SetNo(N++);
		//cout<<"n="<< Node_T->GetValue()-> GetNo() <<endl;

	}
	else {
		Node_T->Tail()->SetNext(new LIST<NODE*>(Node_T->Tail(), node, NULL));
		Node_T->Tail()->GetValue()->SetNo(N++);
	}
}

void NETWORK::AddNode_B(NODE *node) {
	if (Node_B == NULL) {
		Node_B = new LIST<NODE*>(node);
		Node_B->GetValue()->SetNo(M++);

	}
	else {
		Node_B->Tail()->SetNext(new LIST<NODE*>(Node_B->Tail(), node, NULL));
		Node_B->Tail()->GetValue()->SetNo(M++);
	}
}




//-----------------------------------------------------------------------------------------------------//
void NETWORK::AddUBLB() { //NETWORK *network

	int M, N, n, n2, n3, n4;

	LIST<NODE*> *node, *node2;
	LIST<NODE*> *node3, *node4;


	M = this->GetNumberOfNode_B();
	N = this->GetNumberOfNode_T();

	//cout<<"UBLBHERE"<<endl;
	//system("pause");

	//Var----------------------------------------------------------------------初始化------------------------------------------//
	this->Var_LB = new float[M*N + 2 * N]; //Pij,Sj,Hj
	this->Var_UB = new float[M*N + 2 * N];

	for (int i = 0; i< M*N + 2 * N; i++) {
		Var_LB[i] = 0;
		Var_UB[i] = 0;
	}

	//Binary-----------------------------------------------------------------------------------

	this->BinVar_LB = new int[M*N*par]; //   λi,j,np
	this->BinVar_UB = new int[M*N*par];
	for (int i = 0; i< M*N*par; i++) {
		BinVar_LB[i] = 0;
		BinVar_UB[i] = 0;

	}

	//Bilinear-----------------------------------------------------------------------------------

	this->BilVar_LB = new float[2 * M*N*N + 2 * M*N*N*par]; // Zjij' Kji'j' ΔZjij'np ΔKji'j'np
	this->BilVar_UB = new float[2 * M*N*N + 2 * M*N*N*par];

	for (int i = 0; i< 2 * M*N*N + 2 * M*N*N*par; i++) {
		BilVar_LB[i] = 0;
		BilVar_UB[i] = 0;
	}

	//------------------------------------------------------------變數上下限------------------------------------------------//

	//double d;

	node = this->GetNode_B();
	while (node != NULL) {
		n = node->GetValue()->GetNo();//i

		node2 = this->GetNode_T();
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();//j

			if (n == 0) {
				Var_UB[n*N + n2] = Pmax;//tier 1 P
			}
			if (n == 1) {
				Var_UB[n*N + n2] = Pmax2;//tier 2 P
			}
			if (n == 2) {
				Var_UB[n*N + n2] = Pmax3;//tier 2 P
			}
			if (n == 3) {
				Var_UB[n*N + n2] = Pmax4;//tier 2 P
			}
			if (n == 4) {
				Var_UB[n*N + n2] = Pmax5;//tier 2 P
			}
			if (n == 5) {
				Var_UB[n*N + n2] = Pmax6 + 30;//tier 2 P
			}

			Var_UB[M*N + n2] = 10.0;//SINR 


			Var_UB[M*N + N + n2] = log(1 + Var_UB[M*N + n2]);//H   

			node2 = node2->GetNext();
		}
		node = node->GetNext();
	}

	//-----------------------------------------------Binvary----------------------------------------------------------
	node = this->GetNode_B();
	while (node != NULL) {
		n = node->GetValue()->GetNo();

		node2 = this->GetNode_T();
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();

			for (int j = 0; j<par; j++) {
				BinVar_UB[n*N*par + n2 * par + j] = 1;//λ
			}

			node2 = node2->GetNext();
		}
		node = node->GetNext();
	}

	//-------------------------------------------------Biliner Tearm--------------------------------------------------------------------
	//Z
	node = this->GetNode_B();
	while (node != NULL) {
		n = node->GetValue()->GetNo();

		node2 = this->GetNode_T();
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();

			node4 = this->GetNode_T();
			while (node4 != NULL) {
				n4 = node4->GetValue()->GetNo();

				if (n2 != n4) {

					BilVar_UB[n*N*N + n2 * N + n4] = (Var_UB[n*N + n2] * Var_UB[M*N + n4]);  //Z=P*S 
					for (int j = 0; j<par; j++) {
						BilVar_UB[2 * M*N*N + n * N*N*par + n2 * N*par + n4 * par + j] = Var_UB[M*N + n2]; //ΔZ:Sijn	
					}
				}
				node4 = node4->GetNext();
			}
			node2 = node2->GetNext();
		}
		node = node->GetNext();
	}
	//K
	node2 = this->GetNode_T();
	while (node2 != NULL) {
		n2 = node2->GetValue()->GetNo();

		node3 = this->GetNode_B();
		while (node3 != NULL) {
			n3 = node3->GetValue()->GetNo();

			node4 = this->GetNode_T();
			while (node4 != NULL) {
				n4 = node4->GetValue()->GetNo();

				if (n2 != n4) {

					BilVar_UB[M*N*N + n2 * M*N + n3 * N + n4] = (Var_UB[n3*N + n2] * Var_UB[M*N + n4]); //K

					for (int j = 0; j<par; j++) {
						BilVar_UB[2 * M*N*N + M * N*N*par + n2 * M*N*par + n3 * N*par + n4 * par + j] = Var_UB[M*N + n2]; //ΔK:Sijn 
					}
				}

				node4 = node4->GetNext();
			}

			node3 = node3->GetNext();
		}
		node2 = node2->GetNext();
	}

	return;
	system("pause");

}

//----------------------------Add vars 位置---------------------------------------------------------//

void NETWORK::func() { //NETWORK *network

					   //cout<<"Vars Location"<<endl;


	int M, N;

	M = this->GetNumberOfNode_B();
	N = this->GetNumberOfNode_T();

	ILOSTLBEGIN//#define ILOSTLBEGIN using namespace std;


		IloEnv env;//construct an Cplex environment env, which belongs to a handle class IloEnv;
	IloModel model(env);// define the varies "env"
	IloCplex cplex(model);
	IloNumVarArray vars(env), Vars(env), binvars(env), bilvars(env);
	IloNumArray Coes(env), vals(env), startVal(env);



	LIST<NODE*> *node, *node2;// TP
	LIST<NODE*> *node3, *node4;// BS
	LIST<NETWORK*> *network;

	int n, n2, n3, n4, i, j;


	try {

		//if(BNB_DDisplay){
		//	cout << "Adding vars."<<endl;


		char var_name[60];
		char cst_name[60];

		//==============================================================================
		//[0-1]取得每個距離

		double d;
		double distance[200][200];
		double dmin[200];
		//double Ijcounter[200];

		node2 = this->GetNode_T();
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();


			node = this->GetNode_B();
			while (node != NULL) {
				n = node->GetValue()->GetNo();


				d = sqrt(pow(node2->GetValue()->GetX1() - node->GetValue()->GetX2(), 2) + pow(node2->GetValue()->GetY1() - node->GetValue()->GetY2(), 2) + pow((node2->GetValue()->GetZ1()) - node->GetValue()->GetZ2(), 2));
				distance[n][n2] = d;

				node = node->GetNext();
			}
			node2 = node2->GetNext();
		}



		for (int j = 0; j<M; j++) {
			for (int i = 0; i<N; i++) {
				cout << "distance[" << j + 1 << "][" << i + 1 << "] =" << distance[j][i] << endl;
			}
		}

		//[0-2]取距離最近的pico
		node2 = this->GetNode_T();
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();

			double min = distance[1][n2];

			node = this->GetNode_B();
			while (node != NULL) {
				n = node->GetValue()->GetNo();

				if (n != 0) {
					if (distance[n][n2] <= min) {
						min = distance[n][n2];
						dmin[n2] = min;
					}
				}

				node = node->GetNext();
			}
			node2 = node2->GetNext();
		}


		//[0-3]  儲存gain矩陣
		double gain[200][200];
		double Hnorm;

		//double gainall[200];
		node = this->GetNode_B();
		while (node != NULL) {
			n = node->GetValue()->GetNo();

			node2 = this->GetNode_T();
			while (node2 != NULL) {
				n2 = node2->GetValue()->GetNo();

				if (n == 0)
				{
					//					Hnorm = (40*(log10(distance[n][n2]/1000.0))+30*(log10(2000.0))+49.0)/10.0;
					//					Hnorm = (137.4 + 35.2*(log10((distance[n][n2]/1000.0))))/10.0;//MARCO 
					Hnorm = pow(distance[n][n2], afa1);

				}

				if (n != 0)
				{
					//					Hnorm = (37+30*(log10(distance[n][n2]/1000.0))+12)/10.0;
					//					Hnorm = (137.4 + 35.2*(log10((distance[n][n2]/1000.0))))/10.0;//pico 
					Hnorm = pow(distance[n][n2], afa2);

				}

				//				gain[n][n2] = gg/(4.11*pow(10.0,(-14.0)));
				gain[n][n2] = Hnorm;
				//				gainall[n2] = gainall[n2] + gain[n][n2];
				cout << "gain[" << n + 1 << "][" << n2 + 1 << "]" << gain[n][n2] << endl;
				node2 = node2->GetNext();
			}
			node = node->GetNext();
		}

		//[0-4]確定Ij連線 確定好最近的pico然後看他的pint
		int Ij[200][200];
		double gainall[200];
		int comp[200];
		double pin[200];
		node2 = this->GetNode_T();
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();

			node = this->GetNode_B();
			while (node != NULL) {
				n = node->GetValue()->GetNo();

				if (distance[n][n2] == dmin[n2] && n != 0) {

					double pint = (5000.0*gain[0][n2]) / (100.0*gain[n][n2]);
					pin[n2] = pint;

					cout << "pint[" << n2 + 1 << "]=" << pint << endl;
					if (pint >= threshold) {//non-comp marco
						Ij[0][n2] = 1;
						gainall[n2] = gain[0][n2];
					}

					if (pint <= 1) {//non-comp pico
						Ij[n][n2] = 1;
						gainall[n2] = gain[n][n2];
					}

					if (pint <= threshold && pint > 1) {//comp
						Ij[0][n2] = 1;
						Ij[n][n2] = 1;
						gainall[n2] = gain[0][n2] + gain[n][n2];
						comp[n2] = 1;
					}
				}
				node = node->GetNext();
			}
			cout << "comp[" << n2 + 1 << "]=" << comp[n2] << endl;
			node2 = node2->GetNext();
		}

		for (int i = 0; i<N; i++) {
			for (int j = 0; j<M; j++) {
				cout << "Ij[" << j + 1 << "][" << i + 1 << "]=" << Ij[j][i] << endl;
			}
		}

		for (int j = 0; j<N; j++) {
			cout << "gainall[" << j + 1 << "]=" << gainall[j] << endl;
		}
		cout << "pint==============" << endl;
		for (int j = 0; j<N; j++) {
			cout << pin[j] << endl;
		}
		cout << "==================" << endl;


		cout << "gainall==============" << endl;
		for (int j = 0; j<N; j++) {
			cout << gainall[j] << endl;
		}
		cout << "==================" << endl;

		node2 = this->GetNode_T();
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();

			node = this->GetNode_B();
			while (node != NULL) {
				n = node->GetValue()->GetNo();

				if (Ij[n][n2] == 0) {
					Var_UB[n*N + n2] = 0;
				}

				node = node->GetNext();
			}
			node2 = node2->GetNext();
		}

		//============================================================================

		//Add Power
		node = this->GetNode_B();
		while (node != NULL) {
			n = node->GetValue()->GetNo();

			node2 = this->GetNode_T();
			while (node2 != NULL) {
				n2 = node2->GetValue()->GetNo();

				char var_name[60];
				sprintf(var_name, "P_D_%d_%d", n + 1, n2 + 1);
				vars.add(IloNumVar(env, Var_LB[n*N + n2], Var_UB[n*N + n2], var_name));


				node2 = node2->GetNext();
			}//n2 


			node = node->GetNext();
		}


		//Add SINR

		node2 = this->GetNode_T();
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();

			char var_name[60];
			sprintf(var_name, "S_D_%d", n2 + 1);
			vars.add(IloNumVar(env, Var_LB[M*N + n2], Var_UB[M*N + n2], var_name));

			node2 = node2->GetNext();
		}//n2 

		 //Add H_D = 1+SINR 

		node2 = this->GetNode_T();
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();

			char var_name[60];
			sprintf(var_name, "H_D_%d", n2 + 1);
			vars.add(IloNumVar(env, Var_LB[M*N + N + n2], Var_UB[M*N + N + n2], var_name));

			node2 = node2->GetNext();
		}//n2 


		 //------------------------------------------------------Biniary Varible----------------------------------------------------------------------//

		 //Addλ_D 	    //有改
		node = this->GetNode_B();
		while (node != NULL) {
			n = node->GetValue()->GetNo();

			node2 = this->GetNode_T();
			while (node2 != NULL) {
				n2 = node2->GetValue()->GetNo();

				for (int j = 0; j<par; j++) {

					char var_name[60];
					sprintf(var_name, "λ_D_%d_%d_%d ", n + 1, n2 + 1, j + 1);
					binvars.add(IloNumVar(env, BinVar_LB[n*N*par + n2 * par + j], BinVar_UB[n*N*par + n2 * par + j], ILOINT, var_name));//ILOINT代表此變數為整數

				}
				node2 = node2->GetNext();
			}

			node = node->GetNext();
		}

		//-----------------------------------------------------  Biliear Tearm------------------------------------------------------//

		//Add Z=P*S	

		node = this->GetNode_B();
		while (node != NULL) {
			n = node->GetValue()->GetNo();

			node2 = this->GetNode_T();
			while (node2 != NULL) {
				n2 = node2->GetValue()->GetNo();

				node4 = this->GetNode_T();
				while (node4 != NULL) {
					n4 = node4->GetValue()->GetNo();

					char var_name[60];
					sprintf(var_name, "Z_%d_%d_%d", n + 1, n2 + 1, n4 + 1);
					bilvars.add(IloNumVar(env, BilVar_LB[n*N*N + n2 * N + n4], BilVar_UB[n*N*N + n2 * N + n4], var_name));

					node4 = node4->GetNext();
				}//n4 

				node2 = node2->GetNext();
			}//n2 

			node = node->GetNext();
		}


		//Add K=P*S 


		node2 = this->GetNode_T();
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();

			node3 = this->GetNode_B();
			while (node3 != NULL) {
				n3 = node3->GetValue()->GetNo();

				node4 = this->GetNode_T();
				while (node4 != NULL) {
					n4 = node4->GetValue()->GetNo();

					char var_name[60];
					sprintf(var_name, "K_%d_%d_%d", n2 + 1, n3 + 1, n4 + 1);
					bilvars.add(IloNumVar(env, BilVar_LB[M*N*N + n2 * M*N + n3 * N + n4], BilVar_UB[M*N*N + n2 * M*N + n3 * N + n4], var_name));

					node4 = node4->GetNext();
				}//n4 

				node3 = node3->GetNext();
			}//n3 

			node2 = node2->GetNext();
		}//n2 



		 //Add ΔZ: Sijn
		node = this->GetNode_B();
		while (node != NULL) {
			n = node->GetValue()->GetNo();

			node2 = this->GetNode_T();
			while (node2 != NULL) {
				n2 = node2->GetValue()->GetNo();

				node4 = this->GetNode_T();
				while (node4 != NULL) {
					n4 = node4->GetValue()->GetNo();

					for (int j = 0; j<par; j++) {

						char var_name[60];
						sprintf(var_name, "△S_Z_%d_%d_%d_%d", n + 1, n2 + 1, n4 + 1, j + 1);
						bilvars.add(IloNumVar(env, BilVar_LB[2 * M*N*N + n * N*N*par + n2 * N*par + n4 * par + j], BilVar_UB[2 * M*N*N + n * N*N*par + n2 * N*par + n4 * par + j], var_name));

					}

					node4 = node4->GetNext();
				}//n4


				node2 = node2->GetNext();
			}//n2 


			node = node->GetNext();
		}

		//Add ΔK: Sijn  


		node2 = this->GetNode_T();
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();

			node3 = this->GetNode_B();
			while (node3 != NULL) {
				n3 = node3->GetValue()->GetNo();


				node4 = this->GetNode_T();
				while (node4 != NULL) {
					n4 = node4->GetValue()->GetNo();

					for (int j = 0; j<par; j++) {

						char var_name[60];
						sprintf(var_name, "△S_K_%d_%d_%d_%d", n2 + 1, n3 + 1, n4 + 1, j + 1);
						bilvars.add(IloNumVar(env, BilVar_LB[2 * M*N*N + M * N*N*par + n2 * M*N*par + n3 * N*par + n4 * par + j], BilVar_UB[2 * M*N*N + M * N*N*par + n2 * M*N*par + n3 * N*par + n4 * par + j], var_name));

					}

					node4 = node4->GetNext();
				}//n4 


				node3 = node3->GetNext();
			}//n3 


			node2 = node2->GetNext();
		}//n2 


		 //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^---變數宣告結束--------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

		 //-------------------------------------------------目標函式----------------------------
		 //Object function

		 //Cj
		for (int i = 0; i<N; i++) {

			Coes.add(W / log(2.0));
			Vars.add(vars[M*N + N + i]);//H     
										//			Coes.add( 1.0 );  
										//			Vars.add(vars[ M*N  + n2 ]);//H     

		}

		model.add(IloMaximize(env, IloScalProd(Coes, Vars)));

		Coes.clear();
		Vars.clear();



		//----------------------------------目標函式宣告結束-------------------------------
		//----------------------------------------------------------------------------------
		for (int j = 0; j<M; j++) {
			for (int i = 0; i<N; i++) {
				cout << "power[" << j + 1 << "][" << i + 1 << "]=" << vars[j*N + i] << endl;
			}
		}


		//start--------------------------------------------------------------
		//[1] Power 前面上下限已用過
		/*
		node2  = this-> GetNode_T();
		while (node2 != NULL){
		n2 = node2->GetValue()->GetNo();

		node  = this-> GetNode_B();
		while(node != NULL){
		n = node->GetValue()->GetNo();

		if(Ij[n][n2]==0){

		Coes.add(1.0);
		Vars.add(vars[ n*N + n2 ]);//P

		char cst_name[60];
		sprintf(cst_name, "P_%d_%d", n+1 , n2+1);
		IloConstraint cst;
		cst = IloScalProd(Coes,Vars) == 0 ;
		cst.setName(cst_name);
		model.add(cst);
		Coes.clear();
		Vars.clear();


		node  = node  ->GetNext();
		}
		node2  = node2  ->GetNext();
		}
		*/




		//[2-1]Power Sum tier 1

		node = this->GetNode_B();
		while (node != NULL) {
			n = node->GetValue()->GetNo();
			if (n == 0) {

				node2 = this->GetNode_T();
				while (node2 != NULL) {
					n2 = node2->GetValue()->GetNo();

					if (Ij[n][n2] == 1) {
						Coes.add(1.0);
						Vars.add(vars[n*N + n2]);//P
					}
					node2 = node2->GetNext();
				}

				char cst_name[60];
				sprintf(cst_name, "P_Sum_%d", n + 1);
				IloConstraint cst;
				cst = IloScalProd(Coes, Vars) <= 5000.0;
				cst.setName(cst_name);
				model.add(cst);
				Coes.clear();
				Vars.clear();
			}
			node = node->GetNext();
		}
		//[2-2]Power Sum tier 2

		node = this->GetNode_B();
		while (node != NULL) {
			n = node->GetValue()->GetNo();
			if (n != 0) {

				node2 = this->GetNode_T();
				while (node2 != NULL) {
					n2 = node2->GetValue()->GetNo();

					if (Ij[n][n2] == 1) {
						Coes.add(1.0);
						Vars.add(vars[n*N + n2]);//P
					}
					node2 = node2->GetNext();
				}

				char cst_name[60];
				sprintf(cst_name, "P_Sum_%d", n + 1);
				IloConstraint cst;
				cst = IloScalProd(Coes, Vars) <= 100.0;
				cst.setName(cst_name);
				model.add(cst);
				Coes.clear();
				Vars.clear();
			}
			node = node->GetNext();
		}




		//[3-1]SINR  non-comp Marco user
		node2 = this->GetNode_T();
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();

			node = this->GetNode_B();
			while (node != NULL) {
				n = node->GetValue()->GetNo();

				node4 = this->GetNode_T();
				while (node4 != NULL) {
					n4 = node4->GetValue()->GetNo();

					node3 = this->GetNode_B();
					while (node3 != NULL) {
						n3 = node3->GetValue()->GetNo();

						if (Ij[0][n2] == 1 && n == 0 && comp[n2] != 1) {
							if (n3 != n && n2 != n4 && comp[n4] != 1) {
								if (Ij[n3][n4] == 1) {

									Coes.add(gain[n3][n2]);//gi'j
									Vars.add(bilvars[M*N*N + n2 * M*N + n3 * N + n4]);//K
								}
							}
						}

						node3 = node3->GetNext();
					}

					if (Ij[0][n2] == 1 && n == 0 && comp[n2] != 1) {
						if (Ij[0][n4] == 1 && n2 != n4 && comp[n4] != 1) {
							if (gainall[n4] >= gainall[n2]) {

								Coes.add(gain[n][n2]);//gij
								Vars.add(bilvars[n*N*N + n2 * N + n4]);//Z

							}
						}
					}

					node4 = node4->GetNext();
				}

				if (Ij[0][n2] == 1 && n == 0 && comp[n2] != 1) {

					Coes.add(-gain[n][n2]);//gij
					Vars.add(vars[n*N + n2]);//P

				}

				node = node->GetNext();
			}

			if (Ij[0][n2] == 1 && comp[n2] != 1) {

				Coes.add(Pn);//N0
				Vars.add(vars[M*N + n2]);//SINR

										 //		cout << "SINR_D_%d" << n2+1 <<endl;
				sprintf(cst_name, "SINR marco_%d", n2 + 1);
				IloConstraint cst;
				cst = IloScalProd(Coes, Vars) == 0;
				cst.setName(cst_name);
				model.add(cst);
				Vars.clear();
				Coes.clear();
			}
			node2 = node2->GetNext();
		}

		//[3-2]SINR  non-comp pico user
		node2 = this->GetNode_T();
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();

			node = this->GetNode_B();
			while (node != NULL) {
				n = node->GetValue()->GetNo();

				node4 = this->GetNode_T();
				while (node4 != NULL) {
					n4 = node4->GetValue()->GetNo();

					node3 = this->GetNode_B();
					while (node3 != NULL) {
						n3 = node3->GetValue()->GetNo();
						//inter
						if (n != 0 && Ij[n][n2] == 1 && comp[n2] != 1) {
							if (Ij[0][n4] == 1 && n3 == 0 && Ij[n][n4] != 1) {
								if (n2 != n4) {

									Coes.add(gain[0][n2]);//gi'j
									Vars.add(bilvars[M*N*N + n2 * M*N + n3 * N + n4]);//K
								}
							}
						}

						if (n != 0 && Ij[n][n2] == 1 && comp[n2] != 1) {
							if (n3 != n && n3 != 0 && Ij[n3][n4] == 1) {
								if (n2 != n4) {
									Coes.add(gain[n3][n2]);//gi'j
									Vars.add(bilvars[M*N*N + n2 * M*N + n3 * N + n4]);//K
								}
							}
						}

						node3 = node3->GetNext();
					}
					//intra
					if (n != 0 && Ij[n][n2] == 1 && comp[n2] != 1) {
						if (n2 != n4 && Ij[n][n4] == 1 && comp[n4] != 1) {
							if (gainall[n4] >= gainall[n2]) {

								Coes.add(gain[n][n2]);//gij
								Vars.add(bilvars[n*N*N + n2 * N + n4]);//Z

							}
						}
					}

					node4 = node4->GetNext();
				}

				if (n != 0 && Ij[n][n2] == 1 && comp[n2] != 1) {

					Coes.add(-gain[n][n2]);//gij
					Vars.add(vars[n*N + n2]);//P

				}

				node = node->GetNext();
			}

			if (comp[n2] != 1 && Ij[0][n2] != 1) {

				Coes.add(Pn);//N0
				Vars.add(vars[M*N + n2]);//SINR

				sprintf(cst_name, "SINR pico_%d", n2 + 1);
				IloConstraint cst;
				cst = IloScalProd(Coes, Vars) == 0;
				cst.setName(cst_name);
				model.add(cst);
				Vars.clear();
				Coes.clear();
			}
			node2 = node2->GetNext();
		}

		//[3-3]SINR  comp user

		node2 = this->GetNode_T();
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();

			node = this->GetNode_B();
			while (node != NULL) {
				n = node->GetValue()->GetNo();

				node4 = this->GetNode_T();
				while (node4 != NULL) {
					n4 = node4->GetValue()->GetNo();

					node3 = this->GetNode_B();
					while (node3 != NULL) {
						n3 = node3->GetValue()->GetNo();
						//inter
						if (comp[n2] == 1 && comp[n4] == 1 && n == n3 && n2 != n4) {
							if (Ij[n][n2] != 1 && Ij[n][n4] == 1) {
								if (gainall[n4] >= gainall[n2]) {
									Coes.add(gain[n3][n2]);//gi'j
									Vars.add(bilvars[M*N*N + n2 * M*N + n3 * N + n4]);//K
								}
							}
						}

						if (comp[n2] == 1 && comp[n4] != 1 && n == n3 && n2 != n4) {
							if (Ij[n][n2] != 1 && Ij[n][n4] == 1) {
								Coes.add(gain[n3][n2]);//gi'j
								Vars.add(bilvars[M*N*N + n2 * M*N + n3 * N + n4]);//K
							}
						}

						//intra
						if (comp[n2] == 1 && comp[n4] == 1 && n == n3 && n2 != n4) {
							if (Ij[n][n2] == 1 && Ij[n][n4] == 1) {
								if (gainall[n4] >= gainall[n2]) {
									Coes.add(gain[n][n2]);//gij
									Vars.add(bilvars[n*N*N + n2 * N + n4]);//Z
								}
							}
						}

						//intra 
						if (comp[n2] == 1 && comp[n4] != 1 && n == n3 && n2 != n4) {
							if (Ij[n][n2] == 1 && Ij[n][n4] == 1) {
								Coes.add(gain[n][n2]);//gij
								Vars.add(bilvars[n*N*N + n2 * N + n4]);//Z
							}
						}


						node3 = node3->GetNext();
					}

					node4 = node4->GetNext();
				}

				if (comp[n2] == 1 && Ij[n][n2] == 1) {

					Coes.add(-gain[n][n2]);//gij
					Vars.add(vars[n*N + n2]);//Pij
				}

				node = node->GetNext();
			}
			if (comp[n2] == 1) {
				Coes.add(Pn);//N0
				Vars.add(vars[M*N + n2]);//SINR

				sprintf(cst_name, "SINR_%d", n2 + 1);
				IloConstraint cst;
				cst = IloScalProd(Coes, Vars) == 0;
				cst.setName(cst_name);
				model.add(cst);
				Vars.clear();
				Coes.clear();
			}
			node2 = node2->GetNext();
		}

		//[4-1]

		double cc = (1 - pow(2.0, Rmin / W));//C算式前面乘上的常數
		cout << "cc=" << cc << endl;
		node2 = this->GetNode_T();
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();

			node = this->GetNode_B();
			while (node != NULL) {
				n = node->GetValue()->GetNo();

				node4 = this->GetNode_T();
				while (node4 != NULL) {
					n4 = node4->GetValue()->GetNo();

					node3 = this->GetNode_B();
					while (node3 != NULL) {
						n3 = node3->GetValue()->GetNo();

						if (Ij[0][n2] == 1 && n == 0 && comp[n2] != 1) {
							if (n3 != n && n2 != n4 && comp[n4] != 1) {
								if (Ij[n3][n4] == 1) {
									//				cout<<"cc*gain[n3][n2]="<<cc*gain[n3][n2]<<endl;
									//				cout<<"vars[ n3*N + n4 "<<vars[ n3*N + n4 ]<<endl;
									//				Coes.add(cc*gain[n3][n2]);//cc*gi'j
									Coes.add(cc*gain[n3][n2]);//cc*gi'j
									Vars.add(vars[n3*N + n4]);//Pi'j'
								}
							}
						}
						node3 = node3->GetNext();
					}

					if (Ij[0][n2] == 1 && n == 0 && comp[n2] != 1) {
						if (Ij[0][n4] == 1 && n2 != n4 && comp[n4] != 1) {
							if (gainall[n4] >= gainall[n2]) {

								Coes.add(cc*gain[n][n2]);//cc*gij
								Vars.add(vars[n*N + n4]);//Pij'
							}
						}
					}
					node4 = node4->GetNext();
				}

				if (Ij[n][n2] == 1 && comp[n2] != 1 && n == 0) {

					Coes.add(gain[n][n2]);//gij
					Vars.add(vars[n*N + n2]);//Pij
				}
				node = node->GetNext();
			}

			if (Ij[0][n2] == 1 && comp[n2] != 1) {

				sprintf(cst_name, "C_Marco_%d", n2 + 1);
				IloConstraint cst;
				cst = IloScalProd(Coes, Vars) >= ((-cc)*Pn);
				//		cst = IloScalProd(Coes,Vars) >= -cc*Pn;
				cst.setName(cst_name);
				model.add(cst);
				Vars.clear();
				Coes.clear();
			}
			node2 = node2->GetNext();
		}

		//[4-2]
		node2 = this->GetNode_T();
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();

			node = this->GetNode_B();
			while (node != NULL) {
				n = node->GetValue()->GetNo();

				node4 = this->GetNode_T();
				while (node4 != NULL) {
					n4 = node4->GetValue()->GetNo();

					node3 = this->GetNode_B();
					while (node3 != NULL) {
						n3 = node3->GetValue()->GetNo();

						if (n != 0 && Ij[n][n2] == 1 && comp[n2] != 1) {
							if (Ij[0][n4] == 1 && n3 == 0 && Ij[n][n4] != 1) {
								if (n2 != n4) {

									Coes.add(cc*gain[n3][n2]);//cc*gi'j
									Vars.add(vars[n3*N + n4]);//Pi'j'
								}
							}
						}

						if (n != 0 && Ij[n][n2] == 1 && comp[n2] != 1) {
							if (n3 != n && n3 != 0 && Ij[n3][n4] == 1) {
								if (n2 != n4) {
									Coes.add(cc*gain[n3][n2]);//cc*gi'j
									Vars.add(vars[n3*N + n4]);//Pi'j'
								}
							}
						}

						node3 = node3->GetNext();
					}

					if (n != 0 && Ij[n][n2] == 1 && comp[n2] != 1) {
						if (n2 != n4 && Ij[n][n4] == 1 && comp[n4] != 1) {
							if (gainall[n4] >= gainall[n2]) {

								Coes.add(cc*gain[n][n2]);//cc*gij
								Vars.add(vars[n*N + n4]);//Pij'
							}
						}
					}
					node4 = node4->GetNext();
				}

				if (n != 0 && Ij[n][n2] == 1 && comp[n2] != 1) {

					Coes.add(gain[n][n2]);//gij
					Vars.add(vars[n*N + n2]);//Pij
				}
				node = node->GetNext();
			}
			if (Ij[0][n2] != 1 && comp[n2] != 1) {

				sprintf(cst_name, "C_non-comp pico_%d", n2 + 1);
				IloConstraint cst;
				cst = IloScalProd(Coes, Vars) >= ((-cc)*Pn);
				//		cst = IloScalProd(Coes,Vars) >= -cc*Pn;

				cst.setName(cst_name);
				model.add(cst);
				Vars.clear();
				Coes.clear();
			}
			node2 = node2->GetNext();
		}

		//[4-3]
		node2 = this->GetNode_T();
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();

			node = this->GetNode_B();
			while (node != NULL) {
				n = node->GetValue()->GetNo();

				node4 = this->GetNode_T();
				while (node4 != NULL) {
					n4 = node4->GetValue()->GetNo();

					node3 = this->GetNode_B();
					while (node3 != NULL) {
						n3 = node3->GetValue()->GetNo();
						//inter
						if (comp[n2] == 1 && comp[n4] == 1 && n == n3 && n2 != n4) {
							if (Ij[n][n2] != 1 && Ij[n][n4] == 1) {
								if (gainall[n4] >= gainall[n2]) {

									Coes.add(cc*gain[n3][n2]);//cc*gi'j
									Vars.add(vars[n3*N + n4]);//Pi'j'
								}
							}
						}

						if (comp[n2] == 1 && comp[n4] != 1 && n == n3 && n2 != n4) {
							if (Ij[n][n2] != 1 && Ij[n][n4] == 1) {

								Coes.add(cc*gain[n3][n2]);//cc*gi'j
								Vars.add(vars[n3*N + n4]);//Pi'j'
							}
						}
						//intra
						if (comp[n2] == 1 && comp[n4] == 1 && n == n3 && n2 != n4) {
							if (Ij[n][n2] == 1 && Ij[n][n4] == 1) {
								if (gainall[n4] >= gainall[n2]) {

									Coes.add(cc*gain[n][n2]);//cc*gij
									Vars.add(vars[n*N + n4]);//Pij'
								}
							}
						}
						//intra
						if (comp[n2] == 1 && comp[n4] != 1 && n == n3 && n2 != n4) {
							if (Ij[n][n2] == 1 && Ij[n][n4] == 1) {

								Coes.add(cc*gain[n][n2]);//cc*gij
								Vars.add(vars[n*N + n4]);//Pij'
							}
						}
						node3 = node3->GetNext();
					}

					node4 = node4->GetNext();
				}

				if (Ij[n][n2] == 1 && comp[n2] == 1) {

					Coes.add(gain[n][n2]);//gij
					Vars.add(vars[n*N + n2]);//Pij
				}

				node = node->GetNext();
			}

			if (comp[n2] == 1) {
				sprintf(cst_name, "C_comp_%d", n2 + 1);
				IloConstraint cst;
				cst = IloScalProd(Coes, Vars) >= ((-cc)*Pn);
				//		cst = IloScalProd(Coes,Vars) >= -cc*Pn;

				cst.setName(cst_name);
				model.add(cst);
				Vars.clear();
				Coes.clear();
			}
			node2 = node2->GetNext();
		}

		/*
		//[4]
		node2  = this-> GetNode_T();
		while (node2 != NULL){
		n2 = node2->GetValue()->GetNo();

		Coes.add(W/log(2.0));
		Vars.add(vars[ M*N + N  + n2 ]);
		cout<<"W/log(2.0)"<<W/log(2.0)<<"H="<<vars[ M*N + N  + n2 ]<<endl;
		char cst_name[60];
		sprintf(cst_name, "capacity_%d",n2+1  );
		IloConstraint cst;
		cst = IloScalProd(Coes,Vars) >= Rmin;
		cst.setName(cst_name);
		model.add(cst);
		Coes.clear();
		Vars.clear();
		node2 = node2-> GetNext();
		}
		*/

		//ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
		//[6.1]
		node = this->GetNode_B();
		while (node != NULL) {
			n = node->GetValue()->GetNo();

			node2 = this->GetNode_T();
			while (node2 != NULL) {
				n2 = node2->GetValue()->GetNo();
				for (int j = 0; j<par; j++) {
					Coes.add(1.0);
					Vars.add(binvars[n*N*par + n2 * par + j]);   //浪打
				}
				char cst_name[60];
				sprintf(cst_name, "Sumλ_D_%d.%d ", n + 1, n2 + 1);
				IloConstraint cst;
				cst = IloScalProd(Coes, Vars) == 1.0;
				cst.setName(cst_name);
				model.add(cst);
				Coes.clear();
				Vars.clear();
				//				cout<<"xxxxx"<<endl;

				node2 = node2->GetNext();
			}
			node = node->GetNext();
		}

		//[6.2.1]

		node = this->GetNode_B();
		while (node != NULL) {
			n = node->GetValue()->GetNo();

			node2 = this->GetNode_T();
			while (node2 != NULL) {
				n2 = node2->GetValue()->GetNo();

				node4 = this->GetNode_T();
				while (node4 != NULL) {
					n4 = node4->GetValue()->GetNo();

					Coes.add(1.0);
					Vars.add(vars[M*N + n2]); //SINRj			
					if (n2 != n4) {
						for (int j = 0; j<par; j++) {
							Coes.add(-1.0);
							Vars.add(bilvars[2 * M*N*N + n * N*N*par + n2 * N*par + n4 * par + j]);// Zsinr			 
						}
						char cst_name[60];
						sprintf(cst_name, "SINR_D_%d_%d_%d", n + 1, n2 + 1, n4 + 1);
						IloConstraint cst;
						cst = IloScalProd(Coes, Vars) == 0;
						cst.setName(cst_name);
						model.add(cst);
						Coes.clear();
						Vars.clear();
					}

					node4 = node4->GetNext();
				}//n4 

				node2 = node2->GetNext();
			}//n2 

			node = node->GetNext();
		}

		//[6.2.2]

		node = this->GetNode_B();
		while (node != NULL) {
			n = node->GetValue()->GetNo();

			node2 = this->GetNode_T();
			while (node2 != NULL) {
				n2 = node2->GetValue()->GetNo();

				double Smax = Var_UB[M*N + n2] - Var_LB[M*N + n2];
				node4 = this->GetNode_T();
				while (node4 != NULL) {
					n4 = node4->GetValue()->GetNo();

					if (n2 != n4) {
						for (int j = 0; j<par; j++) {

							Coes.add(-Smax);//SINR_U-SINR_L
							Vars.add(binvars[n*N*par + n4 * par + j]); //浪打

							Coes.add(1.0);
							Vars.add(bilvars[2 * M*N*N + n * N*N*par + n2 * N*par + n4 * par + j]);//Zsinr

							char cst_name[60];
							sprintf(cst_name, "△S_Z_%d_%d_%d_%d", n + 1, n2 + 1, n4 + 1, j + 1);
							IloConstraint cst;
							cst = IloScalProd(Coes, Vars) <= 0;
							cst.setName(cst_name);
							model.add(cst);
							Coes.clear();
							Vars.clear();
						}
					}
					node4 = node4->GetNext();
				}
				node2 = node2->GetNext();
			}
			node = node->GetNext();
		}
		//[6.3.1] RLTZ  
		//RLTZ1		

		node = this->GetNode_B();
		while (node != NULL) {
			n = node->GetValue()->GetNo();

			node2 = this->GetNode_T();
			while (node2 != NULL) {
				n2 = node2->GetValue()->GetNo();

				node4 = this->GetNode_T();
				while (node4 != NULL) {
					n4 = node4->GetValue()->GetNo();

					if (n2 != n4) {
						for (int j = 0; j<par; j++) {
							double a1 = (j + 1) / nsum;
							double a2 = j / nsum;
							double a = (pow(a1, ga) - pow(a2, ga))*(Var_UB[n*N + n4] - Var_LB[n*N + n4]);  //a np Pij'_U-Pij'_L

							Coes.add(-1.0*(Var_LB[n*N + n4] + a * j));
							Vars.add(bilvars[2 * M*N*N + n * N*N*par + n2 * N*par + n4 * par + j]); //Zsinr					
						}
						Coes.add(-1.0*Var_LB[M*N + n2]);  //SINRj_L			
						Vars.add(vars[n*N + n4]);//Pij'

						Coes.add(1.0);
						Vars.add(bilvars[n*N*N + n2 * N + n4]); //Z

						char cst_name[60];
						sprintf(cst_name, "Z:RLT1_%d.%d.%d", n + 1, n2 + 1, n4 + 1);
						IloConstraint cst;
						cst = IloScalProd(Coes, Vars) >= 0;
						cst.setName(cst_name);
						model.add(cst);
						Coes.clear();
						Vars.clear();
					}
					node4 = node4->GetNext();
				}//n4 	  
				node2 = node2->GetNext();
			}//n2 
			node = node->GetNext();
		}

		//[6.3.2]RLTZ2
		node = this->GetNode_B();
		while (node != NULL) {
			n = node->GetValue()->GetNo();

			node2 = this->GetNode_T();
			while (node2 != NULL) {
				n2 = node2->GetValue()->GetNo();


				node4 = this->GetNode_T();
				while (node4 != NULL) {
					n4 = node4->GetValue()->GetNo();

					if (n2 != n4) {
						for (int j = 0; j<par; j++) {
							double a1 = (j + 1) / nsum;
							double a2 = j / nsum;
							double a = (pow(a1, ga) - pow(a2, ga))*(Var_UB[n*N + n4] - Var_LB[n*N + n4]);

							Coes.add(-1.0*(Var_LB[n*N + n4] + a * (j + 1)));
							Vars.add(bilvars[2 * M*N*N + n * N*N*par + n2 * N*par + n4 * par + j]);  //Zsinr
						}
						Coes.add(-1.0*Var_LB[M*N + n2]);//SINRj_L
						Vars.add(vars[n*N + n4]);//Pij'

						Coes.add(1.0);
						Vars.add(bilvars[n*N*N + n2 * N + n4]); //Z

						char cst_name[60];
						sprintf(cst_name, "Z:RLT2_%d.%d.%d", n + 1, n2 + 1, n4 + 1);
						IloConstraint cst;
						cst = IloScalProd(Coes, Vars) <= 0;
						cst.setName(cst_name);
						model.add(cst);
						Coes.clear();
						Vars.clear();
					}
					node4 = node4->GetNext();
				}
				node2 = node2->GetNext();
			}
			node = node->GetNext();
		}

		//[6.3.3]RLTZ3
		node = this->GetNode_B();
		while (node != NULL) {
			n = node->GetValue()->GetNo();

			node2 = this->GetNode_T();
			while (node2 != NULL) {
				n2 = node2->GetValue()->GetNo();

				node4 = this->GetNode_T();
				while (node4 != NULL) {
					n4 = node4->GetValue()->GetNo();

					if (n2 != n4) {
						for (int j = 0; j<par; j++) {
							double a1 = (j + 1) / nsum;
							double a2 = j / nsum;
							double a = (pow(a1, ga) - pow(a2, ga))*(Var_UB[n*N + n4] - Var_LB[n*N + n4]);
							double Interval = Var_LB[n*N + n4] + a * j;

							Coes.add(-1.0*Interval);
							Vars.add(bilvars[2 * M*N*N + n * N*N*par + n2 * N*par + n4 * par + j]);//Zsinr

							Coes.add(Interval*(Var_UB[M*N + n2] - Var_LB[M*N + n2]));//Interval*(SINR_U-SINR_L)
							Vars.add(binvars[n*N*par + n4 * par + j]);//浪打
						}
						Coes.add(-1.0*Var_UB[M*N + n2]);    //SINR_U
						Vars.add(vars[n*N + n4]);//Pij'

						Coes.add(1.0);
						Vars.add(bilvars[n*N*N + n2 * N + n4]);//Z

						char cst_name[60];
						sprintf(cst_name, "Z:RLT3_%d.%d.%d", n + 1, n2 + 1, n4 + 1);
						IloConstraint cst;
						cst = IloScalProd(Coes, Vars) <= 0;
						cst.setName(cst_name);
						model.add(cst);
						Coes.clear();
						Vars.clear();
					}

					node4 = node4->GetNext();
				}//n4 
				node2 = node2->GetNext();
			}//n2 
			node = node->GetNext();
		}

		//[6.3.4]RLTZ4
		node = this->GetNode_B();
		while (node != NULL) {
			n = node->GetValue()->GetNo();

			node2 = this->GetNode_T();
			while (node2 != NULL) {
				n2 = node2->GetValue()->GetNo();

				node4 = this->GetNode_T();
				while (node4 != NULL) {
					n4 = node4->GetValue()->GetNo();

					if (n2 != n4) {
						for (int j = 0; j<par; j++) {
							double a1 = (j + 1) / nsum;
							double a2 = j / nsum;
							double a = (pow(a1, ga) - pow(a2, ga))*(Var_UB[n*N + n4] - Var_LB[n*N + n4]);

							double Interval = Var_LB[n*N + n4] + a * (j + 1);

							Coes.add(-1.0*Interval);
							Vars.add(bilvars[2 * M*N*N + n * N*N*par + n2 * N*par + n4 * par + j]);  //Zsinr

							Coes.add(Interval*(Var_UB[M*N + n2] - Var_LB[M*N + n2]));   //Interval*(SINRj_U-SINR_L)
							Vars.add(binvars[n*N*par + n4 * par + j]);//浪打
						}
						Coes.add(-1.0*Var_UB[M*N + n2]);//SINRj_U
						Vars.add(vars[n*N + n4]);//Pij'

						Coes.add(1.0);
						Vars.add(bilvars[n*N*N + n2 * N + n4]);//Z

						char cst_name[60];
						sprintf(cst_name, "Z:RLT4_%d.%d.%d", n + 1, n2 + 1, n4 + 1);
						IloConstraint cst;
						cst = IloScalProd(Coes, Vars) >= 0;
						cst.setName(cst_name);
						model.add(cst);
						Coes.clear();
						Vars.clear();
					}

					node4 = node4->GetNext();
				}
				node2 = node2->GetNext();
			}
			node = node->GetNext();
		}
		//-------------------------------------------------------Bilinear Tearm :K---------------------------------------//
		//[7.2.1]
		node2 = this->GetNode_T();
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();

			node3 = this->GetNode_B();
			while (node3 != NULL) {
				n3 = node3->GetValue()->GetNo();

				node4 = this->GetNode_T();
				while (node4 != NULL) {
					n4 = node4->GetValue()->GetNo();

					if (n2 != n4) {
						for (int j = 0; j<par; j++) {
							Coes.add(-1.0);
							Vars.add(bilvars[2 * M*N*N + M * N*N*par + n2 * M*N*par + n3 * N*par + n4 * par + j]);  //Ksinr					     	
						}
						Coes.add(1.0);
						Vars.add(vars[M*N + n2]);  //SINRj

						char cst_name[60];
						sprintf(cst_name, "SINR_D_%d_%d_%d", n2 + 1, n3 + 1, n4 + 1);
						IloConstraint cst;
						cst = IloScalProd(Coes, Vars) == 0;
						cst.setName(cst_name);
						model.add(cst);
						Coes.clear();
						Vars.clear();
					}
					node4 = node4->GetNext();
				}
				node3 = node3->GetNext();
			}
			node2 = node2->GetNext();
		}


		//[7.2.2] △S_K(n)
		node2 = this->GetNode_T();
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();

			node3 = this->GetNode_B();
			while (node3 != NULL) {
				n3 = node3->GetValue()->GetNo();

				node4 = this->GetNode_T();
				while (node4 != NULL) {
					n4 = node4->GetValue()->GetNo();

					double Smax = Var_UB[M*N + n2] - Var_LB[M*N + n2];  //有改   SINRj

					if (n2 != n4) {
						for (int j = 0; j<par; j++) {
							Coes.add(1.0);
							Vars.add(bilvars[2 * M*N*N + M * N*N*par + n2 * M*N*par + n3 * N*par + n4 * par + j]); //有改 Ksinr

							Coes.add(-Smax);
							Vars.add(binvars[n3*N*par + n4 * par + j]);  //浪打

							char cst_name[60];
							sprintf(cst_name, "△S_K_%d_%d_%d_%d", n2 + 1, n3 + 1, n4 + 1, j + 1);
							IloConstraint cst;
							cst = IloScalProd(Coes, Vars) <= 0;
							cst.setName(cst_name);
							model.add(cst);
							Coes.clear();
							Vars.clear();
						}
					}
					node4 = node4->GetNext();
				}
				node3 = node3->GetNext();
			}
			node2 = node2->GetNext();
		}

		//[7.3.1]RLTK
		//RLTK1		
		node2 = this->GetNode_T();
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();

			node3 = this->GetNode_B();
			while (node3 != NULL) {
				n3 = node3->GetValue()->GetNo();

				node4 = this->GetNode_T();
				while (node4 != NULL) {
					n4 = node4->GetValue()->GetNo();

					if (n2 != n4) {
						for (int j = 0; j<par; j++) {
							double a1 = (j + 1) / nsum;
							double a2 = (j) / nsum;
							double a = (pow(a1, ga) - pow(a2, ga))*(Var_UB[n3*N + n4] - Var_LB[n3*N + n4]);

							Coes.add(-1.0* (Var_LB[n3*N + n4] + a * j));
							Vars.add(bilvars[2 * M*N*N + M * N*N*par + n2 * M*N*par + n3 * N*par + n4 * par + j]);  //Ksinr
						}
						Coes.add(-1.0*Var_LB[M*N + n2]);//SINRj_L
						Vars.add(vars[n3*N + n4]);//Pi'j'

						Coes.add(1.0);
						Vars.add(bilvars[M*N*N + n2 * M*N + n3 * N + n4]);   //K

						char cst_name[60];
						sprintf(cst_name, "K:RLT1_.%d.%d.%d", n2 + 1, n3 + 1, n4 + 1);
						IloConstraint cst;
						cst = IloScalProd(Coes, Vars) >= 0;
						cst.setName(cst_name);
						model.add(cst);
						Coes.clear();
						Vars.clear();
					}
					node4 = node4->GetNext();
				}
				node3 = node3->GetNext();
			}
			node2 = node2->GetNext();
		}

		//[7.3.2] RLTK2
		node2 = this->GetNode_T();
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();

			node3 = this->GetNode_B();
			while (node3 != NULL) {
				n3 = node3->GetValue()->GetNo();

				node4 = this->GetNode_T();
				while (node4 != NULL) {
					n4 = node4->GetValue()->GetNo();

					if (n2 != n4) {
						for (int j = 0; j<par; j++) {
							double a1 = (j + 1) / nsum;
							double a2 = (j) / nsum;
							double a = (pow(a1, ga) - pow(a2, ga))*(Var_UB[n3*N + n4] - Var_LB[n3*N + n4]);//

							Coes.add(-1.0*(Var_LB[n3*N + n4] + a * (j + 1)));
							Vars.add(bilvars[2 * M*N*N + M * N*N*par + n2 * M*N*par + n3 * N*par + n4 * par + j]); //Ksinr
						}
						Coes.add(-1.0*Var_LB[M*N + n2]);  //SINRj_L
						Vars.add(vars[n3*N + n4]);//Pi'j'

						Coes.add(1.0);
						Vars.add(bilvars[M*N*N + n2 * M*N + n3 * N + n4]);   //K 

						char cst_name[60];
						sprintf(cst_name, "K:RLT2_.%d.%d.%d", n2 + 1, n3 + 1, n4 + 1);
						IloConstraint cst;
						cst = IloScalProd(Coes, Vars) <= 0;
						cst.setName(cst_name);
						model.add(cst);
						Coes.clear();
						Vars.clear();
					}
					node4 = node4->GetNext();
				}
				node3 = node3->GetNext();
			}
			node2 = node2->GetNext();
		}

		//[7.3.3]RLTK3			 
		node2 = this->GetNode_T();
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();

			node3 = this->GetNode_B();
			while (node3 != NULL) {
				n3 = node3->GetValue()->GetNo();

				node4 = this->GetNode_T();
				while (node4 != NULL) {
					n4 = node4->GetValue()->GetNo();

					if (n2 != n4) {
						for (int j = 0; j<par; j++) {
							double a1 = (j + 1) / nsum;
							double a2 = (j) / nsum;
							double a = (pow(a1, ga) - pow(a2, ga))*(Var_UB[n3*N + n4] - Var_LB[n3*N + n4]);//Pi'j' 

							double Interval = Var_LB[n3*N + n4] + a * j;

							Coes.add(-1.0*Interval);
							Vars.add(bilvars[2 * M*N*N + M * N*N*par + n2 * M*N*par + n3 * N*par + n4 * par + j]);  //Ksinr

							Coes.add(Interval*(Var_UB[M*N + n2] - Var_LB[M*N + n2]));   //SINRj_U-SINR_L
							Vars.add(binvars[n3*N*par + n4 * par + j]);//浪打
						}
						Coes.add(-1.0*Var_UB[M*N + n2]);//SINRj_U
						Vars.add(vars[n3*N + n4]);//Pi'j'

						Coes.add(1.0);
						Vars.add(bilvars[M*N*N + n2 * M*N + n3 * N + n4]);//K

						char cst_name[60];
						sprintf(cst_name, "K:RLT3_.%d.%d.%d", n2 + 1, n3 + 1, n4 + 1);
						IloConstraint cst;
						cst = IloScalProd(Coes, Vars) <= 0;
						cst.setName(cst_name);
						model.add(cst);
						Coes.clear();
						Vars.clear();
					}
					node4 = node4->GetNext();
				}
				node3 = node3->GetNext();
			}
			node2 = node2->GetNext();
		}

		//[7.3.4]RLTK4		 
		node2 = this->GetNode_T();
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();

			node3 = this->GetNode_B();
			while (node3 != NULL) {
				n3 = node3->GetValue()->GetNo();

				node4 = this->GetNode_T();
				while (node4 != NULL) {
					n4 = node4->GetValue()->GetNo();

					if (n2 != n4) {
						for (int j = 0; j<par; j++) {
							double a1 = (j + 1) / nsum;
							double a2 = (j) / nsum;
							double a = (pow(a1, ga) - pow(a2, ga))*(Var_UB[n3*N + n4] - Var_LB[n3*N + n4]);//方法1  

							double Interval = Var_LB[n3*N + n4] + a * (j + 1);

							Coes.add(-1.0*Interval);
							Vars.add(bilvars[2 * M*N*N + M * N*N*par + n2 * M*N*par + n3 * N*par + n4 * par + j]);   //有改Ksinr

							Coes.add(Interval*(Var_UB[M*N + n2] - Var_LB[M*N + n2]));   //有改SINRj
							Vars.add(binvars[n3*N*par + n4 * par + j]);   //浪打
						}

						Coes.add(-1.0*Var_UB[M*N + n2]);  //SINRj
						Vars.add(vars[n3*N + n4]);//Pi'j'

						Coes.add(1.0);
						Vars.add(bilvars[M*N*N + n2 * M*N + n3 * N + n4]);  //K

						char cst_name[60];
						sprintf(cst_name, "K:RLT4_.%d.%d.%d", n2 + 1, n3 + 1, n4 + 1);
						IloConstraint cst;
						cst = IloScalProd(Coes, Vars) >= 0;
						cst.setName(cst_name);
						model.add(cst);
						Coes.clear();
						Vars.clear();
					}

					node4 = node4->GetNext();
				}

				node3 = node3->GetNext();
			}

			node2 = node2->GetNext();
		}

		//[8]=========================================================
		double HUB, HLB, slmMinus, beta, R, D, PL;
		//[8.2]		 
		node2 = this->GetNode_T();
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();

			PL = 1.0 + Var_LB[M*N + n2];//1+SINRj_L
			R = log(PL) - 1.0;

			Coes.add(PL);
			Vars.add(vars[M*N + N + n2]);  //hd
			Coes.add(-1.0);
			Vars.add(vars[M*N + n2]);  //SINRj

			sprintf(cst_name, "Capacity1_D_.%d", n2 + 1);
			IloConstraint cst;
			cst = IloScalProd(Coes, Vars) <= PL * R + 1;
			cst.setName(cst_name);
			model.add(cst);
			Coes.clear();
			Vars.clear();

			node2 = node2->GetNext();
		}

		//[8.3]			 
		node2 = this->GetNode_T();
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();

			HLB = 1 + Var_LB[M*N + n2];//1+SINRj_L
			HUB = 1 + Var_UB[M*N + n2];//1+SINRj_U
			slmMinus = Var_UB[M*N + n2] - Var_LB[M*N + n2];  //SINRj_U-SINRj_L

			beta = ((HLB*HUB*(log(HUB) - log(HLB))) / slmMinus) - 1.0;
			PL = 1.0 + beta;
			R = log(PL) - 1.0;

			Coes.add(PL);
			Vars.add(vars[M*N + N + n2]);//H

			Coes.add(-1.0);
			Vars.add(vars[M*N + n2]);  //SINRj

			char cst_name[60];
			sprintf(cst_name, "Capacity2_D_%d", n2 + 1);
			IloConstraint cst;
			cst = IloScalProd(Coes, Vars) <= PL * R + 1;
			cst.setName(cst_name);
			model.add(cst);
			Coes.clear();
			Vars.clear();

			node2 = node2->GetNext();
		}

		//[8.4]			 
		node2 = this->GetNode_T();
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();

			PL = 1.0 + Var_UB[M*N + n2];
			R = log(PL) - 1.0;

			Coes.add(PL);
			Vars.add(vars[M*N + N + n2]);//H

			Coes.add(-1.0);
			Vars.add(vars[M*N + n2]);//SINRj

			char cst_name[60];
			sprintf(cst_name, "Capacity3_D_%d", n2 + 1);
			IloConstraint cst;
			cst = IloScalProd(Coes, Vars) <= PL * R + 1.0;
			cst.setName(cst_name);
			model.add(cst);
			Coes.clear();
			Vars.clear();

			node2 = node2->GetNext();
		}

		//[8.5]			 
		node2 = this->GetNode_T();
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();

			slmMinus = Var_UB[M*N + n2] - Var_LB[M*N + n2];  //SINRj_U-SINRj_L
			HUB = log(1.0 + Var_UB[M*N + n2]);
			HLB = log(1.0 + Var_LB[M*N + n2]);

			R = Var_UB[M*N + n2] * HLB;  //SINRj_U
			D = Var_LB[M*N + n2] * HUB;  //SINRj_L

			Coes.add(slmMinus);
			Vars.add(vars[M*N + N + n2]);//H

			Coes.add(HLB - HUB);
			Vars.add(vars[M*N + n2]);  //SINRj

			char cst_name[60];
			sprintf(cst_name, "Capacity4_D_%d", n2 + 1);
			IloConstraint cst;
			cst = IloScalProd(Coes, Vars) >= R - D;
			cst.setName(cst_name);
			model.add(cst);
			Coes.clear();
			Vars.clear();

			node2 = node2->GetNext();
		}

		//--------------------------------------------------------------------------------

		IloCplex cplex(model);
		//	cplex.setParam(cplex.EpGap,0.5);
		//	cplex.setParam(cplex.EpAGap,1.0);
		//	cplex.setParam(cplex.TiLim,10000); //限制跑的時間
		//cplex.setParam(cplex.NodeSel,2);   //啟用best-estimate search
		//	cplex.setParam(cplex.CutLo,0.0);
		//cplex.setParam(cplex.MIPSearch,1);
		//cplex.setParam(cplex.RelObjDif,0.2); //目前的解不優於最佳解的1%就跳過搜尋
		//cplex.setParam(cplex.ObjDif,10);


		cplex.solve();

		cplex.exportModel("2BS-8TP.lp");
		cout << "model.lp" << endl;


		cout << "object finish." << endl;
		cout << "Solution status = " << cplex.getStatus() << endl;
		cout << "Solution value  = " << cplex.getObjValue() << endl;

		//---------------------------------------------------------------------------------------


		char filename[250];
		FILE *fp;
		sprintf(filename, "BS2-TP8.txt");
		fp = fopen(filename, "w");
		double TT = 0;
		double p_tol = 0;

		//--------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------------------------
		fprintf(fp, "[P_D]\n");
		for (int i = 0; i< M*N; i++) {

			Vars.add(vars[i]);
			cplex.getValues(vals, Vars);
			TT = vals[0];

			if (TT != 0)
				fprintf(fp, "[%s]:\t% .7f\n", vars[i].getName(), TT);

			vals.clear();
			Vars.clear();
		}
		fprintf(fp, "\n\n");


		//----------------------------------------------------------------------------------------------
		//----------------------------------------------------------------------------------------------

		fprintf(fp, "[SINR_D]\n");
		for (int i = M * N; i<M*N + N; i++) {

			Vars.add(vars[i]);
			cplex.getValues(vals, Vars);
			TT = vals[0];

			if (TT != 0)
				fprintf(fp, "[%s]:\t% .5f\n", vars[i].getName(), TT);

			vals.clear();
			Vars.clear();
		}
		fprintf(fp, "\n\n");

		//----------------------------------------------------------------------------------------------------------
		//----------------------------------------------------------------------------------------------

		fprintf(fp, "[H_D]\n");
		for (int i = M * N + N; i<M*N + 2 * N; i++) {

			Vars.add(vars[i]);
			cplex.getValues(vals, Vars);
			TT = vals[0];

			if (TT != 0)
				fprintf(fp, "[%s]:\t% .5f\n", vars[i].getName(), TT);

			vals.clear();
			Vars.clear();
		}
		fprintf(fp, "\n\n");
		//----------------------------------------------------------------------------------------------------------

		fprintf(fp, "[ Z ]\n");

		node = this->GetNode_B();
		while (node != NULL) {
			n = node->GetValue()->GetNo();

			node2 = this->GetNode_T();
			while (node2 != NULL) {
				n2 = node2->GetValue()->GetNo();


				node4 = this->GetNode_T();
				while (node4 != NULL) {
					n4 = node4->GetValue()->GetNo();


					if (n2 != n4) {

						Vars.add(bilvars[n*N*N + n2 * N + n4]);
						cplex.getValues(vals, Vars);
						TT = vals[0];

						if (TT != 0)
							fprintf(fp, "[%s]:\t %.7f\n", bilvars[n*N*N + n2 * N + n4].getName(), TT);

						vals.clear();
						Vars.clear();

					}


					node4 = node4->GetNext();
				}//n4 

				node2 = node2->GetNext();
			}//n2 

			node = node->GetNext();
		}

		fprintf(fp, "\n\n");

		//-------------------------------------------------------------------------
		fprintf(fp, "[K ]\n");


		node2 = this->GetNode_T();
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();


			node3 = this->GetNode_B();
			while (node3 != NULL) {
				n3 = node3->GetValue()->GetNo();


				node4 = this->GetNode_T();
				while (node4 != NULL) {
					n4 = node4->GetValue()->GetNo();


					if (n2 != n4) {

						Vars.add(bilvars[M*N*N + n2 * M*N + n3 * N + n4]);
						cplex.getValues(vals, Vars);
						TT = vals[0];

						if (TT != 0)
							fprintf(fp, "[%s]:\t %.7f\n", bilvars[M*N*N + n2 * M*N + n3 * N + n4].getName(), TT);

						vals.clear();
						Vars.clear();

					}

					node4 = node4->GetNext();
				}//n4 

				node3 = node3->GetNext();
			}//n3 

			node2 = node2->GetNext();
		}//n2 

		fprintf(fp, "\n\n");

		//-----------------------------------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------								


		fprintf(fp, "Time:%2f\n", cplex.getTime());
		fprintf(fp, "epagap:%2f\n", cplex.getParam(cplex.EpAGap));
		fprintf(fp, "epgap:%2f\n", cplex.getParam(cplex.EpGap));

		fclose(fp);


		//cplex.writeSolution("3-2-100.sol");       

		//---------------------------------------------

		cplex.clear();
		vals.clear();
		Vars.clear();
		Coes.clear();
		vars.clear();
		model.end();


		//--------------------------------------------------------------------------------


	}
	catch (IloException& e) {}//try 

	return;

}//func

