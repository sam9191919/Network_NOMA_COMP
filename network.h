#pragma once 
#include "lib.h"
using namespace std;

 
#define CCplex_Presolve				true
#define CCplex_Alg					CPX_ALG_AUTOMATIC
#define haste true	

//class LINK;
class NODE;
class NETWORK;
 

class NODE{
	private:
		float X1, Y1,Z1,WW;	// BS coordinate (x,y)
		float X2, Y2,Z2;	// TP coordinate (x,y)
		float Pmax,Interval,Gain,d ,Flow ,Traffic;//
		int No;				// serial number
		int f , t;
	 
		//friend class LINK;
		friend class NETWORK;

	public:
		
		NODE(float x1, float y1 , float z1 ):X1(x1),Y1(y1),Z1(z1){}
		NODE(float x2 , float y2 ,float z2 ,float ww):X2(x2),Y2(y2),Z2(z2),WW(ww){}

		NODE( float gain  ):  Gain(gain){}
		int	  GetNo(){return No;}
		float GetX1(){return X1;}
		float GetY1(){return Y1;}
		float GetZ1(){return Z1;}
		float GetX2(){return X2;}
		float GetY2(){return Y2;}
		float GetZ2(){return Z2;}
        
		 
		float GetPmax(){return Pmax;}
		void  SetPmax(float p){Pmax = p;}
 
		void  SetGain(float gain){Gain = gain;}
	 	float GetGain(){return Gain;}

		float GetFlow(){return Flow;}
		void  SetFlow(float f){Flow = f;}

		
		float GetTraffic(){return Traffic;}
		void  SetTraffic(float t){Traffic = t;}

		void  SetNo(int no){No = no;}

		
};


class LINK{
	private:
		NODE *Tx, *Rx, *Node;			// refference of Tx and Rx node of this link
		float Gain, Dis;	
		int No;					// serial number
	
	public:
		
		LINK(NODE *tx, NODE *rx):Tx(tx), Rx(rx),Gain(0){}
		LINK(NODE *tx, NODE *rx, float gain, int n):Tx(tx), Rx(rx), Gain(gain){}

		
		void SetGain(float gain){Gain = gain;}
		float GetGain(){return Gain;}
 
		NODE *GetTx(){return Tx;}
		NODE *GetRx(){return Rx;}
		
		int GetNo(){return No;}
		void SetNo(int no){No = no;}

	
};
class NETWORK{
	protected:
		int N;					//Number of TP nodes
		int M;					//Number of BS nodes 
		int L;					//Number of links
	

	 
		float NoiseDensity;		//AWGN(Additive White Generate Noise) density
		float Pmax;				//maximal power
		float Obj;				//Objective value
		float Flow,Gain;
	
		LIST<NODE*>	 *Node_B,*Node_B2;
		LIST<NODE*>	 *Node_T;
		LIST<NODE*>  *node, *node2;
		LIST<NODE*>  *node3,*node4;

		//LIST<LINK*>  *Link;
		LIST<NETWORK*> *network;

		void SetDL_Link();
		void SetUL_Link();
		void SetVars();
 

	public:
		NETWORK():N(0), L(0),M(0),network(NULL),
		Node_T(NULL) ,Node_B(NULL),Node_B2(NULL),
		/*node(NULL) ,node2(NULL),node3(NULL),node4(NULL),*/
		NoiseDensity(4.83e-12), Pmax(40), Obj(0){}
	 
		float *Var_LB;
		float *Var_UB;
		float *BilVar_LB;
		float *BilVar_UB;
		int	  *BinVar_LB;
		int   *BinVar_UB;
	 
		float GetFlow(){return Flow;}
	
		void  SetObj(float obj){Obj = obj;}
		float GetObj(){return Obj;}

		void  SetPmax(float pmax){Pmax = pmax;}
		float GetPmax(){return Pmax;}

	 	float GetGain(){return Gain;}
		void  SetGain(float gain){Gain = gain;}
	 
	 
		LIST<NODE*>    *GetNode_T(){return Node_T;}
		LIST<NODE*>    *GetNode_B(){return Node_B;}
 		//LIST<LINK*>	   *GetLink(){return Link;}
	
		int GetNumberOfNode_T(){return N;}
		int GetNumberOfNode_B(){return M;}
		//int GetNumberOfLink(){return L;}
	
		//void  SetPathLoss(float pathloss){PathLoss = pathloss;}	
		//float GetPathLoss(){return PathLoss;}

		void  SetNoiseDensity(float noisedensity){NoiseDensity = noisedensity;}
		float GetNoiseDensity(){return NoiseDensity;}

		
		void AddNode_T(NODE *node);
		void AddNode_B(NODE *node);
	    void AddNode_B2(NODE *node4); 
		void AddLink(LINK *link);
		
		void SetD_Gain();  
		void AddUBLB(); 		// 宣告變數上下限
		void AddVars(); 		//變數位置
		void func();    		//建立MODEL


};

