#pragma once 
#include "network.h"
#include <vector>
#include <stack>

using namespace std; 

#define DDisplay			true
#define NNoiseDensity		(float)4.113294e-12	//for channel gain


void generate_network(NETWORK *network);	
//void initial_prog();
//void print_var(NETWORK *network);


int main (void){	


	//########## Setup Network ##########

  	NETWORK network;
	generate_network(&network);	//generate network with node.txt,W.txt,session.txt
    
    network.AddUBLB( ); 
    network.func( );

	system("pause");

}

void generate_network(NETWORK *network){


	//if(DDisplay)
	//cout << "Data loading...." << endl 
	//	 << "Set and Clean Links...." << endl
	//	 << "Set Vars...." << endl
	//	 << "Set Interference Table...." << endl;
	  cout<<"Femtell-EE "<<endl;

	FILE *fp;
    

// Load TP node

	fp = fopen("10TPS.txt","r");
		float x1, y1, z1;

	 while(fscanf(fp,"%f %f %f ", &x1, &y1, &z1)==3){
		network->AddNode_T(new NODE (x1, y1, z1 ) );
		cout<<"x1 = " << x1 <<endl;
		cout<<"y1 = " << y1 <<endl;
		cout<<"z1 = " << z1 <<endl;
	//cout <<"t="<<t<<endl;
	 //   cout<<"flow = " << w1 <<endl;
	//	cout<<"power = " << p1 <<endl;
	} 
  	fclose(fp);


// Load BS node

	fp = fopen("BST1.txt","r");
	float x2, y2,z2,ww;


	 while(fscanf(fp,"%f %f %f %f", &x2, &y2,&z2,&ww)==4 ){
		network->AddNode_B(new NODE (x2, y2,z2,ww) );
        cout<<"x2 = " << x2 <<endl;
		cout<<"y2 = " << y2 <<endl;
		cout<<"z2 = " << z2 <<endl;
		/*node = network->GetNode_B()->Tail();
		cout<<node->GetValue()->GetX1()<<"\t"<<node->GetValue()->GetY1()<<endl;*/
		
	} 
  
	fclose(fp); 

 

return;
}


//	 
//	 
//
//// Display Network
//	//if(DDisplay){
//	//	cout<<endl<<"[Network Configuration]"<<endl;
//	//	for(LIST<NODE*> *N=network->GetNode() ; N!=NULL ; N=N->GetNext()){
//	//		LIST<LINK*> *OIL;
//	//		printf("NODE %d : (%3.2f , %3.2f)-----OL:"
//	//			,N->GetValue()->GetNo()+SHIFT
//	//			,N->GetValue()->GetX()
//	//			,N->GetValue()->GetY()
//	//			);
//	//		OIL = N->GetValue()->GetOL();
//	//		while(OIL!=NULL){
//	//			cout << OIL->GetValue()->GetNo() << ", ";
//	//			OIL = OIL->GetNext();
//	//		}
//	//		cout <<"(end)"<< endl << "                        -----IL:";
//	//		OIL = N->GetValue()->GetIL();
//	//		while(OIL!=NULL){
//	//			cout << OIL->GetValue()->GetNo() << ", ";
//	//			OIL = OIL->GetNext();
//	//		}
//	//		cout <<"(end)"<< endl;
//
//	//	}
//	
//		
//	/*	for(LIST<LINK*> *L=network->GetLink() ; L!=NULL ; L=L->GetNext()){
//			fprintf(fp,"LINK%3d: (%3.2f,%3.2f)---(%3.2f,%3.2f) "
//				,L->GetValue()->GetNo()
//				,L->GetValue()->GetTx()->GetX1()
//				,L->GetValue()->GetTx()->GetY1()
//				,L->GetValue()->GetRx()->GetX2()
//				,L->GetValue()->GetRx()->GetY2()
//				);
//			
//			fprintf(fp,", %.2e ",L->GetValue()->GetGain());
//			fprintf(fp,"\n");
//		}*/
//		//fprintf(fp,"\nNetwork: N = %d , L = %d "
//		//	,network->GetNumberOfNode()
//		//	,network->GetNumberOfLink()
//		//		
//		//	);
//		//fprintf(fp,"Vars: %d   , Constrains: %d\n"
//		//	,network->GetNumberOfSession()*network->GetNumberOfLink()+network->GetNumberOfLink()
//		//	,network->GetNumberOfLink()*2+network->GetNumberOfNode()*network->GetNumberOfSession()
//		//	);
		//system("pause");
	





