#pragma once 
#include <iostream>
#include <ilcplex/ilocplex.h>

template <class TYPE>
class LIST{
	private:
		LIST *Prev;			//前一個元素
		LIST *Next;			//後一個元素
		TYPE Value;			//值
	public:
		LIST(){				//建構子：	
			Prev = NULL;
			Next = NULL;
		}

		LIST(TYPE value){	//建構子：值
			Prev = NULL;
			Next = NULL;
			Value = value;
		}

		LIST(LIST<TYPE> *prev, TYPE value, LIST<TYPE> *next){	//建構子：前元素,值,後元素
			Prev = prev;
			Next = next;
			Value = value;
		}

		~LIST(){			//解構子
			if(Prev != NULL)
				Prev->SetNext(Next);
			if(Next != NULL)
				Next->SetPrev(Prev);
		}

		TYPE GetValue(void){return Value;}	//取得：值
			
		TYPE *GetValueRef(void){return &Value;}	//取得：值
		
		LIST *GetPrev(void){return Prev;}	//取得：前一個元素
		
		LIST *GetNext(void){return Next;}	//取得：後一個元素

		void SetValue(TYPE value){Value = value;}	//設定：值
		
		void SetPrev(LIST *prev){Prev = prev;}	//設定：前一個元素
		
		void SetNext(LIST *next){Next = next;}	//設定：後一個元素
		
		LIST *Head(void){		//取得：首項元素
			LIST *current=this;
			while(current->GetPrev() != NULL)
				current = current->GetPrev();
			return current;
		}

		LIST *Tail(void){		//取得：末項元素
			LIST *current=this;
			while(current->GetNext() != NULL)
				current = current->GetNext();		
			return current;
		}

		LIST *Find(TYPE value){	//尋找：值所在的最前端元素
			LIST<TYPE> *current=Head();
			do{
				if(current->GetValue() == value)
					return current;
				else;
				current = current->GetNext();
			}while(current != NULL);
			return NULL;
		}

		int GetLength(void){		//取得：LIST長度(元素數量)
			int i=0;
			LIST<TYPE> *current;
			current = Head();
			while(current!=NULL){
				i++;
				current = current->GetNext();
			}
			return i;
		}

		void Remove(){			//移除：該項元素
			if(Prev != NULL)
				Prev->SetNext(Next);
			if(Next != NULL)
				Next->SetPrev(Prev);
			//cout<<"asdfasdfsdaf"<<endl;
			delete this;
		}

		LIST& operator+(LIST& RHS){	//運算子重載： + 連結兩LIST
			RHS.Head()->SetPrev(this->Tail());
			Tail()->SetNext(RHS.Head());
		}
};


