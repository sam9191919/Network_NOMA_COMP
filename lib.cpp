#pragma once 
#include "lib.h"


/*
LIST& LIST<TYPE> :: operator+(LIST& RHS){	//運算子重載： + 連結兩LIST
	RHS.Head()->SetPrev(this->Tail());
	Tail()->SetNext(RHS.Head());
}
*/