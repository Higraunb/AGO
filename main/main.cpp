#include "Point.h"
#include "iostream"
int main()
{
	TPoint<int, 3> a{ 1 ,2 , 3.6};
	operator << (std::cout, a);
}