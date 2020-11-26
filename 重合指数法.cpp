#include<iostream>
using namespace std;

/*
将维吉尼亚密码分析法中的重合指数法运用到
仿射密码的惟密文分析中，对于猜测出的密钥
（a，b），需要通过计算文段的重合指数，如
果非常接近0.065，即可认为这个密钥（a，b）
是正确的。
*/


/*求和函数，求fi*(fi-1)，从i=1累加到i=26，即求Ic的分子部分*/
int sigma(int a[],int buttom,int top)
{
	int sum = 0;
	for (int i = buttom; i <= top; i++)
	{
		sum += a[i] * (a[i] - 1);
	}
	return sum;
}

int main()
{
	char temp;	//用于读取每个字符
	int tab[27] = { 0 };		//统计各字母出现的频数，第0项是空格个数，1-26项是字母Aa到Zz的个数
	int n = 0;		//统计字长
	while (cin.get(temp))		//循环读取直到读不到为止
	{
		if (temp == ' ')
			tab[0]++;		//若为空格，第0项自加一
		else if (temp >= 97 && temp <= 122)
			temp -= 32;		//若为小写字母，暂转为大写字母，统一计算
		if (temp >= 65 && temp <= 90)
		{
			tab[temp - 64]++;		//对应字母频数自加一
			n++;		//总字长自加一
		}
	}
	double Ic;	
	Ic = sigma(tab, 1, 26) / double(n*(n - 1));		//用公式求出重合指数
	cout << "Ic(x) = " << Ic << endl;

	return 0;
}