#include<iostream>
using namespace std;

/*
��ά����������������е��غ�ָ�������õ�
���������Ω���ķ����У����ڲ²������Կ
��a��b������Ҫͨ�������Ķε��غ�ָ������
���ǳ��ӽ�0.065��������Ϊ�����Կ��a��b��
����ȷ�ġ�
*/


/*��ͺ�������fi*(fi-1)����i=1�ۼӵ�i=26������Ic�ķ��Ӳ���*/
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
	char temp;	//���ڶ�ȡÿ���ַ�
	int tab[27] = { 0 };		//ͳ�Ƹ���ĸ���ֵ�Ƶ������0���ǿո������1-26������ĸAa��Zz�ĸ���
	int n = 0;		//ͳ���ֳ�
	while (cin.get(temp))		//ѭ����ȡֱ��������Ϊֹ
	{
		if (temp == ' ')
			tab[0]++;		//��Ϊ�ո񣬵�0���Լ�һ
		else if (temp >= 97 && temp <= 122)
			temp -= 32;		//��ΪСд��ĸ����תΪ��д��ĸ��ͳһ����
		if (temp >= 65 && temp <= 90)
		{
			tab[temp - 64]++;		//��Ӧ��ĸƵ���Լ�һ
			n++;		//���ֳ��Լ�һ
		}
	}
	double Ic;	
	Ic = sigma(tab, 1, 26) / double(n*(n - 1));		//�ù�ʽ����غ�ָ��
	cout << "Ic(x) = " << Ic << endl;

	return 0;
}