#include <NTL/ZZ.h>
#include <iostream>
#include <time.h>
#include<string>
#include<cstring>
#include<fstream>
NTL_CLIENT
using namespace std;
using namespace NTL;

class ta;
class user;

/*���Բ���*/
bool is_prime(ZZ p, int n = 50)
{
	ZZ a;
	for (int i = 0; i < n; i++)		//Ĭ�ϲ���n=50�Σ����ʹ������
	{
		a = RandomBnd(p - 1);
		a++;
		if (MillerWitness(p, a))		//Miller-Rabin����
			return false;
	}
	return true;
}

/*	ȡԭ��	*/
void get_pri_root(ZZ &a, ZZ p,ZZ p0)
{
	/*p����p0���ɵ�������p=2*p0+1*/
	ZZ t;
	while (1)
	{
		t = RandomBnd(p);
		if (t <= 2)
			continue;
		if (PowerMod(t, (p - 1) / 2, p) == 1)		//������Ϊ2
			continue;
		if (PowerMod(t, (p - 1) / p0, p) == 1)		//������Ϊp0
			continue;
		break;
	}
	a = t;
}

/*��������������*/
ZZ connect(ZZ a, ZZ b)
{
	ZZ t, c, n;
	c = 10;
	n = 1;
	t = c;
	/*�ҵ���b�����С��10����������*/
	while (t <= b)
	{
		t = t*c;
		n++;
	}
	return a*t + b;
}

/*�ַ���ת������*/
ZZ str_to_zz(string s)
{
	int len = s.length();
	ZZ n;
	n = 0;
	/*���ս���ת������*/
	for (int i = len - 1; i >= 0; i--)
	{
		/*�����ַ�����ÿһλ�ַ�*/
		int x = len - 1 - i;
		ZZ t;
		t = 1;
		/*ȷ����ǰλ�ϵ�Ȩ*/
		for (int j = 1; j <= x; j++)
			t = t * 256;
		n = n + t*s[i];
	}
	return n;
}

/*RSA��Կ��*/
struct R_key
{
	struct
	{
		ZZ n, b;
	}r_pub;		//��Կ
	
	struct
	{
		ZZ p, q, a;
	}r_pri;		//˽Կ

	ZZ sig;		//˽Կƴ�ӽ��
	ZZ ver;		//��Կƴ�ӽ��
};
/*ElGamal��Կ��*/
struct E_key
{
	struct
	{
		ZZ p, alpha, beta;
	}e_pub;	//��Կ
	
	struct
	{
		ZZ a;
	}e_pri;		//˽Կ

	ZZ sig;		//˽Կƴ�ӽ��
	ZZ ver;		//��Կƴ�ӽ��
};

/*RSA*/
class RSA
{
public:
	R_key rk;
	

	/*����RSA��Կ��Ĭ��512λ*/
	void GenKey(int bits=512)
	{
		ZZ n, b;		//��Կ
		ZZ p, q, a;	//˽Կ
		ZZ phi;			//ŷ������

		SetSeed(conv<ZZ>(static_cast<long>(time(nullptr))));

		/*	�����������p, q	*/
		GenPrime(p, bits);
		GenPrime(q, bits);

		/*	�������n��phi(n)	*/
		n = p*q;
		phi = (p - 1)*(q - 1);

		/*	���������b������b����Ԫa	*/
		do
			RandomBnd(b, phi);
		while (!(b > 1 && GCD(b, phi) == 1));
		InvMod(a, b, phi);		//b=a^(-1) mod phi

		/*����Կ��ֵ����Ӧ����Կ���������*/
		rk.r_pub.n = n;
		rk.r_pub.b = b;
		rk.ver = connect(n, b);

		rk.r_pri.p = p;
		rk.r_pri.q = q;
		rk.r_pri.a = a;
		rk.sig = connect(connect(p, q), a);
	}
	/*RSAǩ������*/
	void RSA_sig(ZZ x, ZZ &y)
	{
		PowerMod(y, x, rk.r_pri.a, rk.r_pub.n);		//sig(x)=x^a mod n
	}
	/*RSAǩ����֤*/
	bool RSA_ver(ZZ x, ZZ y) 
	{
		ZZ xx;
		PowerMod(xx, y, rk.r_pub.b, rk.r_pub.n);
		/*�ж�x=y^b mod n �Ƿ����*/
		if (xx == x)
			return true;
		else
			return false;
	}
};
/*ElGamal*/
class ElGamal
{
public:
	E_key ek;
	ZZ q, k;	//qΪ��Կp��ǰ��kΪѡȡ�������
	

	/*����ElGamal��Կ��Ĭ�ϳ���1024λ*/
	void GenKey(int bits=1024)
	{
		ZZ p, alpha, beta;	//��Կ
		ZZ a;						//˽Կ
		
		ZZ t1, t2;

		SetSeed(conv<ZZ>(static_cast<long>(time(nullptr))));

		/*	��������p	*/
		do
		{
			GenGermainPrime(q, bits - 1);
			p = 2 * q + 1;
		} while (!is_prime(p));

		/*	���ɱ�ԭԪalpha	*/
		get_pri_root(alpha, p, q);

		/*	���������k, a */
		do
			RandomBnd(k, p - 1);
		while (!(k != 0 && GCD(k, p - 1) == 1));

		do
			RandomBnd(a, p - 1);
		while (a == 0);

		beta = PowerMod(alpha, a, p);		//beta=alpha^a mod p

		/*����Կ��ֵ����Կ���������*/
		ek.e_pub.p = p;
		ek.e_pub.alpha = alpha;
		ek.e_pub.beta = beta;
		ek.ver = connect(connect(p, alpha), beta);

		ek.e_pri.a = a;
		ek.sig = a;
	}
	/*ElGamalǩ������*/
	void ElGamal_sig(ZZ x, ZZ &gamma, ZZ &delta)
	{
		ZZ t1, t2, t3;
		t3=PowerMod(ek.e_pub.alpha, k, ek.e_pub.p);		//gamma=alpha^k (mod p)
		gamma = t3;
		t1=InvMod( k, ek.e_pub.p - 1);		//t1=k^-1 mod p-1
		t2 = x - ek.e_pri.a*t3;			//t2=x-a*gamma
		t2 = t2 % (ek.e_pub.p - 1);			//t2=x-a*gamma mod p-1
		MulMod(delta, t2, t1, ek.e_pub.p - 1);			//delta=(x-a*gamma)*k^-1 (mod p-1)
	}
	/*ElGamalǩ����֤*/
	bool ElGamal_ver(ZZ x, ZZ gamma, ZZ delta)
	{
		ZZ left, right,j,k;
		PowerMod(right, ek.e_pub.alpha, x,ek.e_pub.p);			//�ұ�=alpha^x mod p
		PowerMod(j, ek.e_pub.beta, gamma, ek.e_pub.p);
		PowerMod(k, gamma, delta, ek.e_pub.p);
		MulMod(left, j, k, ek.e_pub.p);			//���=beta^gamma*gamma^delta mod p
		if (compare(left,right)==0)
			return true;
		else
			return false;
	}
};

/*�û�*/
class user
{
	friend ta;

public:
	string name;
	ZZ id;
	int method;	 //0δǩ����1RSA��2ElGamal
	ZZ r;  //��ս
	ZZ aT, bT;		//����ָ�����乫��ֵ

	/*TA��RSAǩ�������䲼��֤��*/
	struct RCert
	{
		ZZ id_ver_s;
		ZZ id_ver;
		ZZ s;
		ZZ bT;
	} rcert;
	
	/*TA��ElGamalǩ�������䲼��֤��*/
	struct ECert
	{
		ZZ gamma, delta;
		ZZ id_ver_s;
		ZZ id_ver;
		ZZ s;
		ZZ bT;
	} ecert;

	RSA rsa;
	ElGamal elg;

	user(string s)
	{
		method = 0;
		name = s;
	}
};
/*����Ȩ������*/
class ta
{
public:
	int method;
	RSA rsa;
	ElGamal elg;

	ta()
	{
		method = 0;
	}

	/*�����û���ID*/
	void GenID(user &p)
	{
		p.id = str_to_zz(p.name);
	}
	/*�����û�֤��*/
	void GenCert(user &p,int bits)
	{
		ZZ s, c,t;
		if (method == 1)
		{
			/*����TA ����Կ*/
			rsa.GenKey(bits);

			/*����֤��*/
			t = connect(p.id, p.rsa.rk.ver);
			rsa.RSA_sig(t, s);		//##################################
			c = connect(t, s);

			/*�䲼֤��*/
			p.rcert.id_ver_s = c;
			p.rcert.id_ver = t;
			p.rcert.s = s;
			p.rcert.bT = p.bT;
		}
		else
		{
			/*����TA����Կ*/
			elg.GenKey(bits);

			/*����֤��*/
			t = connect(p.id, p.elg.ek.ver);
			ZZ gamma, delta;
			elg.ElGamal_sig(t, gamma, delta);
			s = connect(gamma, delta);
			c = connect(t, s);

			/*�䲼֤��*/
			p.ecert.gamma = gamma;
			p.ecert.delta = delta;
			p.ecert.id_ver_s = c;
			p.ecert.id_ver = t;
			p.ecert.s = s;
			p.ecert.bT = p.bT;
		}
	}
	/*����֤���ļ�*/
	void GenCertFile(user p)
	{
		ofstream fout;
		fout.open("Certification.txt", ios::out);
		if (!fout.is_open())
		{
			cout << "���ļ�ʧ��" << endl;
			return;
		}
		if (method == 1)
			fout << p.rcert.id_ver_s;
		else
			fout << p.ecert.id_ver_s;
		fout.close();
	}
	/*֤����֤����*/
	bool Cert_ver(user p)
	{
		if (method == 1)
		{
			return rsa.RSA_ver(p.rcert.id_ver, p.rcert.s);
		}
		else
		{
			return elg.ElGamal_ver(p.ecert.id_ver, p.ecert.gamma, p.ecert.delta);
		}
	}
};

/*RSAǩ������*/
void RSA_test()
{
	user Alice("Alice");
	Alice.method = 1;

	ZZ x, y;

	cout << "��������Ϣ��" << endl;
	cin >> x;

	/*Alice ������Կ*/
	Alice.rsa.GenKey();
	cout << endl<<"������Alice����Կ���밴Enter����������" << endl;
	getchar();
	getchar();

	/*Alice ǩ��*/
	Alice.rsa.RSA_sig(x, y);
	cout << "Alice�����ǩ��"<<endl;
	cout <<endl<< y << endl;
	cout << endl<<"�밴Enter����������" << endl;
	getchar();

}

/*ElGamalǩ������*/
void ElGamal_test()
{
	user Alice("Alice");
	Alice.method = 2;

	ZZ x;
	struct
	{
		ZZ gamma;
		ZZ delta;
	} y;

	cout << "��������Ϣ��" << endl;
	cin >> x;

	/*Alice ������Կ*/
	Alice.elg.GenKey();
	cout << endl<<"������Alice����Կ���밴Enter����������"<<endl;
	getchar();
	getchar();

	/*Alice ǩ��*/
	Alice.elg.ElGamal_sig(x, y.gamma, y.delta);
	cout <<endl<< "Alice�����ǩ��" << endl;
	cout << endl<< y.gamma<<'  '<<y.delta << endl;
	cout << endl << "�밴Enter����������" << endl;
	getchar();

}

/*֤��İ䲼����֤����*/
void Cert_test()
{
	const int bits = 64;

	user Alice("Alice");
	user Bob("Bob");
	ta TA;

	int choice;

	cout << "��ѡ��TA��ǩ����ʽ��1.RSA  2.ElGamal����";
	cin >> choice;
	TA.method = choice;
	cout << endl;

	cout << "����Alice����ѡ����Ҫ��ǩ����ʽ��1.RSA  2.ElGamal����";
	cin >> choice;
	Alice.method = choice;
	cout << endl;

	/*Alice ������Կ*/
	if (Alice.method == 1)
		Alice.rsa.GenKey(bits);
	else
		Alice.elg.GenKey(bits);
	cout << "������Alice����Կ���밴Enter����������";
	getchar();
	getchar();

	/*TA ���� Alice ID*/
	TA.GenID(Alice);
	cout << "������Alice��ID���밴Enter����������";
	getchar();

	/*TA �� Alice �䲼֤�飬 �������ļ�*/
	TA.GenCert(Alice,4*bits);		//##################################
	cout << "������Alice��֤�飬�밴Enter����������";
	getchar();
	TA.GenCertFile(Alice);
	cout << "������Alice��֤���ļ����밴Enter����������";
	getchar();

	/*Bob ���� TA ����֤ Alice ��֤��*/
	bool TorF=TA.Cert_ver(Alice);
	cout << "Bob����֤Alice��֤�飬���Ϊ��";
	if (TorF)
		cout << "True.    ";
	else
		cout << "False.    ";
	cout<<"�밴Enter����������";
	getchar();
	cout << endl;
	return;
}

/*������֤����*/
void Interact_Cert_test()
{
	const int bits = 64;

	user Alice("Alice");
	user Bob("Bob");
	ta TA;

	int choice;

	cout << "��ѡ��TA��ǩ��������1. RSA  2. ElGamal����";
	cin >> choice;
	TA.method = choice;

	cout << "��ѡ��Alice��ǩ��������1. RSA  2. ElGamal����";
	cin >> choice;
	Alice.method = choice;

	cout << "��ѡ��Bob��ǩ��������1. RSA  2. ElGamal����";
	cin >> choice;
	Bob.method = choice;

	cout << endl;

	if (Alice.method == 1)
		Alice.rsa.GenKey(bits);
	else
		Alice.elg.GenKey(bits);
	cout << "������Alice����Կ���밴Enter����������";
	getchar();
	getchar();

	/*Alice ����֤�� Cert(Alice)*/
	TA.GenID(Alice);
	cout << "������Alice��ID���밴Enter����������";
	getchar();
	TA.GenCert(Alice,4*bits);
	cout << "������Alice��֤�飬�밴Enter����������";
	getchar();

	/*Bob ��֤Alice �Ĺ�Կ��������Alice �Ĺ�Կ��֤*/
	bool TorFa = TA.Cert_ver(Alice);
	cout << "Bob����֤Alice�Ĺ�Կ�����Ϊ��";
	if (TorFa)
		cout << "True.    ";
	else
		cout << "False.    ";
	cout << "�밴Enter����������";
	getchar();

	if (Bob.method == 1)
		Bob.rsa.GenKey(bits);
	else
		Bob.elg.GenKey(bits);
	cout << "������Bob����Կ���밴Enter����������";
	getchar();

	/*Bob ����֤�� Cert(Bob)*/
	TA.GenID(Bob);
	cout << "������Bob��ID���밴Enter����������";
	getchar();
	TA.GenCert(Bob,4*bits);
	cout << "������Bob��֤�飬�밴Enter����������";
	getchar();

	

	/*Alice ��֤Bob �Ĺ�Կ��������Bob �Ĺ�Կ��֤*/
	bool TorFb = TA.Cert_ver(Bob);
	cout << "Alice����֤Bob�Ĺ�Կ�����Ϊ��";
	if (TorFb)
		cout << "True.    ";
	else
		cout << "False.    ";
	cout << "�밴Enter����������";
	getchar();

	SetSeed(conv<ZZ>(static_cast<long>(time(nullptr))));
	/*Bob ���������ս r1*/
	ZZ r1;
	RandomBits(r1, bits/8);
	Bob.r = r1;
	cout << "������Bob�������սr1���밴Enter����������";
	getchar();

	/*Alice ���������ս r2*/
	ZZ r2;
	RandomBits(r2, bits/8);
	Alice.r = r2;
	cout << "������Alice�������սr2���밴Enter����������";
	getchar();

	/*Alice ���� y1*/
	struct
	{
		ZZ rsa_y;
		ZZ elg_gamma;
		ZZ elg_delta;
	} y1;
	if (Alice.method == 1)
		Alice.rsa.RSA_sig(connect(connect(Bob.id, Bob.r), Alice.r), y1.rsa_y);
	else
		Alice.elg.ElGamal_sig(connect(connect(Bob.id, Bob.r), Alice.r),y1.elg_gamma,y1.elg_delta);
	cout << "Alice�Ѽ���y1���밴Enter����������";
	getchar();


	bool AorR;
	if (Alice.method == 1)
		AorR=Alice.rsa.RSA_ver(connect(connect(Bob.id, Bob.r), Alice.r), y1.rsa_y);
	else
		AorR = Alice.elg.ElGamal_ver(connect(connect(Bob.id, Bob.r), Alice.r), y1.elg_gamma, y1.elg_delta);
	if (AorR)
		cout << "Bob ���ܣ��밴Enter����������";
	else
		cout << "Bob �ܾ����밴Enter����������";
	getchar();

	/*Bob ���� y2*/
	struct
	{
		ZZ rsa_y;
		ZZ elg_gamma;
		ZZ elg_delta;
	} y2;
	if (Bob.method == 1)
		Bob.rsa.RSA_sig(connect(Alice.id, Alice.r), y2.rsa_y);
	else
		Bob.elg.ElGamal_sig(connect(Alice.id, Alice.r), y2.elg_gamma, y2.elg_delta);
	cout << "Bob�Ѽ���y2���밴Enter����������";
	getchar();

	if (Bob.method == 1)
		AorR=Bob.rsa.RSA_ver(connect(Alice.id, Alice.r), y2.rsa_y);
	else
		AorR=Bob.elg.ElGamal_ver(connect(Alice.id, Alice.r), y2.elg_gamma, y2.elg_delta);
	if (AorR)
		cout << "Alice ���ܣ��밴Enter����������";
	else
		cout << "Alice �ܾ����밴Enter����������";
	getchar();
	cout << endl;
}

/*MTI/A0��ԿЭ�̷�������*/
void MTI_A0_test()
{

	const int bits = 1024;
	ZZ p, a, n;

	user U("U"), V("V");
	ta TA;
	int choice;
	
	cout << "��ѡ��TA��ǩ��������1. RSA  2. ElGamal����";
	cin >> choice;
	TA.method = choice;

	cout << "��ѡ��U��ǩ��������1. RSA  2. ElGamal����";
	cin >> choice;
	U.method = choice;

	cout << "��ѡ��V��ǩ��������1. RSA  2. ElGamal����";
	cin >> choice;
	V.method = choice;

	SetSeed(conv<ZZ>(static_cast<long>(time(nullptr))));
	
	/*	��������p*/
	ZZ q;
	do
	{
		GenGermainPrime(q, bits - 1);
		p = 2 * q + 1;
	} while (!is_prime(p));

	/*	���ɱ�ԭԪa	*/
	get_pri_root(a, p, q);
	n = p - 1;

	cout << endl;
	if (U.method == 1)
		U.rsa.GenKey(bits/8);
	else
		U.elg.GenKey(bits/8);
	cout << "������U����Կ���밴Enter����������";
	getchar();
	getchar();

	if (V.method == 1)
		V.rsa.GenKey(bits/8);
	else
		V.elg.GenKey(bits/8);
	cout << "������V����Կ���밴Enter����������";
	getchar();

	SetSeed(conv<ZZ>(static_cast<long>(time(nullptr))));
	/*Uѡȡ����ָ���������㹫��ֵ������֤��*/
	RandomBnd(U.aT, n);
	PowerMod(U.bT, a, U.aT, p);
	cout << "������U������ָ�����乫��ֵ���밴Enter����������";
	getchar();

	TA.GenID(U);
	TA.GenCert(U,bits/2);
	cout << "������U��֤�飬�밴Enter����������";
	getchar();

	/*Vѡȡ����ָ���������㹫��ֵ������֤��*/
	RandomBnd(V.aT, n);
	PowerMod(V.bT, a, V.aT, p);
	cout << "������V������ָ�����乫��ֵ���밴Enter����������";
	getchar();

	TA.GenID(V);
	TA.GenCert(V,bits/2);
	cout << "������V��֤�飬�밴Enter����������";
	getchar();

	/*Uѡȡ���ֵrU��������sU*/
	ZZ rU, sU;
	RandomBnd(rU, n);
	PowerMod(sU, a, rU, p); 
	cout << "U��ѡȡ���ֵrU��������˶�Ӧ��sU���밴Enter����������";
	getchar();

	/*Vѡȡ���ֵrV��������sV*/
	ZZ rV, sV;
	RandomBnd(rV, n);
	PowerMod(sV, a, rV, p);
	cout << "V��ѡȡ���ֵrV��������˶�Ӧ��sV���밴Enter����������";
	getchar();


	ZZ ku1, ku2;
	/*U����Ự��Կ*/
	ZZ KU;
	PowerMod(ku1, sV, U.aT, p);
	if (TA.method == 1)
		PowerMod(ku2, V.rcert.bT, rU, p);
	else
		PowerMod(ku2, V.ecert.bT, rU, p);
	MulMod(KU, ku1, ku2, p);
	cout << "U�Ѽ�����Ự��ԿKU���밴Enter����������";
	getchar();

	ZZ kv1, kv2;
	/*V����Ự��Կ*/
	ZZ KV;
	PowerMod(kv1, sU, V.aT, p);
	if (TA.method == 1)
		PowerMod(kv2, U.rcert.bT, rV, p);
	else
		PowerMod(kv2, U.ecert.bT, rV, p);
	MulMod(KV, kv1, kv2, p);
	cout << "V�Ѽ�����Ự��ԿKV���밴Enter����������";
	getchar();

	cout << endl;
	/*�жϻỰ��Կ�Ƿ���ͬ*/
	cout << "�ѶԱ�U��V��������ĻỰ��Կ�����Ϊ��";
	if (compare(KU, KV) == 0)
	{
		cout << "��ͬ.    " << endl;
		cout << endl;
		cout << KU << endl;
		cout << "�밴Enter����������";
	}
	else
	{
		cout << "��ͬ.    ";
		cout << "�밴Enter����������";
	}
	
	
	getchar();
	cout << endl;
}

int main()
{
	system("cls");
	
	while (1)
	{
		cout << "��ѡ��" << endl;
		cout << "1. RSAǩ������" << endl;
		cout << "2. ElGamalǩ������" << endl;
		cout << "3. ֤��İ䲼����֤" << endl;
		cout << "4. ������֤" << endl;
		cout << "5. MTI/A0��ԿЭ�̷���" << endl;
		cout << endl;
		cout << "�������Ӧ����ţ���Q/q�˳�" << endl << endl;
		char choice;
		cin >> choice;

		system("cls");

		/*ѡ���˳�*/
		if (choice == 'Q' || choice == 'q')
		{
			cout << "���˳�" << endl << endl;
			break;
		}

		/*ѡ��1*/
		else if (choice == '1')
		{
			cout << "1. RSAǩ������" << endl << endl;
			RSA_test();
		}

		/*ѡ��2*/
		else if (choice == '2')
		{
			cout << "2. ElGamalǩ������" << endl << endl;
			ElGamal_test();
		}

		/*ѡ��3*/
		else if (choice == '3')
		{
			cout << "3. ֤��İ䲼����֤" << endl << endl;
			Cert_test();
		}

		/*ѡ��4*/
		else if (choice == '4')
		{
			cout << "4. ������֤" << endl << endl;
			Interact_Cert_test();
		}

		/*ѡ��5*/
		else if (choice == '5')
		{
			cout << "5. MTI/A0��ԿЭ�̷���" << endl << endl;
			MTI_A0_test();
		}

		/*ѡ����Ч*/
		else
			continue;

		cout << "�밴Enter���������˵�����";
		getchar();
		system("cls");

	}

	return 0;
}