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

/*素性测试*/
bool is_prime(ZZ p, int n = 50)
{
	ZZ a;
	for (int i = 0; i < n; i++)		//默认测试n=50次，降低错误概率
	{
		a = RandomBnd(p - 1);
		a++;
		if (MillerWitness(p, a))		//Miller-Rabin测试
			return false;
	}
	return true;
}

/*	取原根	*/
void get_pri_root(ZZ &a, ZZ p,ZZ p0)
{
	/*p是由p0生成的素数，p=2*p0+1*/
	ZZ t;
	while (1)
	{
		t = RandomBnd(p);
		if (t <= 2)
			continue;
		if (PowerMod(t, (p - 1) / 2, p) == 1)		//素因子为2
			continue;
		if (PowerMod(t, (p - 1) / p0, p) == 1)		//素因子为p0
			continue;
		break;
	}
	a = t;
}

/*连接两个大整数*/
ZZ connect(ZZ a, ZZ b)
{
	ZZ t, c, n;
	c = 10;
	n = 1;
	t = c;
	/*找到比b大的最小的10的整数次幂*/
	while (t <= b)
	{
		t = t*c;
		n++;
	}
	return a*t + b;
}

/*字符串转大整数*/
ZZ str_to_zz(string s)
{
	int len = s.length();
	ZZ n;
	n = 0;
	/*按照进制转换方法*/
	for (int i = len - 1; i >= 0; i--)
	{
		/*遍历字符串的每一位字符*/
		int x = len - 1 - i;
		ZZ t;
		t = 1;
		/*确定当前位上的权*/
		for (int j = 1; j <= x; j++)
			t = t * 256;
		n = n + t*s[i];
	}
	return n;
}

/*RSA密钥组*/
struct R_key
{
	struct
	{
		ZZ n, b;
	}r_pub;		//公钥
	
	struct
	{
		ZZ p, q, a;
	}r_pri;		//私钥

	ZZ sig;		//私钥拼接结果
	ZZ ver;		//公钥拼接结果
};
/*ElGamal密钥组*/
struct E_key
{
	struct
	{
		ZZ p, alpha, beta;
	}e_pub;	//公钥
	
	struct
	{
		ZZ a;
	}e_pri;		//私钥

	ZZ sig;		//私钥拼接结果
	ZZ ver;		//公钥拼接结果
};

/*RSA*/
class RSA
{
public:
	R_key rk;
	

	/*生成RSA密钥，默认512位*/
	void GenKey(int bits=512)
	{
		ZZ n, b;		//公钥
		ZZ p, q, a;	//私钥
		ZZ phi;			//欧拉函数

		SetSeed(conv<ZZ>(static_cast<long>(time(nullptr))));

		/*	生成随机素数p, q	*/
		GenPrime(p, bits);
		GenPrime(q, bits);

		/*	计算参数n，phi(n)	*/
		n = p*q;
		phi = (p - 1)*(q - 1);

		/*	生成随机数b，计算b的逆元a	*/
		do
			RandomBnd(b, phi);
		while (!(b > 1 && GCD(b, phi) == 1));
		InvMod(a, b, phi);		//b=a^(-1) mod phi

		/*将密钥赋值到对应的密钥组的数据中*/
		rk.r_pub.n = n;
		rk.r_pub.b = b;
		rk.ver = connect(n, b);

		rk.r_pri.p = p;
		rk.r_pri.q = q;
		rk.r_pri.a = a;
		rk.sig = connect(connect(p, q), a);
	}
	/*RSA签名方案*/
	void RSA_sig(ZZ x, ZZ &y)
	{
		PowerMod(y, x, rk.r_pri.a, rk.r_pub.n);		//sig(x)=x^a mod n
	}
	/*RSA签名验证*/
	bool RSA_ver(ZZ x, ZZ y) 
	{
		ZZ xx;
		PowerMod(xx, y, rk.r_pub.b, rk.r_pub.n);
		/*判断x=y^b mod n 是否成立*/
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
	ZZ q, k;	//q为密钥p的前身，k为选取的随机数
	

	/*生成ElGamal密钥，默认长度1024位*/
	void GenKey(int bits=1024)
	{
		ZZ p, alpha, beta;	//公钥
		ZZ a;						//私钥
		
		ZZ t1, t2;

		SetSeed(conv<ZZ>(static_cast<long>(time(nullptr))));

		/*	生成素数p	*/
		do
		{
			GenGermainPrime(q, bits - 1);
			p = 2 * q + 1;
		} while (!is_prime(p));

		/*	生成本原元alpha	*/
		get_pri_root(alpha, p, q);

		/*	生成随机数k, a */
		do
			RandomBnd(k, p - 1);
		while (!(k != 0 && GCD(k, p - 1) == 1));

		do
			RandomBnd(a, p - 1);
		while (a == 0);

		beta = PowerMod(alpha, a, p);		//beta=alpha^a mod p

		/*将密钥赋值到密钥组的数据中*/
		ek.e_pub.p = p;
		ek.e_pub.alpha = alpha;
		ek.e_pub.beta = beta;
		ek.ver = connect(connect(p, alpha), beta);

		ek.e_pri.a = a;
		ek.sig = a;
	}
	/*ElGamal签名方案*/
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
	/*ElGamal签名验证*/
	bool ElGamal_ver(ZZ x, ZZ gamma, ZZ delta)
	{
		ZZ left, right,j,k;
		PowerMod(right, ek.e_pub.alpha, x,ek.e_pub.p);			//右边=alpha^x mod p
		PowerMod(j, ek.e_pub.beta, gamma, ek.e_pub.p);
		PowerMod(k, gamma, delta, ek.e_pub.p);
		MulMod(left, j, k, ek.e_pub.p);			//左边=beta^gamma*gamma^delta mod p
		if (compare(left,right)==0)
			return true;
		else
			return false;
	}
};

/*用户*/
class user
{
	friend ta;

public:
	string name;
	ZZ id;
	int method;	 //0未签名，1RSA，2ElGamal
	ZZ r;  //挑战
	ZZ aT, bT;		//秘密指数及其公开值

	/*TA用RSA签名方案颁布的证书*/
	struct RCert
	{
		ZZ id_ver_s;
		ZZ id_ver;
		ZZ s;
		ZZ bT;
	} rcert;
	
	/*TA用ElGamal签名方案颁布的证书*/
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
/*可信权威机构*/
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

	/*生成用户的ID*/
	void GenID(user &p)
	{
		p.id = str_to_zz(p.name);
	}
	/*生成用户证书*/
	void GenCert(user &p,int bits)
	{
		ZZ s, c,t;
		if (method == 1)
		{
			/*生成TA 的密钥*/
			rsa.GenKey(bits);

			/*生成证书*/
			t = connect(p.id, p.rsa.rk.ver);
			rsa.RSA_sig(t, s);		//##################################
			c = connect(t, s);

			/*颁布证书*/
			p.rcert.id_ver_s = c;
			p.rcert.id_ver = t;
			p.rcert.s = s;
			p.rcert.bT = p.bT;
		}
		else
		{
			/*生成TA的密钥*/
			elg.GenKey(bits);

			/*生成证书*/
			t = connect(p.id, p.elg.ek.ver);
			ZZ gamma, delta;
			elg.ElGamal_sig(t, gamma, delta);
			s = connect(gamma, delta);
			c = connect(t, s);

			/*颁布证书*/
			p.ecert.gamma = gamma;
			p.ecert.delta = delta;
			p.ecert.id_ver_s = c;
			p.ecert.id_ver = t;
			p.ecert.s = s;
			p.ecert.bT = p.bT;
		}
	}
	/*生成证书文件*/
	void GenCertFile(user p)
	{
		ofstream fout;
		fout.open("Certification.txt", ios::out);
		if (!fout.is_open())
		{
			cout << "打开文件失败" << endl;
			return;
		}
		if (method == 1)
			fout << p.rcert.id_ver_s;
		else
			fout << p.ecert.id_ver_s;
		fout.close();
	}
	/*证书验证方法*/
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

/*RSA签名测试*/
void RSA_test()
{
	user Alice("Alice");
	Alice.method = 1;

	ZZ x, y;

	cout << "请输入消息：" << endl;
	cin >> x;

	/*Alice 生成密钥*/
	Alice.rsa.GenKey();
	cout << endl<<"已生成Alice的密钥，请按Enter键继续……" << endl;
	getchar();
	getchar();

	/*Alice 签名*/
	Alice.rsa.RSA_sig(x, y);
	cout << "Alice已完成签名"<<endl;
	cout <<endl<< y << endl;
	cout << endl<<"请按Enter键继续……" << endl;
	getchar();

}

/*ElGamal签名测试*/
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

	cout << "请输入消息：" << endl;
	cin >> x;

	/*Alice 生成密钥*/
	Alice.elg.GenKey();
	cout << endl<<"已生成Alice的密钥，请按Enter键继续……"<<endl;
	getchar();
	getchar();

	/*Alice 签名*/
	Alice.elg.ElGamal_sig(x, y.gamma, y.delta);
	cout <<endl<< "Alice已完成签名" << endl;
	cout << endl<< y.gamma<<'  '<<y.delta << endl;
	cout << endl << "请按Enter键继续……" << endl;
	getchar();

}

/*证书的颁布和验证测试*/
void Cert_test()
{
	const int bits = 64;

	user Alice("Alice");
	user Bob("Bob");
	ta TA;

	int choice;

	cout << "请选择TA的签名方式（1.RSA  2.ElGamal）：";
	cin >> choice;
	TA.method = choice;
	cout << endl;

	cout << "您好Alice，请选择需要的签名方式（1.RSA  2.ElGamal）：";
	cin >> choice;
	Alice.method = choice;
	cout << endl;

	/*Alice 生成密钥*/
	if (Alice.method == 1)
		Alice.rsa.GenKey(bits);
	else
		Alice.elg.GenKey(bits);
	cout << "已生成Alice的密钥，请按Enter键继续……";
	getchar();
	getchar();

	/*TA 生成 Alice ID*/
	TA.GenID(Alice);
	cout << "已生成Alice的ID，请按Enter键继续……";
	getchar();

	/*TA 向 Alice 颁布证书， 并生成文件*/
	TA.GenCert(Alice,4*bits);		//##################################
	cout << "已生成Alice的证书，请按Enter键继续……";
	getchar();
	TA.GenCertFile(Alice);
	cout << "已生成Alice的证书文件，请按Enter键继续……";
	getchar();

	/*Bob 利用 TA 来验证 Alice 的证书*/
	bool TorF=TA.Cert_ver(Alice);
	cout << "Bob已验证Alice的证书，结果为：";
	if (TorF)
		cout << "True.    ";
	else
		cout << "False.    ";
	cout<<"请按Enter键继续……";
	getchar();
	cout << endl;
	return;
}

/*交互认证测试*/
void Interact_Cert_test()
{
	const int bits = 64;

	user Alice("Alice");
	user Bob("Bob");
	ta TA;

	int choice;

	cout << "请选择TA的签名方案（1. RSA  2. ElGamal）：";
	cin >> choice;
	TA.method = choice;

	cout << "请选择Alice的签名方案（1. RSA  2. ElGamal）：";
	cin >> choice;
	Alice.method = choice;

	cout << "请选择Bob的签名方案（1. RSA  2. ElGamal）：";
	cin >> choice;
	Bob.method = choice;

	cout << endl;

	if (Alice.method == 1)
		Alice.rsa.GenKey(bits);
	else
		Alice.elg.GenKey(bits);
	cout << "已生成Alice的密钥，请按Enter键继续……";
	getchar();
	getchar();

	/*Alice 生成证书 Cert(Alice)*/
	TA.GenID(Alice);
	cout << "已生成Alice的ID，请按Enter键继续……";
	getchar();
	TA.GenCert(Alice,4*bits);
	cout << "已生成Alice的证书，请按Enter键继续……";
	getchar();

	/*Bob 验证Alice 的公钥，并利用Alice 的公钥验证*/
	bool TorFa = TA.Cert_ver(Alice);
	cout << "Bob已验证Alice的公钥，结果为：";
	if (TorFa)
		cout << "True.    ";
	else
		cout << "False.    ";
	cout << "请按Enter键继续……";
	getchar();

	if (Bob.method == 1)
		Bob.rsa.GenKey(bits);
	else
		Bob.elg.GenKey(bits);
	cout << "已生成Bob的密钥，请按Enter键继续……";
	getchar();

	/*Bob 生成证书 Cert(Bob)*/
	TA.GenID(Bob);
	cout << "已生成Bob的ID，请按Enter键继续……";
	getchar();
	TA.GenCert(Bob,4*bits);
	cout << "已生成Bob的证书，请按Enter键继续……";
	getchar();

	

	/*Alice 验证Bob 的公钥，并利用Bob 的公钥验证*/
	bool TorFb = TA.Cert_ver(Bob);
	cout << "Alice已验证Bob的公钥，结果为：";
	if (TorFb)
		cout << "True.    ";
	else
		cout << "False.    ";
	cout << "请按Enter键继续……";
	getchar();

	SetSeed(conv<ZZ>(static_cast<long>(time(nullptr))));
	/*Bob 生成随机挑战 r1*/
	ZZ r1;
	RandomBits(r1, bits/8);
	Bob.r = r1;
	cout << "已生成Bob的随机挑战r1，请按Enter键继续……";
	getchar();

	/*Alice 生成随机挑战 r2*/
	ZZ r2;
	RandomBits(r2, bits/8);
	Alice.r = r2;
	cout << "已生成Alice的随机挑战r2，请按Enter键继续……";
	getchar();

	/*Alice 计算 y1*/
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
	cout << "Alice已计算y1，请按Enter键继续……";
	getchar();


	bool AorR;
	if (Alice.method == 1)
		AorR=Alice.rsa.RSA_ver(connect(connect(Bob.id, Bob.r), Alice.r), y1.rsa_y);
	else
		AorR = Alice.elg.ElGamal_ver(connect(connect(Bob.id, Bob.r), Alice.r), y1.elg_gamma, y1.elg_delta);
	if (AorR)
		cout << "Bob 接受，请按Enter键继续……";
	else
		cout << "Bob 拒绝，请按Enter键继续……";
	getchar();

	/*Bob 计算 y2*/
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
	cout << "Bob已计算y2，请按Enter键继续……";
	getchar();

	if (Bob.method == 1)
		AorR=Bob.rsa.RSA_ver(connect(Alice.id, Alice.r), y2.rsa_y);
	else
		AorR=Bob.elg.ElGamal_ver(connect(Alice.id, Alice.r), y2.elg_gamma, y2.elg_delta);
	if (AorR)
		cout << "Alice 接受，请按Enter键继续……";
	else
		cout << "Alice 拒绝，请按Enter键继续……";
	getchar();
	cout << endl;
}

/*MTI/A0密钥协商方案测试*/
void MTI_A0_test()
{

	const int bits = 1024;
	ZZ p, a, n;

	user U("U"), V("V");
	ta TA;
	int choice;
	
	cout << "请选择TA的签名方案（1. RSA  2. ElGamal）：";
	cin >> choice;
	TA.method = choice;

	cout << "请选择U的签名方案（1. RSA  2. ElGamal）：";
	cin >> choice;
	U.method = choice;

	cout << "请选择V的签名方案（1. RSA  2. ElGamal）：";
	cin >> choice;
	V.method = choice;

	SetSeed(conv<ZZ>(static_cast<long>(time(nullptr))));
	
	/*	生成素数p*/
	ZZ q;
	do
	{
		GenGermainPrime(q, bits - 1);
		p = 2 * q + 1;
	} while (!is_prime(p));

	/*	生成本原元a	*/
	get_pri_root(a, p, q);
	n = p - 1;

	cout << endl;
	if (U.method == 1)
		U.rsa.GenKey(bits/8);
	else
		U.elg.GenKey(bits/8);
	cout << "已生成U的密钥，请按Enter键继续……";
	getchar();
	getchar();

	if (V.method == 1)
		V.rsa.GenKey(bits/8);
	else
		V.elg.GenKey(bits/8);
	cout << "已生成V的密钥，请按Enter键继续……";
	getchar();

	SetSeed(conv<ZZ>(static_cast<long>(time(nullptr))));
	/*U选取秘密指数，并计算公开值，生成证书*/
	RandomBnd(U.aT, n);
	PowerMod(U.bT, a, U.aT, p);
	cout << "已生成U的秘密指数及其公开值，请按Enter键继续……";
	getchar();

	TA.GenID(U);
	TA.GenCert(U,bits/2);
	cout << "已生成U的证书，请按Enter键继续……";
	getchar();

	/*V选取秘密指数，并计算公开值，生成证书*/
	RandomBnd(V.aT, n);
	PowerMod(V.bT, a, V.aT, p);
	cout << "已生成V的秘密指数及其公开值，请按Enter键继续……";
	getchar();

	TA.GenID(V);
	TA.GenCert(V,bits/2);
	cout << "已生成V的证书，请按Enter键继续……";
	getchar();

	/*U选取随机值rU，并计算sU*/
	ZZ rU, sU;
	RandomBnd(rU, n);
	PowerMod(sU, a, rU, p); 
	cout << "U已选取随机值rU并计算出了对应的sU，请按Enter键继续……";
	getchar();

	/*V选取随机值rV，并计算sV*/
	ZZ rV, sV;
	RandomBnd(rV, n);
	PowerMod(sV, a, rV, p);
	cout << "V已选取随机值rV并计算出了对应的sV，请按Enter键继续……";
	getchar();


	ZZ ku1, ku2;
	/*U计算会话密钥*/
	ZZ KU;
	PowerMod(ku1, sV, U.aT, p);
	if (TA.method == 1)
		PowerMod(ku2, V.rcert.bT, rU, p);
	else
		PowerMod(ku2, V.ecert.bT, rU, p);
	MulMod(KU, ku1, ku2, p);
	cout << "U已计算出会话密钥KU，请按Enter键继续……";
	getchar();

	ZZ kv1, kv2;
	/*V计算会话密钥*/
	ZZ KV;
	PowerMod(kv1, sU, V.aT, p);
	if (TA.method == 1)
		PowerMod(kv2, U.rcert.bT, rV, p);
	else
		PowerMod(kv2, U.ecert.bT, rV, p);
	MulMod(KV, kv1, kv2, p);
	cout << "V已计算出会话密钥KV，请按Enter键继续……";
	getchar();

	cout << endl;
	/*判断会话密钥是否相同*/
	cout << "已对比U和V所计算出的会话密钥，结果为：";
	if (compare(KU, KV) == 0)
	{
		cout << "相同.    " << endl;
		cout << endl;
		cout << KU << endl;
		cout << "请按Enter键继续……";
	}
	else
	{
		cout << "不同.    ";
		cout << "请按Enter键继续……";
	}
	
	
	getchar();
	cout << endl;
}

int main()
{
	system("cls");
	
	while (1)
	{
		cout << "请选择：" << endl;
		cout << "1. RSA签名方案" << endl;
		cout << "2. ElGamal签名方案" << endl;
		cout << "3. 证书的颁布和验证" << endl;
		cout << "4. 交互认证" << endl;
		cout << "5. MTI/A0密钥协商方案" << endl;
		cout << endl;
		cout << "请输入对应的序号，或按Q/q退出" << endl << endl;
		char choice;
		cin >> choice;

		system("cls");

		/*选择退出*/
		if (choice == 'Q' || choice == 'q')
		{
			cout << "已退出" << endl << endl;
			break;
		}

		/*选择1*/
		else if (choice == '1')
		{
			cout << "1. RSA签名方案" << endl << endl;
			RSA_test();
		}

		/*选择2*/
		else if (choice == '2')
		{
			cout << "2. ElGamal签名方案" << endl << endl;
			ElGamal_test();
		}

		/*选择3*/
		else if (choice == '3')
		{
			cout << "3. 证书的颁布和验证" << endl << endl;
			Cert_test();
		}

		/*选择4*/
		else if (choice == '4')
		{
			cout << "4. 交互认证" << endl << endl;
			Interact_Cert_test();
		}

		/*选择5*/
		else if (choice == '5')
		{
			cout << "5. MTI/A0密钥协商方案" << endl << endl;
			MTI_A0_test();
		}

		/*选择无效*/
		else
			continue;

		cout << "请按Enter键返回主菜单……";
		getchar();
		system("cls");

	}

	return 0;
}