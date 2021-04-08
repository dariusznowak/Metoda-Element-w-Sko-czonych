//etap tworzenia macierzy [H] globalnej
#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>

using namespace std;

void wypisz_macierz(vector<vector<double>> macierz)
{
	for (int i = 0; i < macierz.size(); i++)
	{
		for (int j = 0; j < macierz[i].size(); j++)
		{
			cout << fixed << showpoint << setprecision(3);
			cout << macierz[i][j] << "\t";
		}
		cout << endl;
	}
	cout << endl;
}

void wypisz_wektor(vector<double> wektor)
{
	for (int i = 0; i < wektor.size(); i++)
	{
		cout << fixed << showpoint << setprecision(1);
		cout << wektor[i] << "   ";
		if (i == wektor.size() - 1) cout << endl << endl;
	}
	cout << endl;
}


/******************************STRUKTURY WEZLU i ELEMENTU SKONCZONEGO oraz FUNKCJA GENERUJĄCA SIATKĘ ELEMENTÓW SKONCZONYCH**************************************************************************************************/
struct node
{//struktura pojedynczego wezlu, zawierajaca wspolrzedne (x,y) danego wezlu
	double x, y;
	int BC = 0; // wartosc domyslna 0; jezeli bd trzeba to zostanie zmieniona na 1
	double t_0; // temperatura poczatkowa - podaje ją puki co tutaj, bo nie wiem na wszystkich wezłach mieliśmy daną taką samą wartość t_0		
	double t_1; //temperatura obliczona w danym kroku
};

struct element
{//struktura pojedynczego elementu skonczonego, zawierajaca ID czterech wezłów, z których się składa
	double K;
	int id[4];
	vector<vector<double>> Hl; // lokalna macierz [H]
	vector<vector<double>> Cl; // lokalna macierz [C]
	vector<vector<double>> jacobian2pkt[4]; // jacobian
	vector<vector<double>> jacobianOdwrocony2pkt[4]; // jacobian odwrocony
	vector<vector<double>> jacobian3pkt[9]; // jacobian
	vector<vector<double>> jacobianOdwrocony3pkt[9]; // jacobian odwrocony
	vector<vector<double>> Hbc; // lokalna macierz [Hbc] (H z war. brzegowym)
	vector<double> Pl;
};

struct globalData
{//struktura zawierajaca dane globalne, czyli informacje o CALEJ SIATCE
	double t_0; //temp. poczatkowa
	double sim_time; //czas symulacji
	double dT; //d"tau" krok czasowy
	double t_alfa; //temp. otoczenia, uzywana w obliczaniu macierzy p dla kazdego boku danego elementu skonczonego
	double alfa; //wsp. przenikania ciepla (do macierzy Hbc z warunkiem brzegowym)
	double h; //wysokosc siatki
	double w; //szerokosc siatki
	double nH; //l. wezlow po wysokosci
	double nW; //l. wezlow po szerokosci
	double cp; //cieplo wlasciwe
	double K; //wsp. przewodzenia ciepla
	double ro; //gestosc
	int sc; //ktory schemat calkowania
	double nE = (nH - 1) * (nW - 1); //l. elementow
	double nN = nH * nW; //l. wezlow
};

void setElementsAndNodes(element* Elem, node* node, globalData dane)
{
	double krokWBok = dane.w / (dane.nW - 1);
	double krokWGore = dane.h / (dane.nH - 1);
	const int temp = (int)dane.nH - 1;

	/*WYZNACZENIE TABLIC INDEKSOW*/
	int i = 1; //licznik wezlow
	int j = 0; //wskaznik ktore id

	for (int nrElem = 1; nrElem <= dane.nE; nrElem++)
	{
		if (nrElem != 1 && (nrElem % temp) == 1) i++; //sprawdzamy czy jestesmy na dole nowej kolumny - jesli tak to trzeba zwiekszyc licznik wezlow o 1

		Elem[nrElem - 1].id[j] = i;
		Elem[nrElem - 1].id[j + 1] = i + dane.nH;
		Elem[nrElem - 1].id[j + 2] = Elem[nrElem - 1].id[j + 1] + 1;
		Elem[nrElem - 1].id[j + 3] = Elem[nrElem - 1].id[j] + 1;

		i++;
	}

	/***********wyznaczanie wspolrzednych wezlow**********/

	i = 0;
	for (double x = 0.0; x <= dane.w; x += krokWBok)
	{
		for (double y = 0.0; y <= dane.h; y += krokWGore)
		{
			//@@@ tymczasowy fragment @@@// ogrzewanie prawej i dolnej sciany
			//if (y == 0.0 || x == dane.w) node[i].BC = 1;
			if (y == 0.0 || x == dane.w || x == 0 || y == dane.h) node[i].BC = 1;

			//@@@  @@@//

			node[i].x = x;
			node[i].y = y;
			//cout << "x = " << x << "    " << "y = " << y << "  BC = " << node[i].BC << endl;
			i++;
		}
	}
	cout << endl;

}

/***************************************************************************************************************************************************************/
struct SoE
{
	vector<vector<double>> HG;
	vector<vector<double>> CG;
	vector<double> PG;
};

struct Elem4
{	//punkty calkowania w ukladzie "ksi-eta"
	double ksi2pkt[4] = { -1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3), -1 / sqrt(3) };	                 							  // pkt calkowania
	double eta2pkt[4] = { -1 / sqrt(3), -1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3) };	 				     						  // 
	double ksi3pkt[9] = { -sqrt(3 / 5.0), 0, sqrt(3 / 5.0), -sqrt(3 / 5.0), 0, sqrt(3 / 5.0), -sqrt(3 / 5.0), 0, sqrt(3 / 5.0) }; // 
	double eta3pkt[9] = { -sqrt(3 / 5.0), -sqrt(3 / 5.0), -sqrt(3 / 5.0), 0, 0, 0, sqrt(3 / 5.0), sqrt(3 / 5.0), sqrt(3 / 5.0) }; // 

	double wagi3pkt[3] = { 5 / 9.0, 8 / 9.0, 5 / 9.0 };

	double ksiPoPowierzchni_2pkt[8] = { -1 / sqrt(3), 1 / sqrt(3), 1, 1, 1 / sqrt(3), -1 / sqrt(3), -1, -1 }; // pkt calkowania po powierzchni (schemat 2pkt), co BARDZO WAZNE: <------------
	double etaPoPowierzchni_2pkt[8] = { -1, -1, -1 / sqrt(3), 1 / sqrt(3), 1, 1, 1 / sqrt(3), -1 / sqrt(3) }; // kolejno: 2 na dół, 2 na prawy bok, 2 na góre i 2 na lewo
																										 //--------------------------------------------------------------

	double ksiPoPowierzchni_3pkt[12] = { -sqrt(3 / 5.0), 0, sqrt(3 / 5.0), 1, 1, 1, sqrt(3 / 5.0), 0, -sqrt(3 / 5.0), -1, -1, -1 }; // pkt calkowania po powierzchni (schemat 3pkt) , co BARDZO WAZNE: <------------
	double etaPoPowierzchni_3pkt[12] = { -1, -1, -1,  -sqrt(3 / 5.0), 0, sqrt(3 / 5.0), 1, 1, 1, sqrt(3 / 5.0), 0, -sqrt(3 / 5.0) }; // kolejno: 3 na dół, 3 na prawy bok, 3 na góre i 3 na lewo
																									 //--------------------------------------------------------------
	vector<vector<double>>pochodnePoKsi2pkt, pochodnePoEta2pkt, pochodnePoKsi3pkt, pochodnePoEta3pkt;
	vector<vector<double>> macierzN2pkt, macierzN3pkt, macierzN_poPowierzchni_2pkt, macierzN_poPowierzchni_3pkt; //wartosci funkcji ksztaltu w poszczegolnych punktach calkowania

	//pochodne funkcji ksztaltu dla kazdej wspolrzednej (tutaj po zmiennej "ksi")
	Elem4()
	{
		oblicz_pochodne2pkt();
		oblicz_pochodne3pkt();
		oblicz_wart_f_ksztaltu_2pkt();
		oblicz_wart_f_ksztaltu_3pkt();
		oblicz_wart_f_ksztaltu_na_powierzchni_2pkt();
		oblicz_wart_f_ksztaltu_na_powierzchni_3pkt();
	}

	//pochodne fukcji ksztaltu w danym punkcie (ksi,eta)
	void oblicz_pochodne2pkt()
	{

		for (int nrPunktu = 0; nrPunktu < 4; nrPunktu++)
		{
			vector<double> tempKsi;
			vector<double> tempEta;
			tempKsi.push_back(-0.25 * (1 - eta2pkt[nrPunktu]));
			tempKsi.push_back(0.25 * (1 - eta2pkt[nrPunktu]));
			tempKsi.push_back(0.25 * (1 + eta2pkt[nrPunktu]));
			tempKsi.push_back(-0.25 * (1 + eta2pkt[nrPunktu]));
			pochodnePoKsi2pkt.push_back(tempKsi);

			tempEta.push_back(-0.25 * (1 - ksi2pkt[nrPunktu]));
			tempEta.push_back(-0.25 * (1 + ksi2pkt[nrPunktu]));
			tempEta.push_back(0.25 * (1 + ksi2pkt[nrPunktu]));
			tempEta.push_back(0.25 * (1 - ksi2pkt[nrPunktu]));
			pochodnePoEta2pkt.push_back(tempEta);
		}
	}

	void oblicz_pochodne3pkt()
	{
		for (int nrPunktu = 0; nrPunktu < 9; nrPunktu++)
		{
			vector<double> tempKsi;
			vector<double> tempEta;
			tempKsi.push_back(-0.25 * (1 - eta3pkt[nrPunktu]));
			tempKsi.push_back(0.25 * (1 - eta3pkt[nrPunktu]));
			tempKsi.push_back(0.25 * (1 + eta3pkt[nrPunktu]));
			tempKsi.push_back(-0.25 * (1 + eta3pkt[nrPunktu]));
			pochodnePoKsi3pkt.push_back(tempKsi);

			tempEta.push_back(-0.25 * (1 - ksi3pkt[nrPunktu]));
			tempEta.push_back(-0.25 * (1 + ksi3pkt[nrPunktu]));
			tempEta.push_back(0.25 * (1 + ksi3pkt[nrPunktu]));
			tempEta.push_back(0.25 * (1 - ksi3pkt[nrPunktu]));
			pochodnePoEta3pkt.push_back(tempEta);
		}
	}

	void oblicz_wart_f_ksztaltu_2pkt()
	{
		for (int nrPunktu = 0; nrPunktu < 4; nrPunktu++)
		{
			vector<double> temp;
			temp.push_back(0.25 * (1 - ksi2pkt[nrPunktu]) * (1 - eta2pkt[nrPunktu]));
			temp.push_back(0.25 * (1 + ksi2pkt[nrPunktu]) * (1 - eta2pkt[nrPunktu]));
			temp.push_back(0.25 * (1 + ksi2pkt[nrPunktu]) * (1 + eta2pkt[nrPunktu]));
			temp.push_back(0.25 * (1 - ksi2pkt[nrPunktu]) * (1 + eta2pkt[nrPunktu]));
			macierzN2pkt.push_back(temp);
		}
		//wypisz_macierz(macierzN2pkt);
	}

	void oblicz_wart_f_ksztaltu_3pkt()
	{
		for (int nrPunktu = 0; nrPunktu < 9; nrPunktu++)
		{
			vector<double> temp;
			temp.push_back(0.25 * (1 - ksi3pkt[nrPunktu]) * (1 - eta3pkt[nrPunktu]));
			temp.push_back(0.25 * (1 + ksi3pkt[nrPunktu]) * (1 - eta3pkt[nrPunktu]));
			temp.push_back(0.25 * (1 + ksi3pkt[nrPunktu]) * (1 + eta3pkt[nrPunktu]));
			temp.push_back(0.25 * (1 - ksi3pkt[nrPunktu]) * (1 + eta3pkt[nrPunktu]));
			macierzN3pkt.push_back(temp);
		}
		//wypisz_macierz(macierzN3pkt);

	}

	void oblicz_wart_f_ksztaltu_na_powierzchni_2pkt()
	{
		for (int nrPunktu = 0; nrPunktu < 8; nrPunktu++)
		{
			vector<double> temp;
			temp.push_back(0.25 * (1 - ksiPoPowierzchni_2pkt[nrPunktu]) * (1 - etaPoPowierzchni_2pkt[nrPunktu]));
			temp.push_back(0.25 * (1 + ksiPoPowierzchni_2pkt[nrPunktu]) * (1 - etaPoPowierzchni_2pkt[nrPunktu]));
			temp.push_back(0.25 * (1 + ksiPoPowierzchni_2pkt[nrPunktu]) * (1 + etaPoPowierzchni_2pkt[nrPunktu]));
			temp.push_back(0.25 * (1 - ksiPoPowierzchni_2pkt[nrPunktu]) * (1 + etaPoPowierzchni_2pkt[nrPunktu]));
			macierzN_poPowierzchni_2pkt.push_back(temp);
		}
	}

	void oblicz_wart_f_ksztaltu_na_powierzchni_3pkt()
	{
		for (int nrPunktu = 0; nrPunktu < 12; nrPunktu++)
		{
			vector<double> temp;
			temp.push_back(0.25 * (1 - ksiPoPowierzchni_3pkt[nrPunktu]) * (1 - etaPoPowierzchni_3pkt[nrPunktu]));
			temp.push_back(0.25 * (1 + ksiPoPowierzchni_3pkt[nrPunktu]) * (1 - etaPoPowierzchni_3pkt[nrPunktu]));
			temp.push_back(0.25 * (1 + ksiPoPowierzchni_3pkt[nrPunktu]) * (1 + etaPoPowierzchni_3pkt[nrPunktu]));
			temp.push_back(0.25 * (1 - ksiPoPowierzchni_3pkt[nrPunktu]) * (1 + etaPoPowierzchni_3pkt[nrPunktu]));
			macierzN_poPowierzchni_3pkt.push_back(temp);
		}
		//wypisz_macierz(macierzN_poPowierzchni_3pkt);
	}

	void wypisz_pochodne2pkt()
	{
		cout << "Pochodne funkcji ksztaltu po \"ksi\"" << endl;
		wypisz_macierz(pochodnePoKsi2pkt);
		cout << endl << "Pochodne funkcji ksztaltu po \"eta\"" << endl;
		wypisz_macierz(pochodnePoEta2pkt);
	}
	void wypisz_pochodne3pkt()
	{
		cout << "Pochodne funkcji ksztaltu po \"ksi\"" << endl;
		wypisz_macierz(pochodnePoKsi3pkt);
		cout << endl << "Pochodne funkcji ksztaltu po \"eta\"" << endl;
		wypisz_macierz(pochodnePoEta3pkt);
	}
};

vector<vector<double>> oblicz_jacobian(Elem4 elem, vector<vector<double>> punktyXY, int nrPktCalkowania, const int ktorySchematCalkowania)//funkcja zwracajaca macierz - jacobian; przyjmuje tabelke punktow calkowania (x,y)
{//sa 4 jacobiany, bo sa 4 wezly (trzeba dorobic dla kolejnych punktow ze struktury elem4)
	vector<vector<double>> jacobian;
	vector<double> temp;
	nrPktCalkowania--;

	if (ktorySchematCalkowania == 2)
	{
		temp.push_back((elem.pochodnePoKsi2pkt[nrPktCalkowania][0]) * punktyXY[0][0] + (elem.pochodnePoKsi2pkt[nrPktCalkowania][1]) * punktyXY[1][0] + (elem.pochodnePoKsi2pkt[nrPktCalkowania][2]) * punktyXY[2][0] + (elem.pochodnePoKsi2pkt[nrPktCalkowania][3]) * punktyXY[3][0]);
		temp.push_back((elem.pochodnePoKsi2pkt[nrPktCalkowania][0]) * punktyXY[0][1] + (elem.pochodnePoKsi2pkt[nrPktCalkowania][1]) * punktyXY[1][1] + (elem.pochodnePoKsi2pkt[nrPktCalkowania][2]) * punktyXY[2][1] + (elem.pochodnePoKsi2pkt[nrPktCalkowania][3]) * punktyXY[3][1]);

		jacobian.push_back(temp);
		temp.clear();

		temp.push_back((elem.pochodnePoEta2pkt[nrPktCalkowania][0]) * punktyXY[0][0] + (elem.pochodnePoEta2pkt[nrPktCalkowania][1]) * punktyXY[1][0] + (elem.pochodnePoEta2pkt[nrPktCalkowania][2]) * punktyXY[2][0] + (elem.pochodnePoEta2pkt[nrPktCalkowania][3]) * punktyXY[3][0]);
		temp.push_back((elem.pochodnePoEta2pkt[nrPktCalkowania][0]) * punktyXY[0][1] + (elem.pochodnePoEta2pkt[nrPktCalkowania][1]) * punktyXY[1][1] + (elem.pochodnePoEta2pkt[nrPktCalkowania][2]) * punktyXY[2][1] + (elem.pochodnePoEta2pkt[nrPktCalkowania][3]) * punktyXY[3][1]);

		jacobian.push_back(temp);
		temp.clear();

		//cout << "Jacobian w punkcie (" << elem.ksi2pkt[nrPktCalkowania] << ", " << elem.eta2pkt[nrPktCalkowania] << "):" << endl;
		//wypisz_macierz(jacobian);

		return jacobian;
	}
	else if (ktorySchematCalkowania == 3)
	{
		temp.push_back((elem.pochodnePoKsi3pkt[nrPktCalkowania][0]) * punktyXY[0][0] + (elem.pochodnePoKsi3pkt[nrPktCalkowania][1]) * punktyXY[1][0] + (elem.pochodnePoKsi3pkt[nrPktCalkowania][2]) * punktyXY[2][0] + (elem.pochodnePoKsi3pkt[nrPktCalkowania][3]) * punktyXY[3][0]);
		temp.push_back((elem.pochodnePoKsi3pkt[nrPktCalkowania][0]) * punktyXY[0][1] + (elem.pochodnePoKsi3pkt[nrPktCalkowania][1]) * punktyXY[1][1] + (elem.pochodnePoKsi3pkt[nrPktCalkowania][2]) * punktyXY[2][1] + (elem.pochodnePoKsi3pkt[nrPktCalkowania][3]) * punktyXY[3][1]);

		jacobian.push_back(temp);
		temp.clear();

		temp.push_back((elem.pochodnePoEta3pkt[nrPktCalkowania][0]) * punktyXY[0][0] + (elem.pochodnePoEta3pkt[nrPktCalkowania][1]) * punktyXY[1][0] + (elem.pochodnePoEta3pkt[nrPktCalkowania][2]) * punktyXY[2][0] + (elem.pochodnePoEta3pkt[nrPktCalkowania][3]) * punktyXY[3][0]);
		temp.push_back((elem.pochodnePoEta3pkt[nrPktCalkowania][0]) * punktyXY[0][1] + (elem.pochodnePoEta3pkt[nrPktCalkowania][1]) * punktyXY[1][1] + (elem.pochodnePoEta3pkt[nrPktCalkowania][2]) * punktyXY[2][1] + (elem.pochodnePoEta3pkt[nrPktCalkowania][3]) * punktyXY[3][1]);

		jacobian.push_back(temp);
		temp.clear();

		//cout << "Jacobian w punkcie (" << elem.ksi3pkt[nrPktCalkowania] << ", " << elem.eta3pkt[nrPktCalkowania] << "):" << endl;
		//wypisz_macierz(jacobian);

		return jacobian;
	}
}

vector<vector<double>> odwroc_macierz_jacobiego(vector<vector<double>> jacobian)
{
	vector<vector<double>> temp = jacobian;
	double det_odw_jacobian = 1 / ((jacobian[0][0] * jacobian[1][1]) - (jacobian[0][1] * jacobian[0][1]));

	temp[0][0] = jacobian[1][1] * det_odw_jacobian;
	temp[1][1] = jacobian[0][0] * det_odw_jacobian;
	temp[0][1] = jacobian[0][1] * det_odw_jacobian * (-1);
	temp[1][0] = jacobian[1][0] * det_odw_jacobian * (-1);

	return temp;
}

vector<vector<double>> dNidX_dNidY(vector<vector<double>> jacobian_odw, Elem4 element, int nrPktCalkowania, const int ktorySchematCalkowania)
{
	nrPktCalkowania--;
	vector<vector<double>> result;
	vector<double> temp;

	if (ktorySchematCalkowania == 2)
	{
		for (int nrFunkcjiKsztaltu = 0; nrFunkcjiKsztaltu <= 3; nrFunkcjiKsztaltu++)
			temp.push_back((jacobian_odw[0][0] * element.pochodnePoKsi2pkt[nrPktCalkowania][nrFunkcjiKsztaltu]) + (jacobian_odw[0][1]) * element.pochodnePoEta2pkt[nrPktCalkowania][nrFunkcjiKsztaltu]);
		result.push_back(temp);
		temp.clear();

		for (int nrFunkcjiKsztaltu = 0; nrFunkcjiKsztaltu <= 3; nrFunkcjiKsztaltu++)
			temp.push_back((jacobian_odw[1][0] * element.pochodnePoKsi2pkt[nrPktCalkowania][nrFunkcjiKsztaltu]) + (jacobian_odw[1][1]) * element.pochodnePoEta2pkt[nrPktCalkowania][nrFunkcjiKsztaltu]);
		result.push_back(temp);
		return result;
	}
	else if (ktorySchematCalkowania == 3)
	{
		for (int nrFunkcjiKsztaltu = 0; nrFunkcjiKsztaltu <= 3; nrFunkcjiKsztaltu++)
			temp.push_back((jacobian_odw[0][0] * element.pochodnePoKsi3pkt[nrPktCalkowania][nrFunkcjiKsztaltu]) + (jacobian_odw[0][1]) * element.pochodnePoEta3pkt[nrPktCalkowania][nrFunkcjiKsztaltu]);
		result.push_back(temp);
		temp.clear();

		for (int nrFunkcjiKsztaltu = 0; nrFunkcjiKsztaltu <= 3; nrFunkcjiKsztaltu++)
			temp.push_back((jacobian_odw[1][0] * element.pochodnePoKsi3pkt[nrPktCalkowania][nrFunkcjiKsztaltu]) + (jacobian_odw[1][1]) * element.pochodnePoEta3pkt[nrPktCalkowania][nrFunkcjiKsztaltu]);
		result.push_back(temp);

		return result;
	}
}

vector<vector<double>> dNdxy_dNdxyT(vector<vector<double>> dNdxyi, vector<vector<double>> jacobian, element elem)
{
	double detJ = (jacobian[0][0] * jacobian[1][1]) - (jacobian[0][1] * jacobian[0][1]);
	vector<double> temp1, temp2;
	vector<vector<double>> mnozenie_dNdx_dNdxT, mnozenie_dNdy_dNdyT, result;

	for (int p = 0; p < 4; p++)
		temp1.push_back(dNdxyi[0][0] * dNdxyi[0][p]);
	mnozenie_dNdx_dNdxT.push_back(temp1);
	temp1.clear();

	for (int p = 0; p < 4; p++)
		temp1.push_back(dNdxyi[0][1] * dNdxyi[0][p]);
	mnozenie_dNdx_dNdxT.push_back(temp1);
	temp1.clear();

	for (int p = 0; p < 4; p++)
		temp1.push_back(dNdxyi[0][2] * dNdxyi[0][p]);
	mnozenie_dNdx_dNdxT.push_back(temp1);
	temp1.clear();

	for (int p = 0; p < 4; p++)
		temp1.push_back(dNdxyi[0][3] * dNdxyi[0][p]);
	mnozenie_dNdx_dNdxT.push_back(temp1);
	temp1.clear();


	//--------------------------------------------------------------------//
	for (int p = 0; p < 4; p++)
		temp2.push_back(dNdxyi[1][0] * dNdxyi[1][p]);
	mnozenie_dNdy_dNdyT.push_back(temp2);
	temp2.clear();

	for (int p = 0; p < 4; p++)
		temp2.push_back(dNdxyi[1][1] * dNdxyi[1][p]);
	mnozenie_dNdy_dNdyT.push_back(temp2);
	temp2.clear();

	for (int p = 0; p < 4; p++)
		temp2.push_back(dNdxyi[1][2] * dNdxyi[1][p]);
	mnozenie_dNdy_dNdyT.push_back(temp2);
	temp2.clear();

	for (int p = 0; p < 4; p++)
		temp2.push_back(dNdxyi[1][3] * dNdxyi[1][p]);
	mnozenie_dNdy_dNdyT.push_back(temp2);
	temp2.clear();

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			mnozenie_dNdx_dNdxT[i][j] += mnozenie_dNdy_dNdyT[i][j];
	}

	for (int i = 0; i < mnozenie_dNdx_dNdxT.size(); i++)
	{//przemnazanie (dN/dx * dN/dxT) oraz (dN/dy * dN/dyT) przez k(t)=30 i deJ=(zalezy dla ktorego pkt, ale raczej rowny 6)
		for (int j = 0; j < mnozenie_dNdx_dNdxT[i].size(); j++)
			mnozenie_dNdx_dNdxT[i][j] = mnozenie_dNdx_dNdxT[i][j] * elem.K * detJ;
	}

	return mnozenie_dNdx_dNdxT;
}

vector<vector<double>> N_razy_Nt_dla_mac_C(Elem4 elementUniwersalny, vector<vector<double>> jacobian, int nrPktCalkowania, const int ktorySchematCalkowania, globalData dane)
{
	double detJ = (jacobian[0][0] * jacobian[1][1]) - (jacobian[0][1] * jacobian[0][1]);
	vector<double> temp;
	vector<vector<double>> result;
	nrPktCalkowania--;

	if (ktorySchematCalkowania == 2)
	{
		for (int p = 0; p < 4; p++)
			temp.push_back(elementUniwersalny.macierzN2pkt[nrPktCalkowania][0] * elementUniwersalny.macierzN2pkt[nrPktCalkowania][p]);
		result.push_back(temp);
		temp.clear();

		for (int p = 0; p < 4; p++)
			temp.push_back(elementUniwersalny.macierzN2pkt[nrPktCalkowania][1] * elementUniwersalny.macierzN2pkt[nrPktCalkowania][p]);
		result.push_back(temp);
		temp.clear();

		for (int p = 0; p < 4; p++)
			temp.push_back(elementUniwersalny.macierzN2pkt[nrPktCalkowania][2] * elementUniwersalny.macierzN2pkt[nrPktCalkowania][p]);
		result.push_back(temp);
		temp.clear();

		for (int p = 0; p < 4; p++)
			temp.push_back(elementUniwersalny.macierzN2pkt[nrPktCalkowania][3] * elementUniwersalny.macierzN2pkt[nrPktCalkowania][p]);
		result.push_back(temp);
		temp.clear();
	}
	else if (ktorySchematCalkowania == 3)
	{
		for (int p = 0; p < 4; p++)
			temp.push_back(elementUniwersalny.macierzN3pkt[nrPktCalkowania][0] * elementUniwersalny.macierzN3pkt[nrPktCalkowania][p]);
		result.push_back(temp);
		temp.clear();

		for (int p = 0; p < 4; p++)
			temp.push_back(elementUniwersalny.macierzN3pkt[nrPktCalkowania][1] * elementUniwersalny.macierzN3pkt[nrPktCalkowania][p]);
		result.push_back(temp);
		temp.clear();

		for (int p = 0; p < 4; p++)
			temp.push_back(elementUniwersalny.macierzN3pkt[nrPktCalkowania][2] * elementUniwersalny.macierzN3pkt[nrPktCalkowania][p]);
		result.push_back(temp);
		temp.clear();

		for (int p = 0; p < 4; p++)
			temp.push_back(elementUniwersalny.macierzN3pkt[nrPktCalkowania][3] * elementUniwersalny.macierzN3pkt[nrPktCalkowania][p]);
		result.push_back(temp);
		temp.clear();
	}


	for (int i = 0; i < result.size(); i++)
	{
		for (int j = 0; j < result[i].size(); j++)
			result[i][j] = result[i][j] * dane.cp * dane.ro * detJ;
	}

	return result;
}

vector<vector<double>> mnozenie_macierzy_przez_stala(vector<vector<double>> macierz, double stala)
{
	for (int i = 0; i < macierz.size(); i++)
		for (int j = 0; j < macierz[i].size(); j++)
		{
			macierz[i][j] = macierz[i][j] * stala;
		}

	return macierz;
}

vector<vector<double>> suma_macierzy(vector<vector<double>> a, vector<vector<double>> b)
{
	for (int i = 0; i < a.size(); i++)
		for (int j = 0; j < a[i].size(); j++)
			a[i][j] += b[i][j];

	return a;
}

vector<vector<double>> sumowanie_macierzy_do_Hbc_3pkt(vector<vector<double>> a, vector<vector<double>> b, vector<vector<double>> c, vector<vector<double>> d)
{
	vector<vector<double>> result;
	vector<double> temp;

	for (int i = 0; i < a.size(); i++)
	{
		for (int j = 0; j < a[i].size(); j++)
			temp.push_back(a[i][j] + b[i][j] + c[i][j] + d[i][j]);
		result.push_back(temp);
		temp.clear();
	}
	return result;
}

vector<vector<double>> sumowanie_macierzy_3pkt_schemat(Elem4 elementUniwersalny, vector<vector<double>> a, vector<vector<double>> b, vector<vector<double>> c, vector<vector<double>> d, vector<vector<double>> e, vector<vector<double>> f, vector<vector<double>> g, vector<vector<double>> h, vector<vector<double>> k)
{
	vector<vector<double>> result;
	vector<double> temp;
	//double w1 = 5 / 9.0;
	double w1 = elementUniwersalny.wagi3pkt[0];
	double w2 = elementUniwersalny.wagi3pkt[1];
	double w3 = elementUniwersalny.wagi3pkt[2];
	//double w2 = 8 / 9.0;
	//double w3 = 5 / 9.0;

	for (int i = 0; i < a.size(); i++)
	{
		for (int j = 0; j < a[i].size(); j++)
			temp.push_back(a[i][j] * w1 * w1 + b[i][j] * w2 * w1 + c[i][j] * w3 * w1 + d[i][j] * w1 * w2 + e[i][j] * w2 * w2 + f[i][j] * w3 * w2 + g[i][j] * w1 * w3 + h[i][j] * w2 * w3 + k[i][j] * w3 * w3);
		result.push_back(temp);
		temp.clear();
	}

	return result;
}

vector<vector<double>> sumowanie_macierzy_do_Hbc_3pkt(vector<vector<double>> a, vector<vector<double>> b)
{
	vector<vector<double>> result;
	vector<double> temp;

	for (int i = 0; i < a.size(); i++)
	{
		for (int j = 0; j < a[i].size(); j++)
			temp.push_back(a[i][j] + b[i][j]);
		result.push_back(temp);
		temp.clear();
	}

	return result;
}

vector<vector<double>> sumowanie_macierzy_do_Hbc_3pkt(vector<vector<double>> a, vector<vector<double>> b, vector<vector<double>> c, Elem4 elementUniwersalny)
{
	vector<vector<double>> result;
	vector<double> temp;

	for (int i = 0; i < a.size(); i++)
	{
		for (int j = 0; j < a[i].size(); j++)
			temp.push_back(a[i][j] * elementUniwersalny.wagi3pkt[0] + b[i][j] * elementUniwersalny.wagi3pkt[1] + c[i][j] * elementUniwersalny.wagi3pkt[2]);
		//temp.push_back(a[i][j] + b[i][j] + c[i][j]);
		result.push_back(temp);
		temp.clear();
	}

	return result;
}

vector<double> sumowanie_macierzy_do_Pl(vector<double> a, vector<double> b)
{
	vector<double> suma;
	for (int i = 0; i < a.size(); i++)
	{
		suma.push_back(a[i] + b[i]);
	}
	return suma;
}

vector<vector<double>> dodaj_macierze(vector<vector<double>> a, vector<vector<double>> b)
{
	for (int i = 0; i < a.size(); i++)
		for (int j = 0; j < a[i].size(); j++)
			a[i][j] += b[i][j];

	return a;
}

void generuj_jacobiany_i_jac_odwrocone(Elem4 elementUniwersalny, element* tabElementowSkonczonych, node* Nd, globalData Dane)
{
	vector<double> temp;
	vector<vector<double>> punktyXY;

	for (int i = 0; i < Dane.nE; i++)
	{
		for (int j = 0; j < 4; j++)
		{ //stworzenie vectora 2d na wzór puntyXY z maina()
			temp.push_back(Nd[tabElementowSkonczonych[i].id[j] - 1].x);  // w skrocie - wczytujemy sobie wspolrzedne wezelow,
			temp.push_back(Nd[tabElementowSkonczonych[i].id[j] - 1].y);	 // tworzacych dany element skonczony
			punktyXY.push_back(temp);
			//cout << "element nr " << jakisElement.id[i] << ": " << "x[" << i << "] = " << punktyXY[i][0] << "   ";
			//cout << "y[" << i << "] = " << punktyXY[i][1] << endl;
			temp.clear();
		}

		if (Dane.sc == 2)
		{
			tabElementowSkonczonych[i].jacobian2pkt[0] = oblicz_jacobian(elementUniwersalny, punktyXY, 1, Dane.sc); // na podstawie obliczonych pochodnych funkcji ksztaltu w obiekcie "element"
			tabElementowSkonczonych[i].jacobian2pkt[1] = oblicz_jacobian(elementUniwersalny, punktyXY, 2, Dane.sc); // obliczamy macierze Jacobiego dla kazdego dla kazdego pkt calkowania;
			tabElementowSkonczonych[i].jacobian2pkt[2] = oblicz_jacobian(elementUniwersalny, punktyXY, 3, Dane.sc); // tak wyglada kazda macierz Jacobiego: [dx/dKsi, dy/dKsi]
			tabElementowSkonczonych[i].jacobian2pkt[3] = oblicz_jacobian(elementUniwersalny, punktyXY, 4, Dane.sc); //									  [dx/dEta, dy/dEta]

			tabElementowSkonczonych[i].jacobianOdwrocony2pkt[0] = odwroc_macierz_jacobiego(tabElementowSkonczonych[i].jacobian2pkt[0]);	//
			tabElementowSkonczonych[i].jacobianOdwrocony2pkt[1] = odwroc_macierz_jacobiego(tabElementowSkonczonych[i].jacobian2pkt[1]);	// odwracanie kazdej
			tabElementowSkonczonych[i].jacobianOdwrocony2pkt[2] = odwroc_macierz_jacobiego(tabElementowSkonczonych[i].jacobian2pkt[2]);	// z 4 macierzy Jacobiego;
			tabElementowSkonczonych[i].jacobianOdwrocony2pkt[3] = odwroc_macierz_jacobiego(tabElementowSkonczonych[i].jacobian2pkt[3]);	//

		}
		else if (Dane.sc == 3)
		{
			tabElementowSkonczonych[i].jacobian3pkt[0] = oblicz_jacobian(elementUniwersalny, punktyXY, 1, Dane.sc); // na podstawie obliczonych pochodnych funkcji ksztaltu w obiekcie "element"
			tabElementowSkonczonych[i].jacobian3pkt[1] = oblicz_jacobian(elementUniwersalny, punktyXY, 2, Dane.sc); // obliczamy macierze Jacobiego dla kazdego dla kazdego pkt calkowania;
			tabElementowSkonczonych[i].jacobian3pkt[2] = oblicz_jacobian(elementUniwersalny, punktyXY, 3, Dane.sc); // tak wyglada kazda macierz Jacobiego: [dx/dKsi, dy/dKsi]
			tabElementowSkonczonych[i].jacobian3pkt[3] = oblicz_jacobian(elementUniwersalny, punktyXY, 4, Dane.sc); //									  [dx/dEta, dy/dEta]
			tabElementowSkonczonych[i].jacobian3pkt[4] = oblicz_jacobian(elementUniwersalny, punktyXY, 5, Dane.sc);
			tabElementowSkonczonych[i].jacobian3pkt[5] = oblicz_jacobian(elementUniwersalny, punktyXY, 6, Dane.sc);
			tabElementowSkonczonych[i].jacobian3pkt[6] = oblicz_jacobian(elementUniwersalny, punktyXY, 7, Dane.sc);
			tabElementowSkonczonych[i].jacobian3pkt[7] = oblicz_jacobian(elementUniwersalny, punktyXY, 8, Dane.sc);
			tabElementowSkonczonych[i].jacobian3pkt[8] = oblicz_jacobian(elementUniwersalny, punktyXY, 9, Dane.sc);

			tabElementowSkonczonych[i].jacobianOdwrocony3pkt[0] = odwroc_macierz_jacobiego(tabElementowSkonczonych[i].jacobian3pkt[0]);	//
			tabElementowSkonczonych[i].jacobianOdwrocony3pkt[1] = odwroc_macierz_jacobiego(tabElementowSkonczonych[i].jacobian3pkt[1]);	// odwracanie kazdej
			tabElementowSkonczonych[i].jacobianOdwrocony3pkt[2] = odwroc_macierz_jacobiego(tabElementowSkonczonych[i].jacobian3pkt[2]);	// z 4 macierzy Jacobiego;
			tabElementowSkonczonych[i].jacobianOdwrocony3pkt[3] = odwroc_macierz_jacobiego(tabElementowSkonczonych[i].jacobian3pkt[3]);	//
			tabElementowSkonczonych[i].jacobianOdwrocony3pkt[4] = odwroc_macierz_jacobiego(tabElementowSkonczonych[i].jacobian3pkt[4]);
			tabElementowSkonczonych[i].jacobianOdwrocony3pkt[5] = odwroc_macierz_jacobiego(tabElementowSkonczonych[i].jacobian3pkt[5]);
			tabElementowSkonczonych[i].jacobianOdwrocony3pkt[6] = odwroc_macierz_jacobiego(tabElementowSkonczonych[i].jacobian3pkt[6]);
			tabElementowSkonczonych[i].jacobianOdwrocony3pkt[7] = odwroc_macierz_jacobiego(tabElementowSkonczonych[i].jacobian3pkt[7]);
			tabElementowSkonczonych[i].jacobianOdwrocony3pkt[8] = odwroc_macierz_jacobiego(tabElementowSkonczonych[i].jacobian3pkt[8]);

			//wypisz_macierz(tabElementowSkonczonych[i].jacobianOdwrocony3pkt[0]);
			//wypisz_macierz(tabElementowSkonczonych[i].jacobianOdwrocony3pkt[1]);
			//wypisz_macierz(tabElementowSkonczonych[i].jacobianOdwrocony3pkt[2]);
			//wypisz_macierz(tabElementowSkonczonych[i].jacobianOdwrocony3pkt[3]);
			//wypisz_macierz(tabElementowSkonczonych[i].jacobianOdwrocony3pkt[4]);
			//wypisz_macierz(tabElementowSkonczonych[i].jacobianOdwrocony3pkt[5]);
			//wypisz_macierz(tabElementowSkonczonych[i].jacobianOdwrocony3pkt[6]);
			//wypisz_macierz(tabElementowSkonczonych[i].jacobianOdwrocony3pkt[7]);
			//wypisz_macierz(tabElementowSkonczonych[i].jacobianOdwrocony3pkt[8]);
		}
	}
	punktyXY.clear();
}

vector<vector<double>> generuj_macierz_H_lokalna(Elem4 elementUniwersalny, element elementSkonczony, node* Nd, globalData Dane) //generowanie macierzy [H] lokalnej - czyli dla danego elementu
{
	//cout << "generowanie H lokalne" << endl;

	if (Dane.sc == 2)
	{
		vector<vector<double>> dNdX_dNdY1 = dNidX_dNidY(elementSkonczony.jacobianOdwrocony2pkt[0], elementUniwersalny, 1, Dane.sc); // obliczanie dla kazdego punktu (x,y) pochodnych dNi/dx i dNi/dy;
		vector<vector<double>> dNdX_dNdY2 = dNidX_dNidY(elementSkonczony.jacobianOdwrocony2pkt[1], elementUniwersalny, 2, Dane.sc); // wyniki sa zapisane jako macierze 2x8, gdzie pierwszy wiersz to:
		vector<vector<double>> dNdX_dNdY3 = dNidX_dNidY(elementSkonczony.jacobianOdwrocony2pkt[2], elementUniwersalny, 3, Dane.sc); // [dN1/dx, dN2/dx, dN3/dx, dN4/dx] 
		vector<vector<double>> dNdX_dNdY4 = dNidX_dNidY(elementSkonczony.jacobianOdwrocony2pkt[3], elementUniwersalny, 4, Dane.sc); // [dN1/dy, dN2/dy, dN3/dy, dN4/dy]

		vector<vector<double>> dNdxy1Ndxy1T = dNdxy_dNdxyT(dNdX_dNdY1, elementSkonczony.jacobian2pkt[0], elementSkonczony); // obliczenie sum: dN/dx * dN/dxT + dN/dy * dN/dyT, pomnozonych przez		//podajemy jeszcze "jakiselement",
		vector<vector<double>> dNdxy2Ndxy2T = dNdxy_dNdxyT(dNdX_dNdY2, elementSkonczony.jacobian2pkt[1], elementSkonczony); // wyznaczniki jacobianow w poszczegolnych punktach oraz k(t)=30;			//który zawiera wsp. K = 30
		vector<vector<double>> dNdxy3Ndxy3T = dNdxy_dNdxyT(dNdX_dNdY3, elementSkonczony.jacobian2pkt[2], elementSkonczony); // wzor: [H] = f(x,y) * w1 * w2 * det[J], 								    //
		vector<vector<double>> dNdxy4Ndxy4T = dNdxy_dNdxyT(dNdX_dNdY4, elementSkonczony.jacobian2pkt[3], elementSkonczony); // gdzie f(x,y) to: k(t) * [dN/dx*dN/dxT + dN/dy*dN/dyT]					//

		vector<vector<double>> macierzH = sumowanie_macierzy_do_Hbc_3pkt(dNdxy1Ndxy1T, dNdxy2Ndxy2T, dNdxy3Ndxy3T, dNdxy4Ndxy4T);

		return macierzH;
	}
	else if (Dane.sc == 3)
	{
		vector<vector<double>> dNdX_dNdY1 = dNidX_dNidY(elementSkonczony.jacobianOdwrocony3pkt[0], elementUniwersalny, 1, Dane.sc); // obliczanie dla kazdego punktu calkowania pochodnych dNi/dx i dNi/dy;
		vector<vector<double>> dNdX_dNdY2 = dNidX_dNidY(elementSkonczony.jacobianOdwrocony3pkt[1], elementUniwersalny, 2, Dane.sc); // wyniki sa zapisane jako macierze 2x8, gdzie pierwszy wiersz to:
		vector<vector<double>> dNdX_dNdY3 = dNidX_dNidY(elementSkonczony.jacobianOdwrocony3pkt[2], elementUniwersalny, 3, Dane.sc); // [dN1/dx, dN2/dx, dN3/dx, dN4/dx] 
		vector<vector<double>> dNdX_dNdY4 = dNidX_dNidY(elementSkonczony.jacobianOdwrocony3pkt[3], elementUniwersalny, 4, Dane.sc); // [dN1/dy, dN2/dy, dN3/dy, dN4/dy]
		vector<vector<double>> dNdX_dNdY5 = dNidX_dNidY(elementSkonczony.jacobianOdwrocony3pkt[4], elementUniwersalny, 5, Dane.sc);
		vector<vector<double>> dNdX_dNdY6 = dNidX_dNidY(elementSkonczony.jacobianOdwrocony3pkt[5], elementUniwersalny, 6, Dane.sc);
		vector<vector<double>> dNdX_dNdY7 = dNidX_dNidY(elementSkonczony.jacobianOdwrocony3pkt[6], elementUniwersalny, 7, Dane.sc);
		vector<vector<double>> dNdX_dNdY8 = dNidX_dNidY(elementSkonczony.jacobianOdwrocony3pkt[7], elementUniwersalny, 8, Dane.sc);
		vector<vector<double>> dNdX_dNdY9 = dNidX_dNidY(elementSkonczony.jacobianOdwrocony3pkt[8], elementUniwersalny, 9, Dane.sc);

		//wypisz_macierz(dNdX_dNdY1);
		//wypisz_macierz(dNdX_dNdY2);
		//wypisz_macierz(dNdX_dNdY3);
		//wypisz_macierz(dNdX_dNdY4);
		//wypisz_macierz(dNdX_dNdY5);
		//wypisz_macierz(dNdX_dNdY6);
		//wypisz_macierz(dNdX_dNdY7);
		//wypisz_macierz(dNdX_dNdY8);
		//wypisz_macierz(dNdX_dNdY9);

		//	cout << "------------------------------------------------------------------------------------------------------------------------" << endl;
	//	cout << "dN/dX1 i dN/dY1 (kolejno dla N1, N2, N3, N4)" << endl;
	//	wypisz_macierz(dNdX1_dNdY1);
	//	cout << "dN/dX2 i dN/dY2 (kolejno dla N1, N2, N3, N4)" << endl;
	//	wypisz_macierz(dNdX2_dNdY2);
	//	cout << "dN/dX3 i dN/dY3 (kolejno dla N1, N2, N3, N4)" << endl;
	//	wypisz_macierz(dNdX3_dNdY3);
	//	cout << "dN/dX4 i dN/dY4 (kolejno dla N1, N2, N3, N4)" << endl;
	//	wypisz_macierz(dNdX4_dNdY4);
	//	cout << "------------------------------------------------------------------------------------------------------------------------" << endl;

		vector<vector<double>> dNdxyNdxy1T = dNdxy_dNdxyT(dNdX_dNdY1, elementSkonczony.jacobian3pkt[0], elementSkonczony); // obliczenie sum: dN/dx * dN/dxT + dN/dy * dN/dyT, pomnozonych przez		//podajemy jeszcze "elementSkonczony",
		vector<vector<double>> dNdxyNdxy2T = dNdxy_dNdxyT(dNdX_dNdY2, elementSkonczony.jacobian3pkt[1], elementSkonczony); // wyznaczniki jacobianow w poszczegolnych punktach oraz k(t)=30;			//który zawiera wsp. K = 30
		vector<vector<double>> dNdxyNdxy3T = dNdxy_dNdxyT(dNdX_dNdY3, elementSkonczony.jacobian3pkt[2], elementSkonczony); // wzor: [H] = f(x,y) * w1 * w2 * det[J], 								//
		vector<vector<double>> dNdxyNdxy4T = dNdxy_dNdxyT(dNdX_dNdY4, elementSkonczony.jacobian3pkt[3], elementSkonczony); // gdzie f(x,y) to: k(t) * [dN/dx*dN/dxT + dN/dy*dN/dyT]					//
		vector<vector<double>> dNdxyNdxy5T = dNdxy_dNdxyT(dNdX_dNdY5, elementSkonczony.jacobian3pkt[4], elementSkonczony);
		vector<vector<double>> dNdxyNdxy6T = dNdxy_dNdxyT(dNdX_dNdY6, elementSkonczony.jacobian3pkt[5], elementSkonczony);
		vector<vector<double>> dNdxyNdxy7T = dNdxy_dNdxyT(dNdX_dNdY7, elementSkonczony.jacobian3pkt[6], elementSkonczony);
		vector<vector<double>> dNdxyNdxy8T = dNdxy_dNdxyT(dNdX_dNdY8, elementSkonczony.jacobian3pkt[7], elementSkonczony);
		vector<vector<double>> dNdxyNdxy9T = dNdxy_dNdxyT(dNdX_dNdY9, elementSkonczony.jacobian3pkt[8], elementSkonczony);

		//wypisz_macierz(dNdxy1Ndxy1T);
		//wypisz_macierz(dNdxy2Ndxy2T);
		//wypisz_macierz(dNdxy3Ndxy3T);
		//wypisz_macierz(dNdxy4Ndxy4T);

		vector<vector<double>> macierzH = sumowanie_macierzy_3pkt_schemat(elementUniwersalny, dNdxyNdxy1T, dNdxyNdxy4T, dNdxyNdxy7T, dNdxyNdxy2T, dNdxyNdxy5T, dNdxyNdxy8T, dNdxyNdxy3T, dNdxyNdxy6T, dNdxyNdxy9T);
		return macierzH;
	}


}

vector<vector<double>> generuj_macierz_C_lokalna(Elem4 elementUniwersalny, element elementSkonczony, node* Nd, globalData Dane) //generowanie macierzy [H] lokalnej - czyli dla danego elementu
{
	if (Dane.sc == 2)
	{
		vector<vector<double>> CL1 = N_razy_Nt_dla_mac_C(elementUniwersalny, elementSkonczony.jacobian2pkt[0], 1, Dane.sc, Dane); // obliczanie dla kazdego pkt calkowania wektora 
		vector<vector<double>> CL2 = N_razy_Nt_dla_mac_C(elementUniwersalny, elementSkonczony.jacobian2pkt[1], 2, Dane.sc, Dane); // wartosci funkcji ksztaltu
		vector<vector<double>> CL3 = N_razy_Nt_dla_mac_C(elementUniwersalny, elementSkonczony.jacobian2pkt[2], 3, Dane.sc, Dane); // pomnozonych przez ich wektor transponowany
		vector<vector<double>> CL4 = N_razy_Nt_dla_mac_C(elementUniwersalny, elementSkonczony.jacobian2pkt[3], 4, Dane.sc, Dane); // oraz przemnozonych przez stale: "c" i "ro"

		vector<vector<double>> macierzC = sumowanie_macierzy_do_Hbc_3pkt(CL1, CL2, CL3, CL4);

		return macierzC;
	}
	else if (Dane.sc == 3)
	{
		vector<vector<double>> CL1 = N_razy_Nt_dla_mac_C(elementUniwersalny, elementSkonczony.jacobian3pkt[0], 1, Dane.sc, Dane); // analogicznie do schematu 2pkt: obliczanie dla kazdego pkt calkowania wektora 
		vector<vector<double>> CL2 = N_razy_Nt_dla_mac_C(elementUniwersalny, elementSkonczony.jacobian3pkt[1], 2, Dane.sc, Dane); // wartosci funkcji ksztaltu
		vector<vector<double>> CL3 = N_razy_Nt_dla_mac_C(elementUniwersalny, elementSkonczony.jacobian3pkt[2], 3, Dane.sc, Dane); // pomnozonych przez ich wektor transponowany
		vector<vector<double>> CL4 = N_razy_Nt_dla_mac_C(elementUniwersalny, elementSkonczony.jacobian3pkt[3], 4, Dane.sc, Dane); // oraz przemnozonych przez stale: "c" i "ro"
		vector<vector<double>> CL5 = N_razy_Nt_dla_mac_C(elementUniwersalny, elementSkonczony.jacobian3pkt[4], 5, Dane.sc, Dane); //
		vector<vector<double>> CL6 = N_razy_Nt_dla_mac_C(elementUniwersalny, elementSkonczony.jacobian3pkt[5], 6, Dane.sc, Dane); //
		vector<vector<double>> CL7 = N_razy_Nt_dla_mac_C(elementUniwersalny, elementSkonczony.jacobian3pkt[6], 7, Dane.sc, Dane); //
		vector<vector<double>> CL8 = N_razy_Nt_dla_mac_C(elementUniwersalny, elementSkonczony.jacobian3pkt[7], 8, Dane.sc, Dane); //
		vector<vector<double>> CL9 = N_razy_Nt_dla_mac_C(elementUniwersalny, elementSkonczony.jacobian3pkt[8], 9, Dane.sc, Dane); //

		//wypisz_macierz(CL1);
		//wypisz_macierz(CL2);
		//wypisz_macierz(CL3);
		//wypisz_macierz(CL4);
		//wypisz_macierz(CL5);
		//wypisz_macierz(CL6);
		//wypisz_macierz(CL7);
		//wypisz_macierz(CL8);
		//wypisz_macierz(CL9);

		vector<vector<double>> macierzC = sumowanie_macierzy_3pkt_schemat(elementUniwersalny, CL1, CL4, CL7, CL2, CL5, CL8, CL3, CL6, CL9);
		//wypisz_macierz(macierzC);

		return macierzC;
	}
}

vector<vector<double>> N_razy_Nt_dla_mac_Hbc(Elem4 elementUniwersalny, element elementSkonczony, node* Nd, int nrPktCalkowania, globalData dane)
{
	vector<vector<double>> result;
	vector<double> temp;
	nrPktCalkowania--;

	if (dane.sc == 2)
	{
		for (int p = 0; p < 4; p++)
			temp.push_back(elementUniwersalny.macierzN_poPowierzchni_2pkt[nrPktCalkowania][0] * elementUniwersalny.macierzN_poPowierzchni_2pkt[nrPktCalkowania][p]);
		result.push_back(temp);
		temp.clear();

		for (int p = 0; p < 4; p++)
			temp.push_back(elementUniwersalny.macierzN_poPowierzchni_2pkt[nrPktCalkowania][1] * elementUniwersalny.macierzN_poPowierzchni_2pkt[nrPktCalkowania][p]);
		result.push_back(temp);
		temp.clear();

		for (int p = 0; p < 4; p++)
			temp.push_back(elementUniwersalny.macierzN_poPowierzchni_2pkt[nrPktCalkowania][2] * elementUniwersalny.macierzN_poPowierzchni_2pkt[nrPktCalkowania][p]);
		result.push_back(temp);
		temp.clear();

		for (int p = 0; p < 4; p++)
			temp.push_back(elementUniwersalny.macierzN_poPowierzchni_2pkt[nrPktCalkowania][3] * elementUniwersalny.macierzN_poPowierzchni_2pkt[nrPktCalkowania][p]);
		result.push_back(temp);
		temp.clear();

		return result;
	}

	else if (dane.sc == 3)
	{
		for (int p = 0; p < 4; p++)
			temp.push_back(elementUniwersalny.macierzN_poPowierzchni_3pkt[nrPktCalkowania][0] * elementUniwersalny.macierzN_poPowierzchni_3pkt[nrPktCalkowania][p]);
		result.push_back(temp);
		temp.clear();

		for (int p = 0; p < 4; p++)
			temp.push_back(elementUniwersalny.macierzN_poPowierzchni_3pkt[nrPktCalkowania][1] * elementUniwersalny.macierzN_poPowierzchni_3pkt[nrPktCalkowania][p]);
		result.push_back(temp);
		temp.clear();

		for (int p = 0; p < 4; p++)
			temp.push_back(elementUniwersalny.macierzN_poPowierzchni_3pkt[nrPktCalkowania][2] * elementUniwersalny.macierzN_poPowierzchni_3pkt[nrPktCalkowania][p]);
		result.push_back(temp);
		temp.clear();

		for (int p = 0; p < 4; p++)
			temp.push_back(elementUniwersalny.macierzN_poPowierzchni_3pkt[nrPktCalkowania][3] * elementUniwersalny.macierzN_poPowierzchni_3pkt[nrPktCalkowania][p]);
		result.push_back(temp);
		temp.clear();

		return result;
	}


}

vector<vector<double>> generuj_macierz_Hbc_lokalna(Elem4 elementUniwersalny, element elementSkonczony, node* Nd, globalData dane)
{
	double detJ;

	if (dane.sc == 2)
	{
		vector<vector<double>> Hbc_bok_dolny_pkt_1 = N_razy_Nt_dla_mac_Hbc(elementUniwersalny, elementSkonczony, Nd, 1, dane);
		vector<vector<double>> Hbc_bok_dolny_pkt_2 = N_razy_Nt_dla_mac_Hbc(elementUniwersalny, elementSkonczony, Nd, 2, dane);
		vector<vector<double>> Hbc_bok_prawy_pkt_1 = N_razy_Nt_dla_mac_Hbc(elementUniwersalny, elementSkonczony, Nd, 3, dane);
		vector<vector<double>> Hbc_bok_prawy_pkt_2 = N_razy_Nt_dla_mac_Hbc(elementUniwersalny, elementSkonczony, Nd, 4, dane);
		vector<vector<double>> Hbc_bok_gorny_pkt_1 = N_razy_Nt_dla_mac_Hbc(elementUniwersalny, elementSkonczony, Nd, 5, dane);
		vector<vector<double>> Hbc_bok_gorny_pkt_2 = N_razy_Nt_dla_mac_Hbc(elementUniwersalny, elementSkonczony, Nd, 6, dane);
		vector<vector<double>> Hbc_bok_lewy_pkt_1 = N_razy_Nt_dla_mac_Hbc(elementUniwersalny, elementSkonczony, Nd, 7, dane);
		vector<vector<double>> Hbc_bok_lewy_pkt_2 = N_razy_Nt_dla_mac_Hbc(elementUniwersalny, elementSkonczony, Nd, 8, dane);

		vector<vector<double>> sumaMacierzyHbcOgrzewanychBokow;
		for (int i = 0; i < 4; i++)
		{
			vector<double> temp;
			for (int j = 0; j < 4; j++)
			{
				temp.push_back(0);
			}
			sumaMacierzyHbcOgrzewanychBokow.push_back(temp);
			temp.clear();
		}

		//zsumowanie macierzy dla pkt calkowania lezacych na ogrzewanych bokach//
		if (Nd[elementSkonczony.id[0] - 1].BC == 1 && Nd[elementSkonczony.id[1] - 1].BC == 1)
		{
			vector<vector<double>> temp = sumaMacierzyHbcOgrzewanychBokow;
			temp = sumowanie_macierzy_do_Hbc_3pkt(Hbc_bok_dolny_pkt_1, Hbc_bok_dolny_pkt_2);
			detJ = dane.w / (dane.nW - 1) / 2;
			temp = mnozenie_macierzy_przez_stala(temp, detJ);
			sumaMacierzyHbcOgrzewanychBokow = sumowanie_macierzy_do_Hbc_3pkt(sumaMacierzyHbcOgrzewanychBokow, temp);
		}
		if (Nd[elementSkonczony.id[1] - 1].BC == 1 && Nd[elementSkonczony.id[2] - 1].BC == 1)
		{
			vector<vector<double>> temp = sumaMacierzyHbcOgrzewanychBokow;
			temp = sumowanie_macierzy_do_Hbc_3pkt(Hbc_bok_prawy_pkt_1, Hbc_bok_prawy_pkt_2);
			detJ = dane.h / (dane.nH - 1) / 2;
			temp = mnozenie_macierzy_przez_stala(temp, detJ);
			sumaMacierzyHbcOgrzewanychBokow = sumowanie_macierzy_do_Hbc_3pkt(sumaMacierzyHbcOgrzewanychBokow, temp);
		}
		if (Nd[elementSkonczony.id[2] - 1].BC == 1 && Nd[elementSkonczony.id[3] - 1].BC == 1)
		{
			vector<vector<double>> temp = sumaMacierzyHbcOgrzewanychBokow;
			temp = sumowanie_macierzy_do_Hbc_3pkt(Hbc_bok_gorny_pkt_1, Hbc_bok_gorny_pkt_2);
			detJ = dane.w / (dane.nW - 1) / 2;
			temp = mnozenie_macierzy_przez_stala(temp, detJ);
			sumaMacierzyHbcOgrzewanychBokow = sumowanie_macierzy_do_Hbc_3pkt(sumaMacierzyHbcOgrzewanychBokow, temp);
		}
		if (Nd[elementSkonczony.id[3] - 1].BC == 1 && Nd[elementSkonczony.id[0] - 1].BC == 1)
		{
			vector<vector<double>> temp = sumaMacierzyHbcOgrzewanychBokow;
			temp = sumowanie_macierzy_do_Hbc_3pkt(Hbc_bok_lewy_pkt_1, Hbc_bok_lewy_pkt_2);
			detJ = dane.h / (dane.nH - 1) / 2;
			temp = mnozenie_macierzy_przez_stala(temp, detJ);
			sumaMacierzyHbcOgrzewanychBokow = sumowanie_macierzy_do_Hbc_3pkt(sumaMacierzyHbcOgrzewanychBokow, temp);
		}
		//przemnazanie przez alfe
		sumaMacierzyHbcOgrzewanychBokow = mnozenie_macierzy_przez_stala(sumaMacierzyHbcOgrzewanychBokow, dane.alfa);

		return sumaMacierzyHbcOgrzewanychBokow;
	}
	else if (dane.sc == 3)
	{
		vector<vector<double>> Hbc_bok_dolny_pkt_1 = N_razy_Nt_dla_mac_Hbc(elementUniwersalny, elementSkonczony, Nd, 1, dane);
		vector<vector<double>> Hbc_bok_dolny_pkt_2 = N_razy_Nt_dla_mac_Hbc(elementUniwersalny, elementSkonczony, Nd, 2, dane);
		vector<vector<double>> Hbc_bok_dolny_pkt_3 = N_razy_Nt_dla_mac_Hbc(elementUniwersalny, elementSkonczony, Nd, 3, dane);
		vector<vector<double>> Hbc_bok_prawy_pkt_1 = N_razy_Nt_dla_mac_Hbc(elementUniwersalny, elementSkonczony, Nd, 4, dane);
		vector<vector<double>> Hbc_bok_prawy_pkt_2 = N_razy_Nt_dla_mac_Hbc(elementUniwersalny, elementSkonczony, Nd, 5, dane);
		vector<vector<double>> Hbc_bok_prawy_pkt_3 = N_razy_Nt_dla_mac_Hbc(elementUniwersalny, elementSkonczony, Nd, 6, dane);
		vector<vector<double>> Hbc_bok_gorny_pkt_1 = N_razy_Nt_dla_mac_Hbc(elementUniwersalny, elementSkonczony, Nd, 7, dane);
		vector<vector<double>> Hbc_bok_gorny_pkt_2 = N_razy_Nt_dla_mac_Hbc(elementUniwersalny, elementSkonczony, Nd, 8, dane);
		vector<vector<double>> Hbc_bok_gorny_pkt_3 = N_razy_Nt_dla_mac_Hbc(elementUniwersalny, elementSkonczony, Nd, 9, dane);
		vector<vector<double>> Hbc_bok_lewy_pkt_1 = N_razy_Nt_dla_mac_Hbc(elementUniwersalny, elementSkonczony, Nd, 10, dane);
		vector<vector<double>> Hbc_bok_lewy_pkt_2 = N_razy_Nt_dla_mac_Hbc(elementUniwersalny, elementSkonczony, Nd, 11, dane);
		vector<vector<double>> Hbc_bok_lewy_pkt_3 = N_razy_Nt_dla_mac_Hbc(elementUniwersalny, elementSkonczony, Nd, 12, dane);

		vector<vector<double>> sumaMacierzyHbcOgrzewanychBokow;
		for (int i = 0; i < 4; i++)
		{
			vector<double> temp;
			for (int j = 0; j < 4; j++)
			{
				temp.push_back(0);
			}
			sumaMacierzyHbcOgrzewanychBokow.push_back(temp);
			temp.clear();
		}

		//zsumowanie macierzy dla pkt calkowania lezacych na ogrzewanych bokach - kolejno sprawdzamy poszczegolne boki//
		if (Nd[elementSkonczony.id[0] - 1].BC == 1 && Nd[elementSkonczony.id[1] - 1].BC == 1)
		{// dolny bok
			vector<vector<double>> temp = sumaMacierzyHbcOgrzewanychBokow;
			//sumowanie czesci macierzy Hbc z poszczegolnych pkt calkowania
			temp = sumowanie_macierzy_do_Hbc_3pkt(Hbc_bok_dolny_pkt_1, Hbc_bok_dolny_pkt_2, Hbc_bok_dolny_pkt_3, elementUniwersalny);
			detJ = dane.w / (dane.nW - 1) / 2;
			temp = mnozenie_macierzy_przez_stala(temp, detJ);
			sumaMacierzyHbcOgrzewanychBokow = sumowanie_macierzy_do_Hbc_3pkt(sumaMacierzyHbcOgrzewanychBokow, temp);
		}
		if (Nd[elementSkonczony.id[1] - 1].BC == 1 && Nd[elementSkonczony.id[2] - 1].BC == 1)
		{// prawy bok
			vector<vector<double>> temp = sumaMacierzyHbcOgrzewanychBokow;
			temp = sumowanie_macierzy_do_Hbc_3pkt(Hbc_bok_prawy_pkt_1, Hbc_bok_prawy_pkt_2, Hbc_bok_prawy_pkt_3, elementUniwersalny);
			detJ = dane.h / (dane.nH - 1) / 2;
			temp = mnozenie_macierzy_przez_stala(temp, detJ);
			sumaMacierzyHbcOgrzewanychBokow = sumowanie_macierzy_do_Hbc_3pkt(sumaMacierzyHbcOgrzewanychBokow, temp);
		}
		if (Nd[elementSkonczony.id[2] - 1].BC == 1 && Nd[elementSkonczony.id[3] - 1].BC == 1)
		{// gorny bok
			vector<vector<double>> temp = sumaMacierzyHbcOgrzewanychBokow;
			temp = sumowanie_macierzy_do_Hbc_3pkt(Hbc_bok_gorny_pkt_1, Hbc_bok_gorny_pkt_2, Hbc_bok_gorny_pkt_3, elementUniwersalny);
			detJ = dane.w / (dane.nW - 1) / 2;
			temp = mnozenie_macierzy_przez_stala(temp, detJ);
			sumaMacierzyHbcOgrzewanychBokow = sumowanie_macierzy_do_Hbc_3pkt(sumaMacierzyHbcOgrzewanychBokow, temp);
		}
		if (Nd[elementSkonczony.id[3] - 1].BC == 1 && Nd[elementSkonczony.id[0] - 1].BC == 1)
		{// lewy bok
			vector<vector<double>> temp = sumaMacierzyHbcOgrzewanychBokow;
			temp = sumowanie_macierzy_do_Hbc_3pkt(Hbc_bok_lewy_pkt_1, Hbc_bok_lewy_pkt_2, Hbc_bok_lewy_pkt_3, elementUniwersalny);
			detJ = dane.h / (dane.nH - 1) / 2;
			temp = mnozenie_macierzy_przez_stala(temp, detJ);
			sumaMacierzyHbcOgrzewanychBokow = sumowanie_macierzy_do_Hbc_3pkt(sumaMacierzyHbcOgrzewanychBokow, temp);
		}

		//przemnazanie przez alfe
		sumaMacierzyHbcOgrzewanychBokow = mnozenie_macierzy_przez_stala(sumaMacierzyHbcOgrzewanychBokow, dane.alfa);

		//wypisz_macierz(sumaMacierzyHbcOgrzewanychBokow);
		return sumaMacierzyHbcOgrzewanychBokow;
	}

}

vector<double> generuj_macierz_P_lokalna(Elem4 elementUniwersalny, element elementSkonczony, node* Nd, globalData dane)
{
	vector<double> P_local;
	for (int i = 0; i < 4; i++)
		P_local.push_back(0);
	double detJ;

	if (dane.sc == 2)
	{
		if (Nd[elementSkonczony.id[0] - 1].BC == 1 && Nd[elementSkonczony.id[1] - 1].BC == 1)
		{//dolny bok
			vector<double> temp;
			detJ = dane.w / (dane.nW - 1) / 2;
			temp.push_back(((-1) * dane.alfa * (elementUniwersalny.macierzN_poPowierzchni_2pkt[0][0] + elementUniwersalny.macierzN_poPowierzchni_2pkt[1][0])) * dane.t_alfa * detJ);
			temp.push_back(((-1) * dane.alfa * (elementUniwersalny.macierzN_poPowierzchni_2pkt[0][1] + elementUniwersalny.macierzN_poPowierzchni_2pkt[1][1])) * dane.t_alfa * detJ);
			temp.push_back(((-1) * dane.alfa * (elementUniwersalny.macierzN_poPowierzchni_2pkt[0][2] + elementUniwersalny.macierzN_poPowierzchni_2pkt[1][2])) * dane.t_alfa * detJ);
			temp.push_back(((-1) * dane.alfa * (elementUniwersalny.macierzN_poPowierzchni_2pkt[0][3] + elementUniwersalny.macierzN_poPowierzchni_2pkt[1][3])) * dane.t_alfa * detJ);
			P_local = sumowanie_macierzy_do_Pl(P_local, temp);
		}
		if (Nd[elementSkonczony.id[1] - 1].BC == 1 && Nd[elementSkonczony.id[2] - 1].BC == 1)
		{//prawy bok
			vector<double> temp;
			detJ = dane.h / (dane.nH - 1) / 2;
			temp.push_back(((-1) * dane.alfa * (elementUniwersalny.macierzN_poPowierzchni_2pkt[2][0] + elementUniwersalny.macierzN_poPowierzchni_2pkt[3][0])) * dane.t_alfa * detJ);
			temp.push_back(((-1) * dane.alfa * (elementUniwersalny.macierzN_poPowierzchni_2pkt[2][1] + elementUniwersalny.macierzN_poPowierzchni_2pkt[3][1])) * dane.t_alfa * detJ);
			temp.push_back(((-1) * dane.alfa * (elementUniwersalny.macierzN_poPowierzchni_2pkt[2][2] + elementUniwersalny.macierzN_poPowierzchni_2pkt[3][2])) * dane.t_alfa * detJ);
			temp.push_back(((-1) * dane.alfa * (elementUniwersalny.macierzN_poPowierzchni_2pkt[2][3] + elementUniwersalny.macierzN_poPowierzchni_2pkt[3][3])) * dane.t_alfa * detJ);
			P_local = sumowanie_macierzy_do_Pl(P_local, temp);
		}
		if (Nd[elementSkonczony.id[2] - 1].BC == 1 && Nd[elementSkonczony.id[3] - 1].BC == 1)
		{//gorny bok
			vector<double> temp;
			detJ = dane.w / (dane.nW - 1) / 2;
			temp.push_back(((-1) * dane.alfa * (elementUniwersalny.macierzN_poPowierzchni_2pkt[4][0] + elementUniwersalny.macierzN_poPowierzchni_2pkt[5][0])) * dane.t_alfa * detJ);
			temp.push_back(((-1) * dane.alfa * (elementUniwersalny.macierzN_poPowierzchni_2pkt[4][1] + elementUniwersalny.macierzN_poPowierzchni_2pkt[5][1])) * dane.t_alfa * detJ);
			temp.push_back(((-1) * dane.alfa * (elementUniwersalny.macierzN_poPowierzchni_2pkt[4][2] + elementUniwersalny.macierzN_poPowierzchni_2pkt[5][2])) * dane.t_alfa * detJ);
			temp.push_back(((-1) * dane.alfa * (elementUniwersalny.macierzN_poPowierzchni_2pkt[4][3] + elementUniwersalny.macierzN_poPowierzchni_2pkt[5][3])) * dane.t_alfa * detJ);
			P_local = sumowanie_macierzy_do_Pl(P_local, temp);
		}
		if (Nd[elementSkonczony.id[3] - 1].BC == 1 && Nd[elementSkonczony.id[0] - 1].BC == 1)
		{//lewy bok
			vector<double> temp;
			detJ = dane.h / (dane.nH - 1) / 2;
			temp.push_back(((-1) * dane.alfa * (elementUniwersalny.macierzN_poPowierzchni_2pkt[6][0] + elementUniwersalny.macierzN_poPowierzchni_2pkt[7][0])) * dane.t_alfa * detJ);
			temp.push_back(((-1) * dane.alfa * (elementUniwersalny.macierzN_poPowierzchni_2pkt[6][1] + elementUniwersalny.macierzN_poPowierzchni_2pkt[7][1])) * dane.t_alfa * detJ);
			temp.push_back(((-1) * dane.alfa * (elementUniwersalny.macierzN_poPowierzchni_2pkt[6][2] + elementUniwersalny.macierzN_poPowierzchni_2pkt[7][2])) * dane.t_alfa * detJ);
			temp.push_back(((-1) * dane.alfa * (elementUniwersalny.macierzN_poPowierzchni_2pkt[6][3] + elementUniwersalny.macierzN_poPowierzchni_2pkt[7][3])) * dane.t_alfa * detJ);
			P_local = sumowanie_macierzy_do_Pl(P_local, temp);
		}
	}
	else if (dane.sc == 3)
	{
		if (Nd[elementSkonczony.id[0] - 1].BC == 1 && Nd[elementSkonczony.id[1] - 1].BC == 1)
		{//dolny bok
			vector<double> temp;
			detJ = dane.w / (dane.nW - 1) / 2;
			temp.push_back(((-1) * dane.alfa * (elementUniwersalny.macierzN_poPowierzchni_3pkt[0][0] * elementUniwersalny.wagi3pkt[0] + elementUniwersalny.macierzN_poPowierzchni_3pkt[1][0] * elementUniwersalny.wagi3pkt[1] + elementUniwersalny.macierzN_poPowierzchni_3pkt[2][0] * elementUniwersalny.wagi3pkt[2])) * dane.t_alfa * detJ);
			temp.push_back(((-1) * dane.alfa * (elementUniwersalny.macierzN_poPowierzchni_3pkt[0][1] * elementUniwersalny.wagi3pkt[0] + elementUniwersalny.macierzN_poPowierzchni_3pkt[1][1] * elementUniwersalny.wagi3pkt[1] + elementUniwersalny.macierzN_poPowierzchni_3pkt[2][1] * elementUniwersalny.wagi3pkt[2])) * dane.t_alfa * detJ);
			temp.push_back(((-1) * dane.alfa * (elementUniwersalny.macierzN_poPowierzchni_3pkt[0][2] * elementUniwersalny.wagi3pkt[0] + elementUniwersalny.macierzN_poPowierzchni_3pkt[1][2] * elementUniwersalny.wagi3pkt[1] + elementUniwersalny.macierzN_poPowierzchni_3pkt[2][2] * elementUniwersalny.wagi3pkt[2])) * dane.t_alfa * detJ);
			temp.push_back(((-1) * dane.alfa * (elementUniwersalny.macierzN_poPowierzchni_3pkt[0][3] * elementUniwersalny.wagi3pkt[0] + elementUniwersalny.macierzN_poPowierzchni_3pkt[1][3] * elementUniwersalny.wagi3pkt[1] + elementUniwersalny.macierzN_poPowierzchni_3pkt[2][3] * elementUniwersalny.wagi3pkt[2])) * dane.t_alfa * detJ);
			P_local = sumowanie_macierzy_do_Pl(P_local, temp);
		}
		if (Nd[elementSkonczony.id[1] - 1].BC == 1 && Nd[elementSkonczony.id[2] - 1].BC == 1)
		{//prawy bok
			vector<double> temp;
			detJ = dane.h / (dane.nH - 1) / 2;
			temp.push_back(((-1) * dane.alfa * (elementUniwersalny.macierzN_poPowierzchni_3pkt[3][0] * elementUniwersalny.wagi3pkt[0] + elementUniwersalny.macierzN_poPowierzchni_3pkt[4][0] * elementUniwersalny.wagi3pkt[1] + elementUniwersalny.macierzN_poPowierzchni_3pkt[5][0] * elementUniwersalny.wagi3pkt[2])) * dane.t_alfa * detJ);
			temp.push_back(((-1) * dane.alfa * (elementUniwersalny.macierzN_poPowierzchni_3pkt[3][1] * elementUniwersalny.wagi3pkt[0] + elementUniwersalny.macierzN_poPowierzchni_3pkt[4][1] * elementUniwersalny.wagi3pkt[1] + elementUniwersalny.macierzN_poPowierzchni_3pkt[5][1] * elementUniwersalny.wagi3pkt[2])) * dane.t_alfa * detJ);
			temp.push_back(((-1) * dane.alfa * (elementUniwersalny.macierzN_poPowierzchni_3pkt[3][2] * elementUniwersalny.wagi3pkt[0] + elementUniwersalny.macierzN_poPowierzchni_3pkt[4][2] * elementUniwersalny.wagi3pkt[1] + elementUniwersalny.macierzN_poPowierzchni_3pkt[5][2] * elementUniwersalny.wagi3pkt[2])) * dane.t_alfa * detJ);
			temp.push_back(((-1) * dane.alfa * (elementUniwersalny.macierzN_poPowierzchni_3pkt[3][3] * elementUniwersalny.wagi3pkt[0] + elementUniwersalny.macierzN_poPowierzchni_3pkt[4][3] * elementUniwersalny.wagi3pkt[1] + elementUniwersalny.macierzN_poPowierzchni_3pkt[5][3] * elementUniwersalny.wagi3pkt[2])) * dane.t_alfa * detJ);
			P_local = sumowanie_macierzy_do_Pl(P_local, temp);
		}
		if (Nd[elementSkonczony.id[2] - 1].BC == 1 && Nd[elementSkonczony.id[3] - 1].BC == 1)
		{//gorny bok																																										 													      
			vector<double> temp;
			detJ = dane.w / (dane.nW - 1) / 2;
			temp.push_back(((-1) * dane.alfa * (elementUniwersalny.macierzN_poPowierzchni_3pkt[6][0] * elementUniwersalny.wagi3pkt[0] + elementUniwersalny.macierzN_poPowierzchni_3pkt[7][0] * elementUniwersalny.wagi3pkt[1] + elementUniwersalny.macierzN_poPowierzchni_3pkt[8][0] * elementUniwersalny.wagi3pkt[2])) * dane.t_alfa * detJ);
			temp.push_back(((-1) * dane.alfa * (elementUniwersalny.macierzN_poPowierzchni_3pkt[6][1] * elementUniwersalny.wagi3pkt[0] + elementUniwersalny.macierzN_poPowierzchni_3pkt[7][1] * elementUniwersalny.wagi3pkt[1] + elementUniwersalny.macierzN_poPowierzchni_3pkt[8][1] * elementUniwersalny.wagi3pkt[2])) * dane.t_alfa * detJ);
			temp.push_back(((-1) * dane.alfa * (elementUniwersalny.macierzN_poPowierzchni_3pkt[6][2] * elementUniwersalny.wagi3pkt[0] + elementUniwersalny.macierzN_poPowierzchni_3pkt[7][2] * elementUniwersalny.wagi3pkt[1] + elementUniwersalny.macierzN_poPowierzchni_3pkt[8][2] * elementUniwersalny.wagi3pkt[2])) * dane.t_alfa * detJ);
			temp.push_back(((-1) * dane.alfa * (elementUniwersalny.macierzN_poPowierzchni_3pkt[6][3] * elementUniwersalny.wagi3pkt[0] + elementUniwersalny.macierzN_poPowierzchni_3pkt[7][3] * elementUniwersalny.wagi3pkt[1] + elementUniwersalny.macierzN_poPowierzchni_3pkt[8][3] * elementUniwersalny.wagi3pkt[2])) * dane.t_alfa * detJ);
			P_local = sumowanie_macierzy_do_Pl(P_local, temp);
		}
		if (Nd[elementSkonczony.id[3] - 1].BC == 1 && Nd[elementSkonczony.id[0] - 1].BC == 1)
		{//lewy bok
			vector<double> temp;
			detJ = dane.h / (dane.nH - 1) / 2;
			temp.push_back(((-1) * dane.alfa * (elementUniwersalny.macierzN_poPowierzchni_3pkt[9][0] * elementUniwersalny.wagi3pkt[0] + elementUniwersalny.macierzN_poPowierzchni_3pkt[10][0] * elementUniwersalny.wagi3pkt[1] + elementUniwersalny.macierzN_poPowierzchni_3pkt[11][0] * elementUniwersalny.wagi3pkt[2])) * dane.t_alfa * detJ);
			temp.push_back(((-1) * dane.alfa * (elementUniwersalny.macierzN_poPowierzchni_3pkt[9][1] * elementUniwersalny.wagi3pkt[0] + elementUniwersalny.macierzN_poPowierzchni_3pkt[10][1] * elementUniwersalny.wagi3pkt[1] + elementUniwersalny.macierzN_poPowierzchni_3pkt[11][1] * elementUniwersalny.wagi3pkt[2])) * dane.t_alfa * detJ);
			temp.push_back(((-1) * dane.alfa * (elementUniwersalny.macierzN_poPowierzchni_3pkt[9][2] * elementUniwersalny.wagi3pkt[0] + elementUniwersalny.macierzN_poPowierzchni_3pkt[10][2] * elementUniwersalny.wagi3pkt[1] + elementUniwersalny.macierzN_poPowierzchni_3pkt[11][2] * elementUniwersalny.wagi3pkt[2])) * dane.t_alfa * detJ);
			temp.push_back(((-1) * dane.alfa * (elementUniwersalny.macierzN_poPowierzchni_3pkt[9][3] * elementUniwersalny.wagi3pkt[0] + elementUniwersalny.macierzN_poPowierzchni_3pkt[10][3] * elementUniwersalny.wagi3pkt[1] + elementUniwersalny.macierzN_poPowierzchni_3pkt[11][3] * elementUniwersalny.wagi3pkt[2])) * dane.t_alfa * detJ);
			P_local = sumowanie_macierzy_do_Pl(P_local, temp);
		}
	}


	return P_local;
}

vector<vector<double>> macierzH_GLOBALNA(element* E, node* Nd, globalData Dane)
{
	vector<vector<double>> HG;
	vector<double> temp;
	for (int i = 0; i < Dane.nN; i++)
	{
		for (int j = 0; j < Dane.nN; j++)
		{
			temp.push_back(0);
		}
		HG.push_back(temp);
		temp.clear();
	}

	for (int k = 0; k < Dane.nE; k++)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				HG[E[k].id[i] - 1][E[k].id[j] - 1] += E[k].Hl[i][j];
			}
		}
	}
	return HG;
}

vector<vector<double>> macierzC_GLOBALNA(element* E, node* Nd, globalData Dane)
{
	vector<vector<double>> CG;
	vector<double> temp;
	for (int i = 0; i < Dane.nN; i++)
	{
		for (int j = 0; j < Dane.nN; j++)
		{
			temp.push_back(0);
		}
		CG.push_back(temp);
		temp.clear();
	}

	for (int k = 0; k < Dane.nE; k++)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				CG[E[k].id[i] - 1][E[k].id[j] - 1] += E[k].Cl[i][j];
			}
		}
	}
	return CG;
}

vector<double> macierzP_GLOBALNA(element* E, node* Nd, globalData Dane)
{
	vector<double> P;
	for (int j = 0; j < Dane.nN; j++)
		P.push_back(0);

	for (int k = 0; k < Dane.nE; k++)
	{
		for (int i = 0; i < 4; i++)
		{
			P[E[k].id[i] - 1] += E[k].Pl[i];
		}
	}
	return P;
}

vector<double> mnozenie_C_tau_t_0_plus_P(SoE soe, globalData Dane, node* Nd)
{
	vector<double> result;
	vector<double> t_0; //musimy sobie zczytac temepratury z wezlow
	vector<vector<double>> C_przez_tau = mnozenie_macierzy_przez_stala(soe.CG, 1 / Dane.dT);

	for (int i = 0; i < Dane.nN; i++)
		t_0.push_back(Nd[i].t_0);

	for (int i = 0; i < C_przez_tau.size(); i++)
	{
		double temp = 0;
		for (int j = 0; j < C_przez_tau[i].size(); j++)
		{
			temp += (-1) * C_przez_tau[i][j] * t_0[j];
		}
		result.push_back(soe.PG[i] + temp);
	}

	for (int i = 0; i < Dane.nN; i++)
		result[i] *= (-1); //od razu przemnazam przez -1, czyli tak jakby przenosze na druga strone rownania

	return result;
}

//funkcje do odwracania macierzy//
void calculateInverse(vector<vector<double>>& A) {
	int n = A.size();

	for (int i = 0; i < n; i++) {
		// Search for maximum in this column
		double maxEl = abs(A[i][i]);
		int maxRow = i;
		for (int k = i + 1; k < n; k++) {
			if (abs(A[k][i]) > maxEl) {
				maxEl = A[k][i];
				maxRow = k;
			}
		}

		// Swap maximum row with current row (column by column)
		for (int k = i; k < 2 * n; k++) {
			double tmp = A[maxRow][k];
			A[maxRow][k] = A[i][k];
			A[i][k] = tmp;
		}

		// Make all rows below this one 0 in current column
		for (int k = i + 1; k < n; k++) {
			double c = -A[k][i] / A[i][i];
			for (int j = i; j < 2 * n; j++) {
				if (i == j) {
					A[k][j] = 0;
				}
				else {
					A[k][j] += c * A[i][j];
				}
			}
		}
	}

	// Solve equation Ax=b for an upper triangular matrix A
	for (int i = n - 1; i >= 0; i--) {
		for (int k = n; k < 2 * n; k++) {
			A[i][k] /= A[i][i];
		}
		// this is not necessary, but the output looks nicer:
		A[i][i] = 1;

		for (int rowModify = i - 1; rowModify >= 0; rowModify--) {
			for (int columModify = n; columModify < 2 * n; columModify++) {
				A[rowModify][columModify] -= A[i][columModify]
					* A[rowModify][i];
			}
			// this is not necessary, but the output looks nicer:
			A[rowModify][i] = 0;
		}
	}
}

vector<vector<double>> odwrocMacierz(SoE soe, globalData dane, vector<vector<double>> macierz)
{
	int n = dane.nN;

	vector<double> line(2 * n, 0);
	vector<vector<double>> A(n, line);

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			A[i][j] = macierz[i][j];

	for (int i = 0; i < n; i++)
		A[i][n + i] = 1;

	calculateInverse(A);

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			macierz[i][j] = A[i][j + n];

	return macierz;
}

void obliczenie_t1(SoE soe, globalData dane, node* Nd)
{
	vector<double> result;
	for (int i = 0; i < dane.nN; i++)
		result.push_back(0);


	for (int i = 0; i < dane.nN; i++)
	{
		for (int j = 0; j < dane.nN; j++)
		{
			result[i] += soe.HG[i][j] * soe.PG[j];
		}
		Nd[i].t_1 = result[i];
		Nd[i].t_0 = Nd[i].t_1;
	}

	//wypisz_wektor(result);
}

int main()
{
	Elem4 elementUniwersalny; //jego konstruktor z automatycznie oblicza pochodne
	SoE soe;


// TUTAJ WPISUJEMY DANE WEJŚCIOWE - co jest czym - 58 linijka //
	globalData Dane = { 100, 500, 50, 1200, 300, 0.1, 0.1, 4.0, 4.0, 700, 25, 7800, 2 };

	element* Elem = new element[Dane.nE]; //tablica elementow skonczonch (kazdy element ma 4 ID, bo przypadają na niego 4 wezły)
	node* Nd = new node[Dane.nN]; //tablica wezlow - indeks w tej tabeli to ID danego wezlu
	for (int i = 0; i < Dane.nE; i++) { Elem[i].K = Dane.K; }
	for (int i = 0; i < Dane.nN; i++) { Nd[i].t_0 = Dane.t_0; }

	setElementsAndNodes(Elem, Nd, Dane); //generowanie siatki elementow skonczonych (wspolrzednych wezlow oraz ich ID, które są przypisane do elementów skonczonych)

	for (double t = 1; t <= (Dane.sim_time / Dane.dT); t += 1.0)
	{
		generuj_jacobiany_i_jac_odwrocone(elementUniwersalny, Elem, Nd, Dane);

		for (int i = 0; i < Dane.nE; i++)
		{// generowanie macierzy H, C, Hbc lokalnych dla każdego elementu skończonego z siatki
			Elem[i].Hl = generuj_macierz_H_lokalna(elementUniwersalny, Elem[i], Nd, Dane);
			Elem[i].Cl = generuj_macierz_C_lokalna(elementUniwersalny, Elem[i], Nd, Dane);
			Elem[i].Hbc = generuj_macierz_Hbc_lokalna(elementUniwersalny, Elem[i], Nd, Dane);

			Elem[i].Hl = suma_macierzy(Elem[i].Hl, Elem[i].Hbc); //dodawanie lokalnych macierzy Hbc z war. brzegowym do poczatkowych macierzy lokalnych H
			//cout << endl << " Macierz [H] dla elementu nr " << i + 1 << endl;
			//wypisz_macierz(Elem[i].Hl);

		// generowanie lokalnych macierzy P
			Elem[i].Pl = generuj_macierz_P_lokalna(elementUniwersalny, Elem[i], Nd, Dane);
			//cout << endl << " Macierz [P] dla elementu nr " << i + 1 << endl;
			//wypisz_wektor(Elem[i].Pl);
		}

		/***obliczenia macierzy H, C, Hbc globalnych i wstawienie ich do obiektu struktury "soe"***/
		soe.HG = macierzH_GLOBALNA(Elem, Nd, Dane);
		soe.CG = macierzC_GLOBALNA(Elem, Nd, Dane);
		soe.PG = macierzP_GLOBALNA(Elem, Nd, Dane);
		//cout << "Macierz H globalna (z uwzglednieniem warunku brzegowego)" << endl;
		//wypisz_macierz(soe.HG);
		//cout << endl << "Macierz C globalna" << endl;
		//wypisz_macierz(soe.CG);
		//cout << endl << "Macierz P globalna" << endl;
		//wypisz_wektor(soe.PG);


		/*obliczanie zmian w czasie*/

		soe.HG = dodaj_macierze(soe.HG, mnozenie_macierzy_przez_stala(soe.CG, 1 / Dane.dT)); //obliczenie [H] + [C] / d_tau
		soe.PG = mnozenie_C_tau_t_0_plus_P(soe, Dane, Nd); //obliczenie ([C]/d_tau) + [P] i przeniesienie na prawa str rownania

		//cout << "pierwszy skladnik" << endl;
		//wypisz_macierz(soe.HG);
		//cout << "drugi skladnik" << endl;
		//wypisz_wektor(soe.PG);

		/*******/

		soe.HG = odwrocMacierz(soe, Dane, soe.HG); //odwrocenie pierwszego skladnika H^
		obliczenie_t1(soe, Dane, Nd); //mnozenie_Hodw_razy_Pz_daszkiem; tutaj tak naprawde mamy juz nasze wyniki temperatury

		double max_temp, min_temp;
		min_temp = Nd[0].t_1;
		max_temp = Nd[0].t_1;
		for (int i = 1; i < Dane.nN; i++)
		{
			if (Nd[i].t_1 < min_temp) min_temp = Nd[i].t_1;
			if (Nd[i].t_1 > max_temp) max_temp = Nd[i].t_1;
		}

		cout << "Time[s] = " << t * Dane.dT << "\tMinTemp[s] = " << min_temp << "\tMaxTemp[s] = " << max_temp << endl;
	}

	return 0;
}
