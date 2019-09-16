#ifndef LEGENDREPOLYNOMIAL_H
#define LEGENDREPOLYNOMIAL_H

//-------------------------Legendre����ʽ-----------------------------------
double getLegendrePolynomial(double x, int p)
{
	double x2 = x * x;
	double value;
	switch (p)
	{
	case 0:		value = 1.0;
		break;
	case 1:		value = x;
		break;
	case 2:		value = 3.0 / 2.0*x2 - 1.0 / 2.0;
		break;
	case 3:		value = (-3.0 / 2.0 + 5.0 / 2.0*x2)*x;
		break;
	case 4:		value = 3.0 / 8.0 + (-15.0 / 4.0 + 35.0 / 8.0*x2)*x2;
		break;
	case 5:		value = (15.0 / 8.0 + (-35.0 / 4.0 + 63.0 / 8.0*x2)*x2)*x;
		break;
	case 6:		value = -5.0 / 16.0 + (105.0 / 16.0 + (-315.0 / 16.0 + 231.0 / 16.0*x2)*x2)*x2;
		break;
	case 7:		value = (-35.0 / 16.0 + (315.0 / 16.0 + (-693.0 / 16.0 + 429.0 / 16.0*x2)*x2)*x2)*x;
		break;
	case 8:		value = 35.0 / 128.0 + (-315.0 / 32.0 + (3465.0 / 64.0 + (-3003.0 / 32.0 + 6435.0 / 128.0*x2)*x2)*x2)*x2;
		break;
	case 9:		value = (315.0 / 128.0 + (-1155.0 / 32.0 + (9009.0 / 64.0 + (-6435.0 / 32.0 + 12155.0 / 128.0* x2)*x2)*x2)*x2)*x;
		break;
	case 10:	value = -63.0 / 256.0 + (3465.0 / 256.0 + (-15015.0 / 128.0 + (45045.0 / 128.0 + (-109395.0 /
		256.0 + 46189.0 / 256.0*x2)*x2)*x2)*x2)*x2;
		break;
	case 11:	value = (-693.0 / 256.0 + (15015.0 / 256.0 + (-45045.0 / 128.0 + (109395.0 / 128.0 +
		(-230945.0 / 256.0 + 88179.0 / 256.0*x2)*x2)*x2)*x2)*x2)*x;
		break;
	case 12:	value = 231.0 / 1024.0 + (-9009.0 / 512.0 + (225225.0 / 1024.0 + (-255255.0 / 256.0 +
		(2078505.0 / 1024.0 + (-969969.0 / 512.0 + 676039.0 / 1024.0*x2)*x2)*x2)*x2)*x2)*x2;
		break;
	case 13:	value = (3003.0 / 1024.0 + (-45045.0 / 512.0 + (765765.0 / 1024.0 + (-692835.0 / 256.0 +
		(4849845.0 / 1024.0 + (-2028117.0 / 512.0 + 1300075.0 / 1024.0*x2)*x2)*x2)*x2)*x2)*x2)*x;
		break;
	case 14:	value = -429.0 / 2048.0 + (45045.0 / 2048.0 + (-765765.0 / 2048.0 + (4849845.0 / 2048.0 +
		(-14549535.0 / 2048.0 + (22309287.0 / 2048.0 + (-16900975.0 / 2048.0 + 5014575.0 / 2048.0*x2)*x2)*x2)*x2)*x2)*x2)*x2;
		break;
	case 15:	value = (-6435.0 / 2048.0 + (255255.0 / 2048.0 + (-2909907.0 / 2048.0 +
		(14549535.0 / 2048.0 + (-37182145.0 / 2048.0 + (50702925.0 / 2048.0 + (-35102025.0 / 2048.0 + 9694845.0 /
			2048.0*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x;
		break;
	case 16:	value = 6435.0 / 32768.0 + (-109395.0 / 4096.0 + (4849845.0 / 8192.0 + (-20369349.0 /
		4096.0 + (334639305.0 / 16384.0 + (-185910725.0 / 4096.0 + (456326325.0 / 8192.0 + (-145422675.0 / 4096.0 +
			300540195.0 / 32768.0*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2;
		break;
	case 17:	value = (109395.0 / 32768.0 + (-692835.0 / 4096.0 + (20369349.0 / 8192.0 + (-66927861.0 /
		4096.0 + (929553625.0 / 16384.0 + (-456326325.0 / 4096.0 + (1017958725.0 / 8192.0 + (-300540195.0 / 4096.0 +
			583401555.0 / 32768.0*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x;
		break;
	case 18:	value = -12155.0 / 65536.0 + (2078505.0 / 65536.0 + (-14549535.0 / 16384.0 + (156165009.0
		/ 16384.0 + (-1673196525.0 / 32768.0 + (5019589575.0 / 32768.0 + (-4411154475.0 / 16384.0 +
		(4508102925.0 / 16384.0 + (-9917826435.0 / 65536.0 + 2268783825.0 / 65536.0*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2;
		break;
	case 19:	value = (-230945.0 / 65536.0 + (14549535.0 / 65536.0 + (-66927861.0 / 16384.0 + (557732175.0 /
		16384.0 + (-5019589575.0 / 32768.0 + (13233463425.0 / 32768.0 + (-10518906825.0 / 16384.0 +
		(9917826435.0 / 16384.0 + (-20419054425.0 / 65536.0 + 4418157975.0 / 65536.0*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x;
		break;
	case 20:	value = 46189.0 / 262144.0 + (-4849845.0 / 131072.0 + (334639305.0 / 262144.0 + (-557732175.0 /
		32768.0 + (15058768725.0 / 131072.0 + (-29113619535.0 / 65536.0 + (136745788725.0 / 131072.0 +
		(-49589132175.0 / 32768.0 + (347123925225.0 / 262144.0 + (-83945001525.0 / 131072.0 + 34461632205.0 /
			262144.0*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2;
		break;
	}
	return value;
}

//-------------------------Legendre����ʽ����-----------------------------------
double getLegendrePolynomialDerivative(double x, int p)
{
	double x2 = x * x;
	double value;
	switch (p)
	{
	case 0:		value = 0.0;
		break;
	case 1:		value = 1.0;
		break;
	case 2:		value = 3.0*x;
		break;
	case 3:		value = 15.0 / 2.0*x2 - 3.0 / 2.0;
		break;
	case 4:		value = (-15.0 / 2.0 + 35.0 / 2.0*x2)*x;
		break;
	case 5:		value = 15.0 / 8.0 + (-105.0 / 4.0 + 315.0 / 8.0*x2)*x2;
		break;
	case 6:		value = (105.0 / 8.0 + (-315.0 / 4.0 + 693.0 / 8.0*x2)*x2)*x;
		break;
	case 7:		value = -35.0 / 16.0 + (945.0 / 16.0 + (-3465.0 / 16.0 + 3003.0 / 16.0*x2)*x2)*x2;
		break;
	case 8:		value = (-315.0 / 16.0 + (3465.0 / 16.0 + (-9009.0 / 16.0 + 6435.0 / 16.0*x2)*x2)*x2)*x;
		break;
	case 9:		value = 315.0 / 128.0 + (-3465.0 / 32.0 + (45045.0 / 64.0 + (-45045.0 / 32.0 + 109395.0 / 128.0*x2)*x2)*x2)*x2;
		break;
	case 10:	value = (3465.0 / 128.0 + (-15015.0 / 32.0 + (135135.0 / 64.0 + (-109395.0 / 32.0 + 230945.0 /
		128.0*x2)*x2)*x2)*x2)*x;
		break;
	case 11:	value = -693.0 / 256.0 + (45045.0 / 256.0 + (-225225.0 / 128.0 + (765765.0 / 128.0 + (-2078505.0 /
		256.0 + 969969.0 / 256.0*x2)*x2)*x2)*x2)*x2;
		break;
	case 12:	value = (-9009.0 / 256.0 + (225225.0 / 256.0 + (-765765.0 / 128.0 + (2078505.0 / 128.0 + (-4849845.0 /
		256.0 + 2028117.0 / 256.0*x2)*x2)*x2)*x2)*x2)*x;
		break;
	case 13:	value = 3003.0 / 1024.0 + (-135135.0 / 512.0 + (3828825.0 / 1024.0 + (-4849845.0 / 256.0 + (43648605.0 /
		1024.0 + (-22309287.0 / 512.0 + 16900975.0 / 1024.0*x2)*x2)*x2)*x2)*x2)*x2;
		break;
	case 14:	value = (45045.0 / 1024.0 + (-765765.0 / 512.0 + (14549535.0 / 1024.0 + (-14549535.0 / 256.0 + (111546435.0 /
		1024.0 + (-50702925.0 / 512.0 + 35102025.0 / 1024.0*x2)*x2)*x2)* x2)*x2)*x2)*x;
		break;
	case 15:	value = -6435.0 / 2048.0 + (765765.0 / 2048.0 + (-14549535.0 / 2048.0 + (101846745.0 / 2048.0 +
		(-334639305.0 / 2048.0 + (557732175.0 / 2048.0 + (-456326325.0 / 2048.0 + 145422675.0 /
			2048.0*x2)*x2)*x2)*x2)*x2)*x2)*x2;
		break;
	case 16:	value = (-109395.0 / 2048.0 + (4849845.0 / 2048.0 + (-61108047.0 / 2048.0 + (334639305.0 / 2048.0 +
		(-929553625.0 / 2048.0 + (1368978975.0 / 2048.0 + (-1017958725.0 / 2048.0 +
			300540195.0 / 2048.0*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x;
		break;
	case 17:	value = 109395.0 / 32768.0 + (-2078505.0 / 4096.0 + (101846745.0 / 8192.0 + (-468495027.0 / 4096.0 +
		(8365982625.0 / 16384.0 + (-5019589575.0 / 4096.0 + (13233463425.0 / 8192.0 +
		(-4508102925.0 / 4096.0 + 9917826435.0 / 32768.0*x2)*x2)*x2)*x2)*x2)*x2) *x2)*x2;
		break;
	case 18:	value = (2078505.0 / 32768.0 + (-14549535.0 / 4096.0 + (468495027.0 / 8192.0 + (-1673196525.0 / 4096.0 +
		(25097947875.0 / 16384.0 + (-13233463425.0 / 4096.0 + (31556720475.0 / 8192.0 +
		(-9917826435.0 / 4096.0 + 20419054425.0 / 32768.0*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x;
		break;
	case 19:	value = -230945.0 / 65536.0 + (43648605.0 / 65536.0 + (-334639305.0 / 16384.0 + (3904125225.0 / 16384.0 +
		(-45176306175.0 / 32768.0 + (145568097675.0 / 32768.0 + (-136745788725.0 / 16384.0 +
		(148767396525.0 / 16384.0 + (-347123925225.0 / 65536.0 + 83945001525.0 / 65536.0*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2;
		break;
	case 20:	value = (-4849845.0 / 65536.0 + (334639305.0 / 65536.0 + (-1673196525.0 / 16384.0 + (15058768725.0 /
		16384.0 + (-145568097675.0 / 32768.0 + (410237366175.0 / 32768.0 + (-347123925225.0 / 16384.0 +
		(347123925225.0 / 16384.0 + (-755505013725.0 / 65536.0 + 172308161025.0 / 65536.0*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x;
		break;
	}
	return value;
}

#endif