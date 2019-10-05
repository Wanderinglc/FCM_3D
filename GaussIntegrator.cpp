
#include "GaussIntegrator.h"


void GaussIntegrator::formGaussCoordinates(double* GaussPoints)
{
	
	switch (NIP)
	{
	case  1:
		GaussPoints[0] = 0.00000000000000e-01;
		break;
	case  2:
		GaussPoints[0] = -5.773502691896258e-01;
		GaussPoints[0] = 5.773502691896258e-01;
		break;
	case  3:
		GaussPoints[0] = -7.745966692414834e-01;
		GaussPoints[1] = 0.000000000000000e-01;
		GaussPoints[2] = 7.745966692414834e-01;
		break;
	case  4:
		GaussPoints[0] = -8.611363115940526e-01;
		GaussPoints[1] = -3.399810435848563e-01;
		GaussPoints[2] = 3.399810435848563e-01;
		GaussPoints[3] = 8.611363115940526e-01;
		break;
	case  5:
		GaussPoints[0] = -9.061798459386640e-01;
		GaussPoints[1] = -5.384693101056831e-01;
		GaussPoints[2] = 0.000000000000000e-01;
		GaussPoints[3] = 5.384693101056831e-01;
		GaussPoints[4] = 9.061798459386640e-01;
		break;
	case  6:
		GaussPoints[0] = -9.324695142031520e-01;
		GaussPoints[1] = -6.612093864662645e-01;
		GaussPoints[2] = -2.386191860831969e-01;
		GaussPoints[3] = 2.386191860831969e-01;
		GaussPoints[4] = 6.612093864662645e-01;
		GaussPoints[5] = 9.324695142031520e-01;
		break;
	case  7:
		GaussPoints[0] = -9.491079123427585e-01;
		GaussPoints[1] = -7.415311855993944e-01;
		GaussPoints[2] = -4.058451513773972e-01;
		GaussPoints[3] = 0.000000000000000e-01;
		GaussPoints[4] = 4.058451513773972e-01;
		GaussPoints[5] = 7.415311855993944e-01;
		GaussPoints[6] = 9.491079123427585e-01;
		break;
	case  8:
		GaussPoints[0] = -9.602898564975362e-01;
		GaussPoints[1] = -7.966664774136267e-01;
		GaussPoints[2] = -5.255324099163290e-01;
		GaussPoints[3] = -1.834346424956498e-01;
		GaussPoints[4] = 1.834346424956498e-01;
		GaussPoints[5] = 5.255324099163290e-01;
		GaussPoints[6] = 7.966664774136267e-01;
		GaussPoints[7] = 9.602898564975362e-01;
		break;
	case  9:
		GaussPoints[0] = -9.681602395076261e-01;
		GaussPoints[1] = -8.360311073266358e-01;
		GaussPoints[2] = -6.133714327005904e-01;
		GaussPoints[3] = -3.242534234038089e-01;
		GaussPoints[4] = 0.000000000000000e-01;
		GaussPoints[5] = 3.242534234038089e-01;
		GaussPoints[6] = 6.133714327005904e-01;
		GaussPoints[7] = 8.360311073266358e-01;
		GaussPoints[8] = 9.681602395076261e-01;
		break;
	case  10:
		GaussPoints[0] = -9.739065285171717e-01;
		GaussPoints[1] = -8.650633666889845e-01;
		GaussPoints[2] = -6.794095682990244e-01;
		GaussPoints[3] = -4.333953941292472e-01;
		GaussPoints[4] = -1.488743389816312e-01;
		GaussPoints[5] = 1.488743389816312e-01;
		GaussPoints[6] = 4.333953941292472e-01;
		GaussPoints[7] = 6.794095682990244e-01;
		GaussPoints[8] = 8.650633666889845e-01;
		GaussPoints[9] = 9.739065285171717e-01;
		break;
	case 11:
		GaussPoints[0] = -9.782286581460570e-01;
		GaussPoints[1] = -8.870625997680953e-01;
		GaussPoints[2] = -7.301520055740493e-01;
		GaussPoints[3] = -5.190961292068118e-01;
		GaussPoints[4] = -2.695431559523450e-01;
		GaussPoints[5] = 0.000000000000000e-01;
		GaussPoints[6] = 2.695431559523450e-01;
		GaussPoints[7] = 5.190961292068118e-01;
		GaussPoints[8] = 7.301520055740493e-01;
		GaussPoints[9] = 8.870625997680953e-01;
		GaussPoints[10] = 9.782286581460570e-01;
		break;
	case 12:
		GaussPoints[0] = -9.815606342467193e-01;
		GaussPoints[1] = -9.041172563704749e-01;
		GaussPoints[2] = -7.699026741943047e-01;
		GaussPoints[3] = -5.873179542866174e-01;
		GaussPoints[4] = -3.678314989981802e-01;
		GaussPoints[5] = -1.252334085114689e-01;
		GaussPoints[6] = 1.252334085114689e-01;
		GaussPoints[7] = 3.678314989981802e-01;
		GaussPoints[8] = 5.873179542866174e-01;
		GaussPoints[9] = 7.699026741943047e-01;
		GaussPoints[10] = 9.041172563704749e-01;
		GaussPoints[11] = 9.815606342467193e-01;
		break;

	case 14:
		GaussPoints[0] = -9.862838086968123e-01;
		GaussPoints[1] = -9.284348836635735e-01;
		GaussPoints[2] = -8.272013150697650e-01;
		GaussPoints[3] = -6.872929048116855e-01;
		GaussPoints[4] = -5.152486363581541e-01;
		GaussPoints[5] = -3.191123689278898e-01;
		GaussPoints[6] = -1.080549487073437e-01;
		GaussPoints[7] = 1.080549487073437e-01;
		GaussPoints[8] = 3.191123689278898e-01;
		GaussPoints[9] = 5.152486363581541e-01;
		GaussPoints[10] = 6.872929048116855e-01;
		GaussPoints[11] = 8.272013150697650e-01;
		GaussPoints[12] = 9.284348836635735e-01;
		GaussPoints[13] = 9.862838086968123e-01;
		break;
	case 16:
		GaussPoints[0] = -9.894009349916499e-01;
		GaussPoints[1] = -9.445750230732326e-01;
		GaussPoints[2] = -8.656312023878317e-01;
		GaussPoints[3] = -7.554044083550030e-01;
		GaussPoints[4] = -6.178762444026437e-01;
		GaussPoints[5] = -4.580167776572274e-01;
		GaussPoints[6] = -2.816035507792589e-01;
		GaussPoints[7] = -9.501250983763744e-02;
		GaussPoints[8] = 9.501250983763744e-02;
		GaussPoints[9] = 2.816035507792589e-01;
		GaussPoints[10] = 4.580167776572274e-01;
		GaussPoints[11] = 6.178762444026437e-01;
		GaussPoints[12] = 7.554044083550030e-01;
		GaussPoints[13] = 8.656312023878317e-01;
		GaussPoints[14] = 9.445750230732326e-01;
		GaussPoints[15] = 9.894009349916499e-01;
		break;
	}
	//return GaussPoints;
	
}


void GaussIntegrator::formGaussWeights(double* GaussWeights)
{
	switch (NIP)
	{
	case  1:
		GaussWeights[0] = 2.000000000000000e+00;
		break;
	case  2:
		GaussWeights[0] = 1.000000000000000e+00;
		GaussWeights[1] = 1.000000000000000e+00;
		break;
	case  3:
		GaussWeights[0] = 5.555555555555556e-01;
		GaussWeights[1] = 8.888888888888889e-01;
		GaussWeights[2] = 5.555555555555556e-01;
		break;
	case  4:
		GaussWeights[0] = 3.478548451374539e-01;
		GaussWeights[1] = 6.521451548625461e-01;
		GaussWeights[2] = 6.521451548625461e-01;
		GaussWeights[3] = 3.478548451374539e-01;
		break;
	case  5:
		GaussWeights[0] = 2.369268850561891e-01;
		GaussWeights[1] = 4.786286704993665e-01;
		GaussWeights[2] = 5.688888888888889e-01;
		GaussWeights[3] = 4.786286704993665e-01;
		GaussWeights[4] = 2.369268850561891e-01;
		break;
	case  6:
		GaussWeights[0] = 1.713244923791703e-01;
		GaussWeights[1] = 3.607615730481386e-01;
		GaussWeights[2] = 4.679139345726910e-01;
		GaussWeights[3] = 4.679139345726910e-01;
		GaussWeights[4] = 3.607615730481386e-01;
		GaussWeights[5] = 1.713244923791703e-01;
		break;
	case  7:
		GaussWeights[0] = 1.294849661688697e-01;
		GaussWeights[1] = 2.797053914892767e-01;
		GaussWeights[2] = 3.818300505051189e-01;
		GaussWeights[3] = 4.179591836734694e-01;
		GaussWeights[4] = 3.818300505051189e-01;
		GaussWeights[5] = 2.797053914892767e-01;
		GaussWeights[6] = 1.294849661688697e-01;
		break;
	case  8:
		GaussWeights[0] = 1.012285362903763e-01;
		GaussWeights[1] = 2.223810344533745e-01;
		GaussWeights[2] = 3.137066458778873e-01;
		GaussWeights[3] = 3.626837833783620e-01;
		GaussWeights[4] = 3.626837833783620e-01;
		GaussWeights[5] = 3.137066458778873e-01;
		GaussWeights[6] = 2.223810344533745e-01;
		GaussWeights[7] = 1.012285362903763e-01;
		break;
	case  9:
		GaussWeights[0] = 8.127438836157441e-02;
		GaussWeights[1] = 1.806481606948574e-01;
		GaussWeights[2] = 2.606106964029355e-01;
		GaussWeights[3] = 3.123470770400028e-01;
		GaussWeights[4] = 3.302393550012598e-01;
		GaussWeights[5] = 3.123470770400028e-01;
		GaussWeights[6] = 2.606106964029355e-01;
		GaussWeights[7] = 1.806481606948574e-01;
		GaussWeights[8] = 8.127438836157441e-02;
		break;
	case  10:
		GaussWeights[0] = 6.667134430868814e-02;
		GaussWeights[1] = 1.494513491505806e-01;
		GaussWeights[2] = 2.190863625159820e-01;
		GaussWeights[3] = 2.692667193099964e-01;
		GaussWeights[4] = 2.955242247147529e-01;
		GaussWeights[5] = 2.955242247147529e-01;
		GaussWeights[6] = 2.692667193099964e-01;
		GaussWeights[7] = 2.190863625159820e-01;
		GaussWeights[8] = 1.494513491505806e-01;
		GaussWeights[9] = 6.667134430868814e-02;
		break;
	case 11:
		GaussWeights[0] = 5.566856711617367e-02;
		GaussWeights[1] = 1.255803694649046e-01;
		GaussWeights[2] = 1.862902109277343e-01;
		GaussWeights[3] = 2.331937645919905e-01;
		GaussWeights[4] = 2.628045445102467e-01;
		GaussWeights[5] = 2.729250867779006e-01;
		GaussWeights[6] = 2.628045445102467e-01;
		GaussWeights[7] = 2.331937645919905e-01;
		GaussWeights[8] = 1.862902109277343e-01;
		GaussWeights[9] = 1.255803694649046e-01;
		GaussWeights[10] = 5.566856711617367e-02;
		break;
	case 12:
		GaussWeights[0] = 4.717533638651183e-02;
		GaussWeights[1] = 1.069393259953184e-01;
		GaussWeights[2] = 1.600783285433462e-01;
		GaussWeights[3] = 2.031674267230659e-01;
		GaussWeights[4] = 2.334925365383548e-01;
		GaussWeights[5] = 2.491470458134028e-01;
		GaussWeights[6] = 2.491470458134028e-01;
		GaussWeights[7] = 2.334925365383548e-01;
		GaussWeights[8] = 2.031674267230659e-01;
		GaussWeights[9] = 1.600783285433462e-01;
		GaussWeights[10] = 1.069393259953184e-01;
		GaussWeights[11] = 4.717533638651183e-02;
		break;
	case 14:
		GaussWeights[0] = 3.511946033175186e-02;
		GaussWeights[1] = 8.015808715976021e-02;
		GaussWeights[2] = 1.215185706879032e-01;
		GaussWeights[3] = 1.572031671581935e-01;
		GaussWeights[4] = 1.855383974779378e-01;
		GaussWeights[5] = 2.051984637212956e-01;
		GaussWeights[6] = 2.152638534631578e-01;
		GaussWeights[7] = 2.152638534631578e-01;
		GaussWeights[8] = 2.051984637212956e-01;
		GaussWeights[9] = 1.855383974779378e-01;
		GaussWeights[10] = 1.572031671581935e-01;
		GaussWeights[11] = 1.215185706879032e-01;
		GaussWeights[12] = 8.015808715976021e-02;
		GaussWeights[13] = 3.511946033175186e-02;
		break;
	case 16:
		GaussWeights[0] = 2.715245941175409e-02;
		GaussWeights[1] = 6.225352393864789e-02;
		GaussWeights[2] = 9.515851168249278e-02;
		GaussWeights[3] = 1.246289712555339e-01;
		GaussWeights[4] = 1.495959888165767e-01;
		GaussWeights[5] = 1.691565193950025e-01;
		GaussWeights[6] = 1.826034150449236e-01;
		GaussWeights[7] = 1.894506104550685e-01;
		GaussWeights[8] = 1.894506104550685e-01;
		GaussWeights[9] = 1.826034150449236e-01;
		GaussWeights[10] = 1.691565193950025e-01;
		GaussWeights[11] = 1.495959888165767e-01;
		GaussWeights[12] = 1.246289712555339e-01;
		GaussWeights[13] = 9.515851168249278e-02;
		GaussWeights[14] = 6.225352393864789e-02;
		GaussWeights[15] = 2.715245941175409e-02;
		break;
	}

	//return GaussWeights;
}

