#include "gcmethod_1d.h"

int main()
{
	mesh_1d m("gcmethod_1d.ini");
	gcmethod_1d g(m, "gcmethod_1d.ini");
    //g.calculate();
    g.calculate_diverge_contour();
	return 0;
}

