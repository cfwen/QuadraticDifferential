// main.cpp : Defines the entry point for the console application.
//

#include "QuadraticDifferential.h"

int main(int argc, char * argv[])
{
    if (argc < 3)
    {
        cout << "usage: QuadraticDifferential input-mfile output-mfile" << endl;
        cout << "         - input-mesh: mfile with conformal factor u on vertex" << endl;
        cout << "         - output-mesh: mfile with quadratic differential q on edge" << endl;
        return 0;
    }

    CQuadraticDifferential qd;
    qd.read(argv[1], { "Vertex:u" });
    //qd.setGeometry(Geometry::Hyperbolic);
    qd.compute();
    qd.write(argv[2]);
    
    return 0;
}

