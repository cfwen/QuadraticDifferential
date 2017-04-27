#include "QuadraticDifferential.h"
#include "parser/parser.h"

// fix ARPACK names on win32
#ifdef _WIN32
#define dsaupd_ DSAUPD
#define dseupd_ DSEUPD
#endif
#include "unsupported/Eigen/ArpackSupport"


CQuadraticDifferential::CQuadraticDifferential()
{
    mesh = new CMesh();
    geometry = Unknown;
}


CQuadraticDifferential::~CQuadraticDifferential()
{
}

int CQuadraticDifferential::read(const string & filename, const set<string> & traits)
{
    mesh->read_m(filename, traits);

    // read u
    for (CVertex * v : mesh->vertices())
    {
        MeshLib::CParser parser(v->string());
        vector<MeshLib::CToken*> & tokens = parser.tokens();
        double u = 0.0;
        for (auto pt : tokens)
        {
            if (pt->key() == "u")
            {
                string str = strutil::trim(pt->value(), "()");
                istringstream iss(str);
                iss >> u;                
            }
        }
        v->u() = u;
    }
    return 0;
}

int CQuadraticDifferential::write(const string & filename, const set<string> & traits)
{
    mesh->write_m(filename, traits);
    return 0;
}

int CQuadraticDifferential::compute()
{
    if (geometry == Unknown)
    {
        int ec = calculate_euler_characteristic();
        if (ec < 0) geometry = Hyperbolic;
        else if (ec == 0) geometry = Euclidean;
        else
        {
            geometry = Spherical;
            cout << "Spherical geometry not supported" << endl;
            exit(-1);
        }
    }
    int k = 0;
    for (CVertex * v : mesh->vertices())
    {
        v->index() = k++;
    }
    k = 0;
    for (CFace * f : mesh->faces())
    {
        f->index() = k++;
    }
    k = 0;
    for (CEdge * e : mesh->edges())
    {
        e->index() = k++;
        if (geometry == Euclidean)
        {
            e->length() = exp(e->vertex1()->u()) + exp(e->vertex2()->u());
        }
        else if (geometry == Hyperbolic)
        {
            e->length() = -log(tanh(-e->vertex1()->u() / 2)) - log(tanh(-e->vertex2()->u() / 2));
        }
    }

    for (CVertex * v : mesh->vertices())
    {
        embedLocal(v);
    }

    int nv = mesh->num_vertices();
    int ne = mesh->num_edges();
    int kd = ne - 3 * nv;
    if (kd == 0) kd = 2;
    Eigen::SparseMatrix<double> A(3 * nv, ne);
    typedef Eigen::Triplet<double> T;
    vector<T> triplets;
    triplets.reserve(ne * 6);

    for (CVertex * v : mesh->vertices())
    {
        int vid = v->index();
        Nei & nei = v->nei();
        for (CVertex * vj : v->vertices())
        {
            CEdge * e = mesh->edge(v,vj);
            int eid = e->index();
            int vjid = vj->index();
            triplets.push_back(T(vid, eid, 1));
            CPoint2 dz = nei[vjid];
            triplets.push_back(T(vid + nv, eid, dz[0] / dz.norm2()));
            triplets.push_back(T(vid + nv * 2, eid, -dz[1] / dz.norm2()));
        }
    }
    ofstream ofs("A");
    for (auto t : triplets)
    {
        ofs << t.row() << " " << t.col() << " " << t.value() << endl;
    }
    ofs.close();

    A.setFromTriplets(triplets.begin(), triplets.end());
    A.finalize();

    Eigen::SparseMatrix<double> AA = A.transpose()*A;
    Eigen::ArpackGeneralizedSelfAdjointEigenSolver<Eigen::SparseMatrix<double>, Eigen::SparseLU<Eigen::SparseMatrix<double>>> eigs(AA, kd, "SM");
    int nconv = eigs.getNbrConvergedEigenValues();

    cout << "there should be " << kd << " zero singular values" << endl;
    cout << "number of singular values converged is " << nconv << endl;
    cout << eigs.eigenvalues() << endl;

    if (nconv > 0)
    {
        Eigen::VectorXd evalues = eigs.eigenvalues();
        basis = eigs.eigenvectors();
        for (CEdge * e : mesh->edges())
        {
            int id = e->index();
            Eigen::RowVectorXd ev = basis.row(id);
            e->q() = vector<double>(ev.data(), ev.data() + ev.size());
            ostringstream oss;
            oss << "q=(" << ev << ")";
            e->string() = oss.str();
        }
    }

    return 0;
}

inline void hyper_circle_to_circle(complex<double> & c, double & r)
{
    double d = abs(c)*abs(c);
    double f = (exp(r) - 1.0) / (exp(r) + 1.0);
    f = f*f;
    double a = 1.0 - f*d;
    complex<double> b = 2.0*(f - 1.0)*c / a;
    d = (d - f) / a;
    c = -b / 2.0;
    r = sqrt(abs(b)*abs(b) / 4.0 - d);
}
/*
 * for a triangle, given uv of two vertices, and edge length with third vertex,
 * compute uv of third vertex
 */
inline CPoint2 CQuadraticDifferential::embedVertex3(CPoint2 & p1, CPoint2 & p2, double l13, double l23)
{
    complex<double> c1(p1[0], p1[1]);
    complex<double> c2(p2[0], p2[1]);
    double r1 = l13;
    double r2 = l23;
    if (geometry == Hyperbolic)
    {
        hyper_circle_to_circle(c1, r1);
        hyper_circle_to_circle(c2, r2);
    }
    complex<double> dz = c2 - c1;
    double d = abs(dz);
    if (d > r1 + r2 || d < fabs(r2 - r1))
    {
        cerr << "triangle can't be embedded" << endl;
        exit(-1);
    }

    double a = (r1*r1 - r2*r2 + d*d) / d / 2;
    complex<double> z0 = c1 + dz*a / d;
    double h = sqrt(r1*r1 - a*a);
    complex<double> rz = (dz*h / d)*complex<double>(0, 1);
    complex<double> p3 = z0 + rz;
    complex<double> p4 = z0 - rz;
    complex<double> e1 = c2 - c1;
    complex<double> e2 = p3 - c1;

    bool b = real(e1)*imag(e2) - imag(e1)*real(e2) > 0;
    if (b) return CPoint2(real(p3), imag(p3));
    else return CPoint2(real(p4), imag(p4));
}

/*
 * locally embed v and its neighbors to plain/Poincare disk, always assume v is embedded to origin
 */
int CQuadraticDifferential::embedLocal(CVertex * v)
{
    Nei & nei = v->nei();
    // embed first vertex
    CHalfEdge * he0 = v->most_clw_out_halfedge();
    CVertex * v1 = he0->target();
    CPoint2 p0(0,0);
    double l01 = he0->edge()->length();
    CPoint2 p1(l01, 0);
    if (geometry == Hyperbolic) p1[0] = (exp(l01) - 1.0) / (exp(l01) + 1.0);
    nei[v1->index()] = p1;
    CHalfEdge * he1 = he0;
    CHalfEdge * he2 = he0->ccw_rotate_about_source();

    while (he2 && he2 != he0)
    {
        CVertex * v2 = he2->target();
        double l02 = he2->edge()->length();
        double l12 = he1->next()->edge()->length();
        p1 = embedVertex3(p0, p1, l02, l12);
        nei[v2->index()] = p1;
        he1 = he2;
        he2 = he2->ccw_rotate_about_source();
    }
    
    return 0;
}

