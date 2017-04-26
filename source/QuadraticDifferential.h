#ifndef _QUADRATIC_DIFFERENTIAL_H_
#define _QUADRATIC_DIFFERENTIAL_H_

#include <iostream>
#include <string>
#include <set>
#include <map>
#include "Mesh/mesh.h"
#include "Eigen/Eigen"
//#include "armadillo"

using namespace std;

#ifndef ADD_PROPERTY
#define ADD_PROPERTY(T, x) \
private:\
    T m_##x; \
public:\
    inline T & x() { return m_##x; }
#endif // ! ADD_PROPERTY

#ifndef  PI
#define PI 3.14159265358979323846
#endif // ! PI

enum Geometry {Unknown, Spherical, Euclidean, Hyperbolic };
using CPoint2 = MeshLib::CPoint2;
using CPoint = MeshLib::CPoint;

class CQuadraticDifferential
{
public:
    using Nei = map<int, CPoint2>;
    class CQDVertex
    {
        ADD_PROPERTY(int, index)
        ADD_PROPERTY(Nei, nei)
        ADD_PROPERTY(double, u)
    };
    class CQDEdge
    {
        ADD_PROPERTY(int, index)
        ADD_PROPERTY(double, length)
        ADD_PROPERTY(vector<double>, q)
    };
    class CQDFace
    {
        ADD_PROPERTY(int, index)
    };
    class CQDHalfEdge
    {

    };

    using CMesh = MeshLib::CBaseMesh<CQDVertex, CQDEdge, CQDFace, CQDHalfEdge>;
    using CVertex = typename CMesh::CVertex;
    using CEdge = typename CMesh::CEdge;
    using CFace = typename CMesh::CFace;
    using CHalfEdge = typename CMesh::CHalfEdge;

public:
    CQuadraticDifferential();
    ~CQuadraticDifferential();

    int read(const string & filename, const set<string> & traits = {});
    int write(const string & filename, const set<string> & traits = {});

    int calculate_euler_characteristic()
    {
        return mesh->num_vertices() - mesh->num_edges() + mesh->num_faces();
    }
    void setGeometry(Geometry g) { geometry = g; }
    int compute();

protected:
    int embedLocal(CVertex * v);

private:
    CMesh * mesh;
    Eigen::MatrixXd basis;
    Geometry geometry;
};

#endif // !_QUADRATIC_DIFFERENTIAL_H_
